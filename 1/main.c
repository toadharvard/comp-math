#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))

double APPROXIMATE_ERROR = 0.1;
double ABSOLUTE_ERROR = -1;
size_t BLOCK_SIZE = 1000;
size_t NET_SIZE = 20000;
size_t THREAD_COUNT = 8;
size_t ITERATIONS_COUNT = 1;
double K_PARAM = 1;
bool DEBUG = false;

void set_debug(bool debug)
{
    DEBUG = debug;
}

void set_iterations(size_t iterations_count)
{
    ITERATIONS_COUNT = iterations_count;
}

void set_globals(
    size_t thread_count, double approximate_error, double absolute_error, size_t block_size, size_t net_size, double k_param)
{
    APPROXIMATE_ERROR = approximate_error;
    ABSOLUTE_ERROR = absolute_error;
    BLOCK_SIZE = block_size;
    NET_SIZE = net_size;
    THREAD_COUNT = thread_count;
    K_PARAM = k_param;
}

typedef struct exec_res_t
{
    size_t number_of_iterations;
    double total_time;
} exec_res_t;

typedef double (*fun_xy)(double, double);

typedef struct net_t
{
    size_t size;
    double h;
    double **u;
    double **f;
} net_t;

typedef struct ref_net_t
{
    size_t size;
    double h;
    double **u;
} ref_net_t;

typedef struct u_func_t
{
    fun_xy u;
    fun_xy f;
} u_func_t;

double **create_double_2d_arr(size_t size)
{
    double **res = calloc(size, sizeof(*res));
    for (int i = 0; i < size; i++)
        res[i] = calloc(size, sizeof(*res[i]));
    return res;
}

void free_double_2d_arr(double **arr, size_t size)
{
    for (int i = 0; i < size; i++)
        free(arr[i]);
    return free(arr);
}

net_t *create_net(
    size_t size,
    u_func_t u_func)
{
    net_t *res = malloc(sizeof(*res));
    res->size = size;
    res->h = 1.0 / (size - 1);
    res->u = create_double_2d_arr(size);
    res->f = create_double_2d_arr(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if ((i == 0) || (j == 0) || (i == (size - 1)) || (j == (size - 1)))
            {
                res->u[i][j] = u_func.u(i * res->h, j * res->h);
            }
            else
            {
                res->u[i][j] = 0;
            }
            res->f[i][j] = u_func.f(i * res->h, j * res->h);
        }
    }
    return res;
}

void free_net(net_t *net)
{
    free_double_2d_arr(net->u, net->size);
    free_double_2d_arr(net->f, net->size);
    return free(net);
}

u_func_t create_u_func(fun_xy u, fun_xy f)
{
    u_func_t res;
    res.f = f;
    res.u = u;
    return res;
}

ref_net_t *create_ref_net(size_t size, fun_xy u)
{
    ref_net_t *res = malloc(sizeof(*res));
    double h = 1.0 / (size - 1);
    double **arr = create_double_2d_arr(size);

    res->u = arr;
    res->h = h;
    res->size = size;

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            res->u[i][j] = u(i * h, j * h);
        }
    }
    return res;
}

void free_ref_net(ref_net_t *net)
{
    free_double_2d_arr(net->u, net->size);
    return free(net);
}

double process_block(net_t *net, size_t a, size_t b)
{
    int i0 = 1 + a * BLOCK_SIZE;
    int im = min(i0 + BLOCK_SIZE, net->size - 1);
    int j0 = 1 + b * BLOCK_SIZE;
    int jm = min(j0 + BLOCK_SIZE, net->size - 1);

    double dm = 0;
    for (int i = i0; i < im; i++)
    {
        for (int j = j0; j < jm; j++)
        {
            double temp = net->u[i][j];
            net->u[i][j] = 0.25 * (net->u[i - 1][j] + net->u[i + 1][j] + net->u[i][j - 1] + net->u[i][j + 1] - net->h * net->h * net->f[i][j]);
            double d = fabs(temp - net->u[i][j]);
            dm = max(dm, d);
        }
    }
    return dm;
}

double process_iteration(net_t *net, double *dm, double dmax, size_t number_of_blocks)
{
    dmax = 0;
    for (int nx = 0; nx < number_of_blocks; nx++)
    {
        dm[nx] = 0;

        int i, j;
        double d;

#pragma omp parallel for shared(net, nx, dm) private(i, j, d)
        for (i = 0; i < nx + 1; i++)
        {
            j = nx - i;
            d = process_block(net, i, j);
            dm[i] = max(dm[i], d);
        }
    }

    for (int nx = number_of_blocks - 2; nx > -1; nx--)
    {
        int i, j;
        double d;

#pragma omp parallel for shared(net, nx, dm) private(i, j, d)
        for (i = 0; i < nx + 1; i++)
        {
            j = 2 * (number_of_blocks - 1) - nx - i;
            d = process_block(net, i, j);
            dm[i] = max(dm[i], d);
        }
    }

    for (int i = 0; i < number_of_blocks; i++)
        dmax = max(dm[i], dmax);
    return dmax;
}

bool approximateCond(double dmax)
{
    return dmax > APPROXIMATE_ERROR;
}

bool absoluteCond(net_t *net, u_func_t used_function)
// using mape
{
    ref_net_t *ref = create_ref_net(net->size, used_function.u);
    double sm = 0;
    double mx = -INFINITY;
    double mn = INFINITY;

    for (int i = 0; i < net->size; i++)
    {
        for (int j = 0; j < net->size; j++)
        {
            double diff = (ref->u[i][j] - net->u[i][j]);
            mx = max(mx, diff);
            mn = min(mn, diff);
            sm += fabs(diff);
        }
    }
    free_ref_net(ref);
    return sm / (net->size * net->size) > ABSOLUTE_ERROR;
}

size_t process_net(net_t *net, u_func_t used_function)
{
    size_t number_of_iterations = 0;
    size_t size_without_borders = net->size - 2;
    size_t number_of_blocks = size_without_borders / BLOCK_SIZE;
    number_of_blocks += (BLOCK_SIZE * number_of_blocks != size_without_borders) ? 1 : 0;
    double dmax = 0;
    double *dm = calloc(number_of_blocks, sizeof(*dm));

    if (DEBUG){
        printf("Number of blocks = %zd, blocks/thread = %lf, BLOCK_SIZE = %zd\n", number_of_blocks, number_of_blocks * 1.0 / THREAD_COUNT, BLOCK_SIZE);
    }

    do
    {
        number_of_iterations += 1;
        dmax = process_iteration(net, dm, dmax, number_of_blocks);
        if (ITERATIONS_COUNT != -1 && number_of_iterations % ITERATIONS_COUNT == 0)
        {
            break;
        }
    } while (APPROXIMATE_ERROR != -1 ? approximateCond(dmax) : absoluteCond(net, used_function));
    free(dm);
    return number_of_iterations;
}

// Examples
double u_1(double x, double y) { return 1 * x * x * y * (x - 1) * (y - 1); }
double f_1(double x, double y) { return 1 * (6 * x - 2) * (y - 1) * y + 2 * 1 * (x - 1) * x * x; }

double u_2(double x, double y) { return pow(x, 2.0) + pow(y, 3.0); }
double f_2(double x, double y) { return 6 * y + 2; }

double u_3(double x, double y) { return 7 * pow(x, 2.0) - 10 * pow(x, 3.0) + x + y - 1; }
double f_3(double x, double y) { return 14 - 60 * x; }
double u_4(double x, double y) { return exp(x + 2 * y) + pow(x, 2.0); }
double f_4(double x, double y) { return 5 * exp(x + 2 * y) + 2; }
double u_5(double x, double y) { return 1000 * pow(x, 3) + 2000 * pow(y, 3); }
double f_5(double x, double y) { return 6000 * x + 12000 * y; }

double u_6 /*from book with K param*/ (double x, double y) { return (1 - 2 * y) * (1 - 2 * x) * K_PARAM; }
double f_6 /*from book*/ (double x, double y) { return 0; }

exec_res_t run(u_func_t used_function)
{
    net_t *net = create_net(NET_SIZE, used_function);
    double t1, t2, dt;

    t1 = omp_get_wtime();
    size_t number_of_iterations = process_net(net, used_function);
    t2 = omp_get_wtime();
    dt = t2 - t1;

    free_net(net);
    exec_res_t res = {.number_of_iterations = number_of_iterations, .total_time = dt};
    return res;
}

exec_res_t execute(int function_number)
{
    u_func_t allowed_functions[] = {
        create_u_func(u_1, f_1),
        create_u_func(u_2, f_2),
        create_u_func(u_3, f_3),
        create_u_func(u_4, f_4),
        create_u_func(u_5, f_5),
        create_u_func(u_6, f_6),
    };
    u_func_t used_function;
    if (function_number <= 6 && function_number >= 1)
    {
        used_function = allowed_functions[function_number - 1];
    }
    else
    {
        used_function = allowed_functions[0];
    }
    omp_set_num_threads(THREAD_COUNT);
    exec_res_t res = run(used_function);
    return res;
}

int main(int argc, char *argv[])
{
    exec_res_t res = execute(5);
    printf("Iterations: %zu, Time: %f\n", res.number_of_iterations, res.total_time);
}