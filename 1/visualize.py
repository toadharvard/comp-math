import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact

x = BLOCK_SIZES = [
    1,  # 1 [0: 11]
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    16,
    32,
    64,
    128,  # 2 [11: 21]
    256,
    300,
    400,
    500,
    600,
    700,
    800,
    900,
    1000,
    2000,
    2730,  # 3 [21:]
    3000,
    4000,
    5000,
    6000,
    7000,
    8000,
    9000,
    10000,
    11000,
    12000,
    13000,
    14000,
]
a = np.load('per_thread.npy')

@interact
def plot_graph():
    for (_, y, _) in a:
        plt.plot(x, y)
    plt.gca().set_aspect('auto')
    plt.gca().autoscale(enable=True, axis='both', tight=False)
    plt.gca().autoscale_view()
    plt.show()