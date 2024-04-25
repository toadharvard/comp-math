for image in ./images/*; do
    for c in 1 3; do
        for m in pi numpy rsvd pcafast; do
            svdi compress --in-file="$image" --out-file="/tmp/tmp.svdi" --compression=$c --method="$m" --power-iterations=10 && svdi decompress --in-file="/tmp/tmp.svdi" --out-file="./compressed/${c}.${m}.${image##*/}"
        done
    done
done
