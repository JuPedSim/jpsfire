import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def sort_files(_files):
    """
    sort files wrt.  time
    """
    tmp_files = []
    for _file in _files:
        tmp_files.append(_file.split("t_")[-1].split(".npz")[0])

    return sorted(tmp_files, key=lambda x: int(os.path.splitext(x)[0]))


def plot_npz_dir(_dirname):
    doors = os.path.join(_dirname, "Door*")
    print(doors)
    dirs = glob.glob(doors)
    print(dirs)
    time_sorted_files = []
    for d in dirs:
        files = glob.glob("%s/*.npz"%d)
        time_sorted_files = sort_files(files)
        print(d)
        for i, f in enumerate(time_sorted_files):
            sf = os.path.join(d, "t_%s.npz" % f)
            data = np.load(sf)['smoke_factor_grid_norm']
            plt.imshow(data, origin="lower",
                       interpolation="spline36",
                       vmin=0, vmax=10,
                       cmap=mpl.cm.spectral)
            plt.colorbar()
            Time = f.split(".npz")[0].split("_")[-1]
            plt.title("t = %s"%Time)
            print("Time %15s >> %s/%.3d.png"%(Time, d, i))
            plt.savefig("%s/%.3d.png"%(d, i))
            plt.clf()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit("usage: %s <dirname>/3_sfgrids/EXTINCTION_COEFFICIENT/Z_2.250000"%
                 sys.argv[0])

    dirname = sys.argv[1]
    plot_npz_dir(dirname)
