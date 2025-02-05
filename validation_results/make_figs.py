import os
import glob
import numpy as np
import matplotlib.pyplot as plt

directory = "data"
files = glob.glob(os.path.join(directory, "*"))
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple",
          "tab:brown", "tab:pink", "tab:gray", "tab:olive", "tab:cyan",
          "purple", "teal", "gold", "m", "darkolivegreen", "saddlebrown",
          "firebrick", "dodgerblue","indigo", "darkslategrey",
          "aqua", "tomato", "darksalmon", "tan", "crimson"]
for fname in files:
    f = open(fname, "r")
    header = f.readline().split(";")
    data = np.loadtxt(f, delimiter=";")
    try:
        time = data[0, :]
    except IndexError:
        continue
    fig, ax = plt.subplots(1, 1)
    for i, level in enumerate(header[1:]):
        ax.plot(time, data[i+1, :], colors[i], label="%s mV"%level)
    if "Ca_con" not in fname:
        ax.set_ylabel("Current (nA)")
    else:
        ax.set_ylabel("Ca concentration (mM)")
    ax.set_xlabel("time (ms)")
    ax.set_title(fname)
    ax.legend()
    print(fname[:-4]+".png")
    fig.savefig(fname[:-4]+".png", dpi=100)
    plt.show()
