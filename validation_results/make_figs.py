import os
import glob
import numpy as np
import matplotlib.pyplot as plt

directory = "data"
files = glob.glob(os.path.join(directory, "*"))

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
        ax.plot(time, data[i+1, :], label="%s mV"%level)
    if "calcium" not in fname:
        ax.set_ylabel("Current (nA)")
    else:
        ax.set_ylabel("Ca concentration (mM)")
    ax.set_xlabel("time (ms)")
    ax.set_title(fname)
    ax.legend()
plt.show()
