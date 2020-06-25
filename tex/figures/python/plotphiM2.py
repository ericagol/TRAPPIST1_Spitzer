import matplotlib

# matplotlib.use('PS')
matplotlib.use("agg")

import matplotlib.image as mpimg
import pylab as pl
import numpy as np

pl.rc("font", size=16)
params = {"legend.fontsize": 13}
pl.rcParams.update(params)
# pl.rc("text", usetex=True) <-- this doesn't work on Travis unless we install latex

pl.figure(figsize=(8, 6))


grid = pl.GridSpec(5, 10, wspace=0.05, hspace=0.5)


for p in range(1, 6, 1):

    if p == 1:
        ax = pl.subplot(grid[0, 0:9])
        ax.text(
            0.05,
            0.95,
            "b,c,d",
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        # legend
        pl.plot(
            ([0, 0]), ([0, 0]), lw=2, c="k", label=r"in resonance $d\phi < 45^\circ$"
        )
        pl.plot(
            ([0, 0]), ([0, 0]), lw=2, c="b", label=r"in resonance $d\phi ≥ 45^\circ$"
        )
        pl.plot(([0, 0]), ([0, 0]), lw=2, c="deepskyblue", label="not in resonance")
        # pl.axis('off')
        ax.legend(
            loc="upper center",
            ncol=3,
            markerscale=50.0,
            bbox_to_anchor=(0.5, 1.8),
            fancybox=True,
            shadow=True,
        )
    if p == 2:
        ax = pl.subplot(grid[1, 0:9])
        ax.text(
            0.05,
            0.95,
            "c,d,e",
            horizontalalignment="left",
            color="w",
            verticalalignment="top",
            transform=ax.transAxes,
        )
    if p == 3:
        ax = pl.subplot(grid[2, 0:9])
        ax.text(
            0.05,
            0.95,
            "d,e,f",
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
        )
    if p == 4:
        ax = pl.subplot(grid[3, 0:9])
        ax.text(
            0.05,
            0.95,
            "e,f,g",
            horizontalalignment="left",
            color="w",
            verticalalignment="top",
            transform=ax.transAxes,
        )
    if p == 5:
        ax = pl.subplot(grid[4, 0:9])
        ax.text(
            0.05,
            0.95,
            "f,g,h",
            horizontalalignment="left",
            color="w",
            verticalalignment="top",
            transform=ax.transAxes,
        )
        pl.xlabel("time [Myr]")

    pl.ylabel(r"$\phi[^{\circ}]$")

    # ax.set_xlim(-0.2, 10.0)

    img = mpimg.imread("../../../data/Grimm/tlM%02d.png" % p)
    imgplot = pl.imshow(img, aspect=0.55)

    h, w, c = img.shape

    print(h, w, c)

    xt = (np.arange(0, 11, 1) + 0.2) * w / 10.2
    xl = np.arange(0, 11, 1)

    yt = np.arange(0, 361, 90) * h / 360
    yl = [360, "", 180, "", 0]

    pl.xticks(xt, xl)
    pl.yticks(yt, yl)

    if p != 5:
        ax.xaxis.set_ticklabels([])

bincenter, histT1, histT2, histT3, histT4, histT5 = np.loadtxt(
    "../../../data/Grimm/hist.dat", unpack=True, usecols=(2, 3, 4, 5, 6, 7)
)

for p in range(1, 6, 1):
    if p == 1:
        ax6 = pl.subplot(grid[0, 9])
        pl.xlim([10, 10000000])
        pl.barh(bincenter, histT1, align="center")
    if p == 2:
        ax6 = pl.subplot(grid[1, 9])
        pl.xlim([10, 10000000])
        pl.barh(bincenter, histT2, align="center")
    if p == 3:
        ax6 = pl.subplot(grid[2, 9])
        pl.xlim([10, 10000000])
        pl.barh(bincenter, histT3, align="center")
    if p == 4:
        ax6 = pl.subplot(grid[3, 9])
        pl.xlim([10, 10000000])
        pl.barh(bincenter, histT4, align="center")
    if p == 5:
        ax6 = pl.subplot(grid[4, 9])
        pl.xlim([10, 10000000])
        pl.barh(bincenter, histT5, align="center")
        pl.xlabel(r"$log_{10}$ counts")

    pl.xscale("log")
    ax6.yaxis.set_ticklabels([])
    ax6.set_ylim(0, 360)
    yl = []
    yt = [360, 270, 180, 90, 0]
    pl.yticks(yt, yl)

    if p != 5:
        ax6.xaxis.set_ticklabels([])


name = "../tlM.png"
pl.savefig(name, dpi=300)
pl.clf()
