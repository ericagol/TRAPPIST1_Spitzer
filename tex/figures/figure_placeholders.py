import matplotlib.pyplot as plt
import glob
import os

figures = list(map(os.path.basename, glob.glob("static/*")))

for figure in figures:
    if figure.endswith(".pdf") or figure.endswith(".png"):
        fig, ax = plt.subplots(1)
        ax.axis("off")
        ax.annotate(
            "Figure failed to compile.",
            xy=(0.5, 0.5),
            xycoords="data",
            va="center",
            ha="center",
            fontsize=30,
            color="red",
        )
        fig.savefig(figure, bbox_inches="tight")
        plt.close()
