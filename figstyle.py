"""
Shared figure styling for the cross-species study.

Matches the Lancet-style aesthetic used in the MIMIC-IV HIV-sepsis figures:
serif type, soft grey panels with white gridlines, grey spines, filled area
under ringed (white-centred) markers, rounded panel-label boxes, and a
crimson-to-navy diverging colormap.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

mpl.use("Agg")

PALETTE = {
    "ink": "#1a1a1a",
    "navy": "#00468B",
    "crimson": "#AD002A",
    "green": "#42B540",
    "sky": "#0099B4",
    "purple": "#925E9F",
    "orange": "#FDAF91",
    "red": "#ED0000",
    "muted": "#ADB6B6",
    "grid": "#E6E9ED",
    "panel": "#F0F0F0",
    "bg": "#ffffff",
    "spine": "#B0B0B0",
}
CAT_COLORS = ["#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6"]

CMAP_SEQ = LinearSegmentedColormap.from_list(
    "lancet_seq", ["#ffffff", "#D6E4EF", "#8FB9D6", "#0099B4", "#00468B", "#01254A"]
)
CMAP_DIV = LinearSegmentedColormap.from_list("lancet_div", ["#AD002A", "#f4f4f4", "#00468B"])

_available = {f.name for f in fm.fontManager.ttflist}
SERIF = "Times New Roman" if "Times New Roman" in _available else "DejaVu Serif"


def apply_style():
    mpl.rcParams.update(
        {
            "figure.dpi": 130,
            "savefig.dpi": 300,
            "figure.facecolor": PALETTE["bg"],
            "savefig.facecolor": PALETTE["bg"],
            "savefig.bbox": "tight",
            "font.family": SERIF,
            "font.serif": [SERIF, "Times New Roman", "Times", "DejaVu Serif", "STIXGeneral"],
            "mathtext.fontset": "stix",
            "font.size": 11,
            "axes.titlesize": 13.5,
            "axes.titleweight": "bold",
            "axes.titlepad": 12,
            "axes.labelsize": 11.5,
            "axes.labelcolor": PALETTE["ink"],
            "axes.edgecolor": PALETTE["ink"],
            "axes.linewidth": 0.9,
            "axes.facecolor": PALETTE["bg"],
            "axes.grid": False,
            "grid.color": PALETTE["grid"],
            "grid.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "xtick.color": PALETTE["ink"],
            "ytick.color": PALETTE["ink"],
            "text.color": PALETTE["ink"],
            "legend.frameon": False,
            "legend.fontsize": 10,
            "axes.unicode_minus": False,
        }
    )


def style_panel(ax, grid_axis="both"):
    """Soft grey panel with white gridlines and grey spines."""
    ax.set_facecolor(PALETTE["panel"])
    ax.grid(True, axis=grid_axis, color="#FFFFFF", lw=1.2, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color(PALETTE["spine"])
        spine.set_linewidth(1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10)


def panel_label(ax, letter, title, x=0.0, y=1.02, inside=False):
    """Panel label. By default placed just above the axes (outside the plot box)
    so it never occludes data; set inside=True for an in-panel boxed label."""
    if inside:
        ax.text(
            0.02, 0.98, f"({letter})  {title}", transform=ax.transAxes,
            ha="left", va="top", fontsize=12, fontweight="bold", color=PALETTE["ink"],
            bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#D0D5DD",
                      linewidth=0.8, alpha=0.95), zorder=6,
        )
    else:
        ax.text(
            x, y, f"({letter})  {title}", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=12, fontweight="bold", color=PALETTE["ink"],
            zorder=6,
        )


def ringed(ax, x, y, color, s=70, lw=2.2, zorder=5):
    ax.scatter(x, y, s=s, facecolors="white", edgecolors=color, linewidths=lw, zorder=zorder)


def area_line(ax, x, y, color, lw=2.8, y_floor=None):
    if y_floor is None:
        y_floor = min(y)
    ax.fill_between(x, y_floor, y, alpha=0.30, color=color, zorder=1)
    ax.plot(x, y, "-", color=color, lw=lw, zorder=3, solid_capstyle="round")


def save_fig(fig, name, figdir: Path, formats=("png", "pdf", "eps", "tif")):
    for ax in fig.get_axes():
        items = [ax.title, ax.xaxis.label, ax.yaxis.label]
        items += ax.get_xticklabels() + ax.get_yticklabels()
        for item in items:
            item.set_fontfamily(SERIF)
        leg = ax.get_legend()
        if leg is not None:
            for t in leg.get_texts():
                t.set_fontfamily(SERIF)
    png_path = figdir / f"{name}.png"
    if "png" in formats:
        fig.savefig(png_path)
    if "pdf" in formats:
        fig.savefig(figdir / f"{name}.pdf")
    if "eps" in formats:
        # EPS does not support transparency; flatten to white for journal upload.
        fig.savefig(
            figdir / f"{name}.eps",
            format="eps",
            facecolor="white",
            edgecolor="none",
            transparent=False,
        )
    plt.close(fig)
    if "tif" in formats and png_path.exists():
        try:
            from PIL import Image as PILImage

            im = PILImage.open(png_path).convert("RGB")
            im.save(figdir / f"{name}.tif", format="TIFF", compression="tiff_lzw", dpi=(300, 300))
        except Exception as e:  # pragma: no cover
            print(f"  (tif skipped for {name}: {e})")
