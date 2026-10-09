from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .partition import Results, smooth  # noqa: E402

Y_LABEL = r"Partition ratio  $C_{lipid} / C_{bulk}$"


def _style(ax: plt.Axes) -> None:
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_timeseries(results: Results, window: int, path: str | Path, dpi: int = 300) -> None:
    drugs = sorted({k[1] for k in results})
    ratios = sorted({k[0] for k in results})
    colors = dict(zip(ratios, plt.cm.tab10(np.linspace(0, 1, max(len(ratios), 1)))))

    fig, axes = plt.subplots(1, len(drugs), figsize=(6 * len(drugs), 4.5), sharey=True, squeeze=False)
    for ax, drug in zip(axes[0], drugs):
        for ratio in ratios:
            labelled = False
            for (r, d, _), df in sorted(results.items()):
                if (r, d) != (ratio, drug):
                    continue
                ax.plot(
                    df["time_ps"] / 1000.0,
                    smooth(df["partition_ratio"], window),
                    color=colors[ratio],
                    alpha=0.85,
                    lw=1.8,
                    label=None if labelled else ratio,
                )
                labelled = True
        _style(ax)
        ax.set_xlabel("Time (ns)")
        ax.set_title(drug.upper())
        ax.legend(title="Composition", fontsize=8)
    axes[0][0].set_ylabel(Y_LABEL)
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)


def plot_bars(grouped: pd.DataFrame, tail_fraction: float, path: str | Path, dpi: int = 300) -> None:
    ratios = sorted(grouped["ratio"].unique())
    drugs = sorted(grouped["drug"].unique())
    x = np.arange(len(ratios))
    width = 0.8 / len(drugs)

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, drug in enumerate(drugs):
        subset = grouped[grouped["drug"] == drug].set_index("ratio").reindex(ratios)
        means = subset["mean"].to_numpy()
        errors = subset["std"].fillna(0.0).to_numpy()
        ax.bar(x + i * width, means, width=width, label=drug.upper(), edgecolor="black", linewidth=1.5)
        ax.errorbar(x + i * width, means, yerr=errors, fmt="none", ecolor="k", capsize=4, elinewidth=1.5)
    _style(ax)
    ax.set_xticks(x + width * (len(drugs) - 1) / 2)
    ax.set_xticklabels(ratios)
    ax.set_xlabel("Lipid composition")
    ax.set_ylabel(Y_LABEL + f"  (last {tail_fraction:.0%} avg)")
    ax.legend(title="Drug")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
