"""
Parse and plot the matplotlib stuff
"""

import os
from typing import Tuple

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr

RCP_SCENARIOS = [26, 45, 85]
MODELS = ["DP16", "K14"]


def parse_file(rcp: float, model: str) -> xr.DataArray:
    """
    Parse a single file for a particular RCP scenario
    """

    fname = f"LSLproj_MC_{model}_DEL_24_rcp{rcp}.tsv"
    assert os.path.isfile(fname), f"{fname} is not a valid file"

    df = pd.read_csv(fname, delim_whitespace=True, skiprows=2, header=None).set_index(0)
    ds = xr.DataArray(
        df.values[:, :, np.newaxis, np.newaxis],
        coords={
            "time": df.index.values,
            "sample": df.columns,
            "scenario": [rcp],
            "model": [model],
        },
        dims=["time", "sample", "scenario", "model"],
    )
    ds.name = "lsl"
    return ds


def parse_files() -> xr.Dataset:
    """
    Parse all three RCP scenarios
    """
    fname = "scenarios-decadal.nc"
    try:
        ds = xr.open_dataset(fname)
    except FileNotFoundError:
        ds = xr.concat(
            [
                xr.concat(
                    [parse_file(rcp=rcp, model=model) for rcp in RCP_SCENARIOS],
                    dim="scenario",
                )
                for model in MODELS
            ],
            dim="model",
        ).interp(time=np.arange(2010, 2201), method="cubic")
        ds.to_netcdf(fname, format="NETCDF4")
        ds = xr.open_dataset(fname)
    return ds


def plot_scenarios(n_samples: int) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot some samples of the scenarios
    """

    scenarios = parse_files()
    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(RCP_SCENARIOS) + 1,
        figsize=(12, 3.5),
        sharey=True,
        gridspec_kw={"width_ratios": [1, 1, 1, 0.5]},
    )

    for rcp, ax in zip(RCP_SCENARIOS, axes):

        samples_keep = np.random.choice(
            scenarios["sample"].values, size=n_samples, replace=False
        )

        for model, color in zip(MODELS, ["gray", "blue"]):

            for i in samples_keep:

                scenarios["lsl"].sel(
                    sample=i, scenario=rcp, time=slice(2000, 2200), model=model,
                ).plot(
                    ax=ax, c=color, linewidth=0.1, alpha=0.25,
                )

        ax.set_title("RCP Scenario {}".format(rcp / 10.0))
        ax.set_ylim(top=1000)
        ax.set_ylabel(None)
        ax.set_xlabel(None)

    for rcp in RCP_SCENARIOS:
        df = scenarios["lsl"].sel(time=2200, scenario=rcp).to_pandas()
        sns.distplot(
            df, hist=False, rug=False, ax=axes[-1], vertical=True, label=rcp / 10.0
        )

    axes[-1].set_xlabel(None)
    axes[-1].set_xticks([])
    axes[-1].legend()
    axes[-1].set_title("LSL in 2200")

    lines = [
        Line2D([0], [0], color=c, linewidth=1, linestyle="-")
        for c in ["gray", "blue"]
    ]
    axes[0].legend(lines, MODELS)
    axes[0].set_ylabel("Local Sea Level [mm]")

    fig.tight_layout()
    return fig, axes


if __name__ == "__main__":

    _ = plot_scenarios(n_samples=500)
    plt.savefig("simulations.pdf")
