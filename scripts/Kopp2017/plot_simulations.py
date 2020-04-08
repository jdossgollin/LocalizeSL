"""
Plot the simulations
"""

import matplotlib.pyplot as plt
import xarray as xr


def plot_scenarios(n_samples: int) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot some samples of the scenarios
    """

    scenarios = xr.open_dataset("scenarios-decadal.nc")
    fig, axes = plt.subplots(
        nrows=1, ncols=len(RCP_SCENARIOS), figsize=(10, 3.5), sharey=True, sharex=True
    )

    for rcp, ax in zip(RCP_SCENARIOS, axes):

        samples_keep = np.random.choice(
            scenarios["sample"].values, size=n_samples, replace=False
        )

        for i in samples_keep:

            scenarios.sel(sample=i, scenario=rcp, time=slice(2000, 2200)).plot(
                ax=ax, c="gray", linewidth=0.1, alpha=0.25,
            )

        ax.set_title("RCP Scenario {}".format(rcp / 10.0))
        ax.set_ylim(top=750)

    fig.tight_layout()
    return fig, axes


if __name__ == "__main__":

    fig, axes = plot_scenarios(n_samples=250)
    plt.savefig("simulations.pdf")
