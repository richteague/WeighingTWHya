"""Scripts to read the pickled fits."""

import warnings
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve_fft
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore")


def get_stats(samples, idx1='pa_deg', idx2='inc_deg', kernel=2, bins=None):
    """
    Get the convolved stats.

    Args:
        samples (dict): Dictionary of samples.
        idx1 (Optional[str]): Name of the value for the x-axis.
        idx2 (Optional[str]): Name of the value for the y-axis.
        kernel (Optinal[int]): Size of the Gaussian kernel to smooth the data.
        bins (Optional[list of ndarrays]): A list of the x- and y-axis bins to
            use.

    Returns:
        x (ndarray): Value at the center of the x-axis bins.
        y (ndarray): Value at the center of the y-axis bins.
        H (ndarray): Smoothed 2D histogram of the samples.
        Hx (ndarray): Collapsed histogram along the y-axis.
        Hy (ndarray): Collapsed histogram along the x-axis.
        V (ndarray): The contours of the [16, 50, 84]the percentiles.
    """

    # Default binning.
    if bins is None:
        bins = [np.linspace(samples[idx1].min(), samples[idx1].max(), 101),
                np.linspace(samples[idx2].min(), samples[idx2].max(), 101)]

    # Calculate the histogram.
    H, xe, ye = np.histogram2d(samples[idx1], samples[idx2],
                               bins=bins, normed=True)
    x = np.average([xe[1:], xe[:-1]], axis=0)
    y = np.average([ye[1:], ye[:-1]], axis=0)
    if kernel > 0:
        H = convolve_fft(H, Gaussian2DKernel(x_stddev=kernel))

    # Make the marginalized values.
    Hx = np.sum(H, axis=1)
    Hx /= np.trapz(Hx, x)
    Hy = np.sum(H, axis=0)
    Hy /= np.trapz(Hy, y)

    # Calculating the proper contour levels for a PDF. Stripped from DFM's
    # corner.py. Note: may hang up with incorrect indexing.
    levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5)**2)
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]
    V = np.empty(len(levels))
    for i, v0 in enumerate(levels):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]
    V.sort()
    m = np.diff(V) == 0
    while np.any(m) and np.sum(H) > 0.0:
        V[np.where(m)[0][0]] *= 1.0 - 1e-4
        m = np.diff(V) == 0
    V.sort()

    return x, y, H.T, Hx, Hy, V


def plot_PDF(samples, idx1='pa_deg', idx2='inc_deg', kernel=2, cmaps=None,
             bins=None, labels=None, nbins=300, spread=0.5, return_fig=False,
             xlabel='Position Angle (deg)', ylabel='Inclination (deg)'):
    """
    Plot the posteriors of the provided runs.

    Args:
        samples (list): List of the samples to plot.
        idx1 (Optional[str]): Name of the value for the x-axis.
        idx2 (Optional[str]): Name of the value for the y-axis.
        kernel (Optional[int]): Size of the Gaussian kernel to smooth the data.
        cmaps (Optional[list]): List of the colormaps to use for the samples.
        bins (Optional[list]): A list of the x- and y-axis bins to use.
        labels (Optional[list]): List of strings containing the legend labels.
        nbins (Optional[int]): Number of automatic bins to use if bins is None.
        spread (Optional[float]): Spread of the x- and y-axis relative to the
            PDF ragenes.
        return_fig (Optional[bool]): Whether to return the axes.
        xlabel(Optional[str]): Label of the x-axis.
        ylabel(Optional[str]): Label of the y-axis.

    Returns:
        ax (matplotlib axes): Matplotlib axis if requested.
    """

    # Make the axes.
    gridspec_kw = dict(height_ratios=[1, 4], width_ratios=[4, 1],
                       wspace=0.04, hspace=0.04)
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(3.75, 3.75),
                            gridspec_kw=gridspec_kw)
    axs[0, 1].axis('off')

    # Default colormaps.
    if cmaps is None:
        cmaps = [cm.Reds, cm.Oranges, cm.Greys, cm.Purples, cm.Blues]
        if not len(samples) % 2:
            idx = int((len(cmaps) - 1) / 2)
            cmaps = cmaps[:idx] + cmaps[-idx:]
        idxs = np.linspace(0, len(cmaps) - 1, len(samples)).astype(int)
        cmaps = np.take(cmaps, idxs)

    # Binning for the data.
    if bins is None:
        x_tmp = np.array([[s[idx1].min(), s[idx1].max()] for s in samples])
        x_tmp = [x_tmp[:, 0].min(), x_tmp[:, 1].max()]
        x_tmp = [x_tmp[0] - spread * (x_tmp[1] - x_tmp[0]),
                 x_tmp[1] + spread * (x_tmp[1] - x_tmp[0])]
        xbins = np.linspace(x_tmp[0], x_tmp[1], int(nbins+1))
        y_tmp = np.array([[s[idx2].min(), s[idx2].max()] for s in samples])
        y_tmp = [y_tmp[:, 0].min(), y_tmp[:, 1].max()]
        y_tmp = [y_tmp[0] - spread * (y_tmp[1] - y_tmp[0]),
                 y_tmp[1] + spread * (y_tmp[1] - y_tmp[0])]
        ybins = np.linspace(y_tmp[0], y_tmp[1], int(nbins+1))
        bins = [xbins, ybins]

    # Plot the samples.
    for s, cmap in zip(samples, cmaps):
        x, y, H, Hx, Hy, V = get_stats(s, bins=bins, idx1=idx1, idx2=idx2,
                                       kernel=kernel)
        for v, vv in enumerate(V):
            axs[1, 0].contour(x, y, H, vv, linewidths=1.0,
                              colors=[cmap(((v + 0.5) / (len(V))))],
                              zorder=v-len(V))
        axs[0, 0].plot(x, Hx, color=cmap(0.7), lw=1.0)
        axs[0, 0].fill_between(x, Hx, color=cmap(0.7), alpha=0.3, lw=0.0,
                               zorder=-1)
        axs[1, 1].plot(Hy, y, color=cmap(0.7), lw=1.0)
        axs[1, 1].fill_betweenx(y, Hy, color=cmap(0.7), alpha=0.3, lw=0.0,
                                zorder=-1)

    # Set the limits.
    axs[1, 0].set_xlim(bins[0][0], bins[0][-1])
    axs[1, 0].set_ylim(bins[1][0], bins[1][-1])
    axs[0, 0].set_xlim(bins[0][0], bins[0][-1])
    axs[1, 1].set_ylim(bins[1][0], bins[1][-1])
    axs[0, 0].set_ylim(0, axs[0, 0].get_ylim()[1])
    axs[1, 1].set_xlim(0, axs[1, 1].get_xlim()[1])

    # Remove the spines from the marginalized PDFs.
    axs[1, 0].tick_params(which='both', right=1, top=1)
    axs[0, 0].set_xticklabels([])
    axs[1, 1].set_yticklabels([])
    axs[0, 0].spines['top'].set_visible(False)
    axs[0, 0].spines['left'].set_visible(False)
    axs[0, 0].spines['right'].set_visible(False)
    axs[0, 0].set_yticks([])
    axs[1, 1].spines['top'].set_visible(False)
    axs[1, 1].spines['bottom'].set_visible(False)
    axs[1, 1].spines['right'].set_visible(False)
    axs[1, 1].set_xticks([])

    # Plot labels.
    if labels is not None:
        for l, cmap in zip(labels, cmaps):
            axs[0, 1].plot([np.nan], [np.nan], color=cmap(0.7), label=l)
        axs[0, 1].legend(markerfirst=True, fontsize=6)
    axs[1, 0].set_xlabel(xlabel)
    axs[1, 0].set_ylabel(ylabel)

    if return_fig:
        return fig


def combined_PDF(samples, idx, as_error=True):
    """Get the percentiles of the combined PDFs."""
    combined = samples[0][idx]
    for s in samples[1:]:
        combined = np.concatenate([combined, s[idx]])
    pcnts = np.percentile(combined, [16, 50, 84])
    if not as_error:
        return pcnts
    return pcnts[1], pcnts[1] - pcnts[0], pcnts[2] - pcnts[1]
