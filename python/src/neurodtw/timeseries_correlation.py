# %% *************  CORRELATE TWO TIME SERIES WITH DYNAMIC TIME WARPING  ***************

# %% Imports
import pdb
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from tslearn.metrics import dtw_path
import scipy.stats as stats

# %% Global constants
TITLE_FONTSIZE = 20
SUPLABEL_FONTSIZE = 15
INFO_TXT_FILENAME = "CorrelateAnalysis_Info.txt"


# %% Main function
def timeseries_correlation(
    series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims
):
    """
    Conceptual Approach
    -------------------
    To assess latency correlations between 2 time series:
      1) Get Single Subject (or group) time series
      2) Get N bootstrap time series (w. same bootstrapping
          subject-distributions for the two time series)
      3) DTW of 2) with 1) using AREA measure!
      4) Gives you two numbers for each bootstrap sample -> DTW result of
          time series 1 & 2 latency (lag)
      5) Correlate these with Pearson and Spearman (latter preferred for
          presentations)

    Input Parameters
    ----------------
    1) series_1, series_2
       => Timeseries as two matrices of shape:
          dataindices (e.g., timepoints) x subjects (or repetitions)
    2) name_1, name_2
       => Names of your timeseries
    3) savepath
       => Where to save output to
    4) num_boots
       => Number of bootstraps (we recommend at least 10000)
    5) outlier
       => 1: to exclude +-5SD outliers from your scatterplot
       => 0: don't exclude outliers
    6) try_to_fix_ylims
       => 1: try to fix ylimits of marginals and include lines for easier
             comparison
             - based on num_boots
             - if you want to modify this yourself (e.g., hardcode lines &
               maxima, change the if statement blocks for the two marginals)
       => 0: let matlab handle ylimits of marginals (likely won't allow
             comparisons

    Note
    ----
    1) Number of subjects has to be equal in both time series!
    2) Your full time series are analysed.
       If you would like to assess only a specific interval of your time
       series, use indexing before calling this function.
    3) We fix marginals' x-axis ranges to +-8 SD. This should be okay for
       most data-sets since we zscore DTW area distributions. However, if you
       do not want to standardise like this (which will make your marginals
       harder to compare) or you want to implement different values, change
       the lines that use xlimits & ylimits (see line 180)
    """

    # %% Preparation

    # first check if user tried to compare groups with unequal sample sizes
    if np.shape(series_1)[1] != np.shape(series_2)[1]:
        print("Your time series have a different number of subjects - cancelling!")
        return

    # prepare some vars
    name_1 = str(name_1)
    name_2 = str(name_2)
    plotpath = os.path.join(savepath, "Plots/")
    varpath = os.path.join(savepath, "Variables/")
    if not os.path.exists(plotpath):
        os.makedirs(plotpath)
    if not os.path.exists(varpath):
        os.makedirs(varpath)
    DTW_marginals = {"name_1": [], "name_2": []}
    DTW_stats = dict()

    # %% Run the analysis
    series_1_dtw_area, series_2_dtw_area = run_bootstrap_dtw(
        series_1, series_2, num_boots
    )

    # %% Outlier Rejection
    if outlier:

        # exclude outliers if +-5 SD (change in std_thresh line if wanted)
        series_1_std = np.std(series_1_dtw_area)
        series_2_std = np.std(series_2_dtw_area)
        std_thresh = 5  # Threshold for outlier rejection

        # *****************************************************************
        #                           IMPORTANT
        # This might seem like a strange way of removing outliers
        # ==> However it ensures that we check all 4 directions of the
        #   scatterplot for outliers, save indices as we do so, then just
        #   take those unique indices (removing duplicates) and (IMPORTANTLY)
        #   IN THE END removing these indices FROM BOTH DTW distributions
        # ==> THIS IS IMPORTANT BECAUSE DTW DISTRIBUTIONS HAVE TO STAY
        #   LINKED DUE TO THE BOOTSTRAP
        # *****************************************************************

        outlier_idx = []
        # series 1 - positive outliers
        series_1_pos_outlier_idx = np.where(
            series_1_dtw_area
            > (np.mean(series_1_dtw_area) + (std_thresh * series_1_std))
        )[0]
        if len(series_1_pos_outlier_idx) > 0:
            outlier_idx.extend(series_1_pos_outlier_idx)
        # series 1 - negative outliers
        series_1_neg_outlier_idx = np.where(
            series_1_dtw_area
            < (np.mean(series_1_dtw_area) - (std_thresh * series_1_std))
        )[0]
        if len(series_1_neg_outlier_idx) > 0:
            outlier_idx.extend(series_1_neg_outlier_idx)
        # series 2 - positive outliers
        series_2_pos_outlier_idx = np.where(
            series_2_dtw_area
            > (np.mean(series_2_dtw_area) + (std_thresh * series_2_std))
        )[0]
        if len(series_2_pos_outlier_idx) > 0:
            outlier_idx.extend(series_2_pos_outlier_idx)
        # series 2 - negative outliers
        series_2_neg_outlier_idx = np.where(
            series_2_dtw_area
            < (np.mean(series_2_dtw_area) - (std_thresh * series_2_std))
        )[0]
        if len(series_2_neg_outlier_idx) > 0:
            outlier_idx.extend(series_1_neg_outlier_idx)

        # remove duplicates
        outlier_idx = np.unique(outlier_idx)

        # remove outliers from BOTH series!
        if len(outlier_idx) > 0:
            series_1_dtw_area = np.delete(series_1_dtw_area, outlier_idx)
            series_2_dtw_area = np.delete(series_2_dtw_area, outlier_idx)

    # %% Save var, skew & kurtosis values of DTW Area Distributions
    DTW_stats["series_1_var"] = np.var(series_1_dtw_area)
    DTW_stats["series_1_skew"] = stats.skew(series_1_dtw_area)
    DTW_stats["series_1_kurt"] = stats.kurtosis(series_1_dtw_area)
    DTW_stats["series_2_var"] = np.var(series_2_dtw_area)
    DTW_stats["series_2_skew"] = stats.skew(series_2_dtw_area)
    DTW_stats["series_2_kurt"] = stats.kurtosis(series_2_dtw_area)
    if outlier:
        picklepath = os.path.join(
            varpath, ("%iSD OutlierRejected Marginal Stats.pkl" % (std_thresh))
        )
        with open(picklepath, "wb") as file:
            pickle.dump(DTW_stats, file)
    else:
        picklepath = os.path.join(varpath, ("Marginal Stats.pkl"))
        with open(picklepath, "wb") as file:
            pickle.dump(DTW_stats, file)

    # %% Standardise DTW Area Distributions & correlate
    # So line of linear fit in scatter & marginal distributions reflect correlation
    # values better (since correlations have internal standardisation)

    # zscore
    z_x = stats.zscore(series_1_dtw_area)
    z_y = stats.zscore(series_2_dtw_area)

    # pearson corr
    pears_res = stats.pearsonr(z_x, z_y)
    pears_r = pears_res[0]
    pears_p = pears_res[1]

    # spearman corr
    spear_res = stats.spearmanr(z_x, z_y)
    spear_r = spear_res[0]
    spear_p = spear_res[1]

    # prepare strings
    str_df = str(num_boots - 2)  # degrees of freedom
    pears_str_p = ""
    spear_str_p = ""
    if pears_p < 0.0001:
        pears_str_p = "p < .0001"
    else:
        pears_str_p = "p = " + str(pears_p.round(3))
    if spear_p < 0.0001:
        spear_str_p = "p < .0001"
    else:
        spear_str_p = "p = " + str(spear_p.round(3))

    # %% Plot one large figure
    f, ax_global = plt.subplots()
    ax_global.set_xticks([])
    ax_global.set_yticks([])
    for spine in ax_global.spines.values():  # despine outmost axis
        spine.set_visible(False)
    gs = f.add_gridspec(4, 4, wspace=0, hspace=0)  # prepare the grid
    scatter_color = "#c79fef"  # lavender
    line_color = "#8ab8fe"  # carolina blue
    xlimits = [-8, 8]
    xticks = np.squeeze([-8, -4, 0, 4, 8])
    ylimits = [-8, 8]
    yticks = xticks.copy()
    if try_to_fix_ylims:
        ylimmax = num_boots / 8
        ref_val_1 = round(ylimmax * 0.25)
        ref_val_2 = round(ylimmax * 0.5)
        ref_val_3 = round(ylimmax * 0.75)

    # %% first, scatterplot
    ax = f.add_subplot(gs[0:3, 1:])  # set wanted gridcells as ax
    ax.scatter(
        z_x, z_y, c="white", edgecolor=scatter_color, s=7, linewidth=0.5, alpha=0.6
    )
    ax.set_ylim(ylimits)
    ax.set_xlim(xlimits)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    series_1_xlimits = ax.get_xlim()  # for plotting marginals later
    series_1_xticks = ax.get_xticks()
    series_2_xlimits = ax.get_ylim()
    series_2_xticks = ax.get_yticks()
    for spine in ax.spines.values():
        spine.set_visible(False)
    # ax.set_title("DTW area-diff correlations of %s & %s.\nPearson: r(%s) = %s: %s & Spearman: r(%s) = %s: %s" % (name_1, name_2, str_df, pears_str_r, pears_str_p, str_df, spear_str_r, spear_str_p))
    ax.set_title(
        "DTW area-diff correlations of "
        + name_1
        + " & "
        + name_2
        + "\nPearson: r("
        + str(str_df)
        + ") = "
        + str(pears_r.round(2))
        + ": "
        + pears_str_p
        + "\nSpearman: r("
        + str(str_df)
        + ") = "
        + str(spear_r.round(2))
        + ": "
        + spear_str_p
    )
    # have to create the line of least squares manually
    coefficients = np.polyfit(z_x, z_y, 1)  # 1 == linear (fit coefficients to orig z_x)
    slope, intercept = coefficients
    # create y_fit using new z_x that spans over the full xlimits (so line isnt cut-off)
    z_x_new = np.arange(xlimits[0], xlimits[1] + 1, 1)
    y_fit = intercept + (slope * z_x_new)
    ax.plot(z_x_new, y_fit, linewidth=1.5, color=line_color)

    # %% second, series_2 marginal distribution
    ax_series_2 = f.add_subplot(gs[:-1, 0])
    ax_series_2.hist(
        z_y,
        100,
        orientation="horizontal",
        facecolor=scatter_color,
        edgecolor=line_color,
        linewidth=0.15,
    )
    ax_series_2.set_ylim(series_2_xlimits)  # set ylim because horizontal hist!
    ax_series_2.set_yticks(series_2_xticks)
    ax_series_2.invert_xaxis()
    ax_series_2.set_ylabel(name_2 + " Marginal Distribution")
    for side in ax_series_2.spines.keys():  # side is left/right/top/bottom
        ax_series_2.spines[side].set_linewidth(0.2)
        if side == "top":
            ax_series_2.spines[side].set_visible(False)
    if try_to_fix_ylims:
        ax_series_2.set_xlim(left=ylimmax)
        ref_val_strings = []
        for val in [ref_val_1, ref_val_2, ref_val_3]:
            ax_series_2.axvline(val, color="#029386", linewidth=0.2, linestyle="--")
            ref_val_strings.append(str(val))
        ax_series_2.set_xticks(
            [ref_val_1, ref_val_2, ref_val_3], labels=ref_val_strings
        )

    # %% third, series_1 marginal distribution
    ax_series_1 = f.add_subplot(gs[-1, 1:])
    ax_series_1.hist(
        z_x, 100, facecolor=scatter_color, edgecolor=line_color, linewidth=0.15
    )
    ax_series_1.set_xlim(series_1_xlimits)
    ax_series_1.set_xticks(series_1_xticks)
    ax_series_1.set_xlabel(name_1 + " Marginal Distribution")
    for side in ax_series_1.spines.keys():  # side is left/right/top/bottom
        ax_series_1.spines[side].set_linewidth(0.2)
        if side == "left":
            ax_series_1.spines[side].set_visible(False)
    if try_to_fix_ylims:
        ax_series_1.set_ylim(top=ylimmax)
        ref_val_strings = []
        for val in [ref_val_1, ref_val_2, ref_val_3]:
            ax_series_1.axhline(val, color="#029386", linewidth=0.2, linestyle="--")
            ref_val_strings.append(str(val))
        ax_series_1.set_yticks(
            [ref_val_1, ref_val_2, ref_val_3], labels=ref_val_strings
        )
    ax_series_1.yaxis.set_ticks_position("right")

    # %% save figures & marginal distributions
    # figures
    if outlier:
        fig_filename = (
            "DTW Latency Correlation - "
            + name_1
            + " & "
            + name_2
            + " - "
            + str(std_thresh)
            + "std outlier removed"
        )
        marginals_filename = (
            "Marginal Distributions - " + str(std_thresh) + "std outlier removed"
        )
    else:
        fig_filename = "DTW Latency Correlation - " + name_1 + " & " + name_2
        marginals_filename = "Marginal Distributions"
    f.savefig(
        os.path.join(plotpath, fig_filename + ".png"),
        bbox_inches="tight",
        dpi=300,
    )
    f.savefig(
        os.path.join(plotpath, fig_filename + ".svg"),
        bbox_inches="tight",
        dpi=300,
    )
    # marginal distributions
    series_2_marginals = series_2_dtw_area
    series_1_marginals = series_1_dtw_area
    z_series_2_marginals = z_y
    z_series_1_marginals = z_x
    vars_to_save = {
        "series_1_marginals": series_1_marginals,
        "series_2_marginals": series_2_marginals,
        "z_series_1_marginals": z_series_1_marginals,
        "z_series_2_marginals": z_series_2_marginals,
    }

    picklepath = os.path.join(varpath, (marginals_filename + ".pkl"))
    with open(picklepath, "wb") as file:
        pickle.dump(vars_to_save, file)
    plt.show()


# %% Local functions
def run_bootstrap_dtw(series_1, series_2, num_boots):
    """Run the main bootstrap DTW analysis"""
    # compute standardised (subjects / reps) averages of both time series
    avg_series_1 = series_1.mean(axis=1)
    avg_series_2 = series_2.mean(axis=1)
    z_avg_series_1 = stats.zscore(avg_series_1)
    z_avg_series_2 = stats.zscore(avg_series_2)

    # bootstrap manually using randint & run dtw between bootstrap & observed z_avg
    num_subs = np.shape(series_1)[1]
    series_1_dtw_area = np.zeros(num_boots)
    series_2_dtw_area = np.zeros(num_boots)

    # bootstrap loop
    for i in range(num_boots):
        print("\nboots num #%s" % (i))

        # get num_subj bootstrap samples (using randint == we bootstrap with duplicates)
        this_trap = np.random.randint(0, num_subs, num_subs)
        boot_series_1 = series_1[:, this_trap]
        boot_series_2 = series_2[:, this_trap]

        # average & standardise trapped series
        avg_boot_series_1 = boot_series_1.mean(axis=1)
        avg_boot_series_2 = boot_series_2.mean(axis=1)
        z_avg_boot_series_1 = stats.zscore(avg_boot_series_1)
        z_avg_boot_series_2 = stats.zscore(avg_boot_series_2)

        # ***** dtw between z-scored averages of series & bootstrapped series *****

        # series_1
        WP_1, score_1 = dtw_path(z_avg_boot_series_1, z_avg_series_1)
        ix_1 = []
        iy_1 = []
        for WP_idxs in WP_1:
            ix_1.append(WP_idxs[0])
            iy_1.append(WP_idxs[1])
        xmax_1 = ix_1[-1] + 1  # because this is used with range
        DIAG_1 = np.arange(xmax_1)
        areaDIAG_1 = np.trapz(DIAG_1)
        # IMPORTANT NOTE ABOUT NP.TRAPZ - TAKES IN (Y) & (Y,X) WHEREAS MATLAB'S TRAPZ
        # TAKES (Y) & (X,Y) !!!!!
        areaWP_1 = np.trapz(iy_1, ix_1)
        series_1_dtw_area[i] = (areaDIAG_1 - areaWP_1) / areaDIAG_1

        # series_2
        WP_2, score_2 = dtw_path(z_avg_boot_series_2, z_avg_series_2)
        ix_2 = []
        iy_2 = []
        for WP_idxs in WP_2:
            ix_2.append(WP_idxs[0])
            iy_2.append(WP_idxs[1])
        xmax_2 = ix_2[-1] + 1  # because this is used with range
        DIAG_2 = np.arange(xmax_2)
        areaDIAG_2 = np.trapz(DIAG_2)
        # IMPORTANT NOTE ABOUT NP.TRAPZ - TAKES IN (Y) & (Y,X) WHEREAS MATLAB'S TRAPZ
        # TAKES (Y) & (X,Y) !!!!!
        areaWP_2 = np.trapz(iy_2, ix_2)
        series_2_dtw_area[i] = (areaDIAG_2 - areaWP_2) / areaDIAG_2

    # bootstrap loop end
    return series_1_dtw_area, series_2_dtw_area
