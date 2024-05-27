# %% ************  LATENCY DIFFERENCES WITH DYNAMIC TIME WARPING (IN MS)  **************

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


# %% Main function
def latency_difference(
    analysis_design,
    query,
    reference,
    name_query,
    name_reference,
    units,
    sampling_rate,
    filepath,
    permutations,
):
    """
    Conceptual Approach
    -------------------
    To assess a latency difference between 2 time series:
      1) Average across subjects & compute DTW between the averages
      2) The median of the distribution of the warping path's x/y distances
          is used to compute the latency difference in milliseconds
      3) The area between the warping path and the main diagonal is used
          in a permutation procedure for assigning a p value to the latency
          difference

    Note
    ----
    1) This function allows assessment of latency differences in within-
        or between-subjects designs. Set the "analysis_design" input parameter
        accordingly. This affects permutation procedure, adopting the
        permutation equivalents of paired- or independent-sample t-tests
        respectively.
    2) This whole function is in line with the proposition of
        Zoumpoulaki et al. (2015)
        Latency as a region contrast: Measuring ERP latency differences with
        Dynamic Time Warping
    3) Your full time series are analysed.
        If you would like to assess only a specific interval of your time
        series, use indexing before calling this function.

    Input Parameters
    ----------------
    1) analysis_design
        => Character array or string informing about your analysis design
        => Has to be either "within" or "between", based on whether you have
          a within- or between-subjects design
    2) query
        => Query time series (matrix w. dimensions: datapoints x subjects)
    3) reference
        => Reference time series (matrix w. dimensions: datapoints x subjects)
    4) name_query
        => Name of the query time series
    5) name_reference
        => Name of the reference time series
    6) units
        => In what units was your data measured?
          - For microvolts enter '\muV'
    7) sampling_rate
        => What was the sampling rate of your signal? In Hertz
        => Used to compute the latency difference in milliseconds
    8) filepath
        => Path with files
    9) permutations
        => Number of permutations to adopt in statistical testing
          - We recommend at least 10000
    """

    # %% Preparation

    analysis_design = str(analysis_design)
    units = str(units)
    name_query = str(name_query)
    name_reference = str(name_reference)
    filepath = str(filepath)

    # Sanity check on analysis_design value
    if analysis_design not in ["between", "within"]:
        print('analysis_design must be either "between" or "within"! Fix & re-run')
        return

    # Sanity check if series are identical for within and have same length
    try:
        query = np.array(query)
        reference = np.array(reference)
    except:
        raise ValueError(
            "Query & Ref input vars must be convertible to np arrays - fix & rerun!"
        )
    # => can have different number of subjects if between
    if analysis_design == "within":
        if query.shape != reference.shape:
            print(
                "Note - query & reference don't have identical shapes - fix & re-run!"
            )
            return
    else:  # has to be 'between' now (
        if query.shape[0] != reference.shape[0]:
            print(
                "Note - query & reference do not have same datapoint-length - fix & "
                + "re-run!"
            )
            return

    # %% MAIN DTW COMPUTATION
    # We need indices of x & y coordinates of warping path!
    #   ORIGINAL FORM: WP, score = dtw(x,y) (we don't need score)
    #   ==> QUERY IS ON X-AXIS OF WP, SO FIRST INPUT
    #   ==> REFERENCE IS ON Y-AXIS OF WP, SO SECOND INPUT!
    #   ==> WP is a list of lists with [ix, iy] coordinates of WP
    #   ==> WP of this will be below main diagonal if QUERY is
    #       LATER than REFERENCE
    #   ==> I.e. AREA will be positive in this case because area = areaDIAG -
    #       areaWP (and hence areaDIAG will be larger than areaWP if WP is below
    #       DIAG [and QUERY is indeed LATER than REFERENCE!]

    # average across subjects
    plot_avg_query = query.mean(axis=1)
    plot_avg_reference = reference.mean(axis=1)

    # DTW has to be done on standardised averages
    z_avg_query = stats.zscore(plot_avg_query)
    z_avg_reference = stats.zscore(plot_avg_reference)
    WP, score = dtw_path(z_avg_query, z_avg_reference)

    # extract indices of WP's x & y coordinates
    i_query_x = []
    i_reference_y = []
    for WP_idxs in WP:
        i_query_x.append(WP_idxs[0])
        i_reference_y.append(WP_idxs[1])
    warped_query = z_avg_query[i_query_x]  # apply indices to get aligned time series
    warped_reference = z_avg_reference[i_reference_y]

    # plot standardised & aligned time series
    plot_standardised_and_aligned_series(
        z_avg_query, z_avg_reference, warped_query, warped_reference, units, filepath
    )

    # Latency in MS is distribution of point-wise distances:
    # (distance = ix - iy) / sampling rate (Hz) / 1000
    distance = np.array(i_query_x) - np.array(i_reference_y)
    latency = np.median(distance)
    lat_in_ms = latency / (sampling_rate / 1000)

    print("\n *** Latency Difference: %i ms! *** \n" % (lat_in_ms))

    # Plot WP compared to main diagonal
    plotWP(
        i_query_x,
        i_reference_y,
        plot_avg_query,
        plot_avg_reference,
        name_query,
        name_reference,
        units,
        lat_in_ms,
        filepath,
    )

    # %% PERMUTATION TEST ON AREA MEASURE
    # ==> have this for median(distance), i.e., lat_in_ms_measure, too

    # compute true observed area
    xmax = i_query_x[-1] + 1  # because this is used with range
    DIAG = np.arange(xmax)
    areaDIAG = np.trapz(DIAG)
    # IMPORTANT NOTE ABOUT NP.TRAPZ - TAKES IN (Y) & (Y,X) WHEREAS MATLAB'S TRAPZ
    # TAKES (Y) & (X,Y) !!!!!
    areaWP = np.trapz(i_reference_y, i_query_x)
    area = (areaDIAG - areaWP) / areaDIAG

    # prepare vars for permuting
    original_query = query.copy()  # copy so it's not just a reference
    original_reference = reference.copy()
    perm_x = [[] for _ in range(permutations)]
    perm_y = [[] for _ in range(permutations)]
    perm_area = [[] for _ in range(permutations)]
    datapoints = np.shape(query)[0]
    subjects = np.shape(query)[1]

    # permutation loop
    for perms in range(permutations):
        # for a within-subject contrast, do the coin flipping for permuting class labels
        if analysis_design == "within":
            perm_query = np.zeros((datapoints, subjects))
            perm_reference = np.zeros((datapoints, subjects))
            for s in range(subjects):
                coin = np.random.randint(1, 3, subjects)  # 3 bc. indexing exclusive
                if coin[s] == 1:
                    perm_query[:, s] = original_query[:, s]
                    perm_reference[:, s] = original_reference[:, s]
                elif coin[s] == 2:
                    perm_query[:, s] = original_reference[:, s]
                    perm_reference[:, s] = original_query[:, s]
        # for a between-subject contrast, randomly assign subjects to groups, preserving
        # the original group sizes of query & reference
        else:
            subjects_query = np.shape(query)[1]
            subjects_reference = np.shape(reference)[1]
            subjects_total = subjects_query + subjects_reference
            this_perm = np.random.permutation(subjects_total)
            data_total = np.concatenate((query, reference), axis=1)
            perm_query = data_total[:, this_perm[:subjects_query]]
            perm_reference = data_total[:, this_perm[subjects_query:]]

        # get averages across subjects
        perm_query_avg = perm_query.mean(axis=1)
        perm_reference_avg = perm_reference.mean(axis=1)

        # DTW has to be done on standardised averages
        perm_query_z_avg = stats.zscore(perm_query_avg)
        perm_reference_z_avg = stats.zscore(perm_reference_avg)
        perm_WP, perm_score = dtw_path(perm_query_z_avg, perm_reference_z_avg)

        # extract indices of permuted WP's x & y coordinates - store in lists
        perm_i_query_x = []
        perm_i_reference_y = []
        for WP_idxs in perm_WP:
            perm_i_query_x.append(WP_idxs[0])
            perm_i_reference_y.append(WP_idxs[1])
        perm_x[perms] = perm_i_query_x
        perm_y[perms] = perm_i_reference_y

        # compute permuted area - store results in list
        perm_xmax = perm_i_query_x[-1] + 1  # because this is used with range
        perm_DIAG = np.arange(perm_xmax)
        perm_areaDIAG = np.trapz(perm_DIAG)
        # IMPORTANT NOTE ABOUT NP.TRAPZ - TAKES IN (Y) & (Y,X) WHEREAS MATLAB'S TRAPZ
        # TAKES (Y) & (X,Y) !!!!!
        perm_areaWP = np.trapz(perm_i_reference_y, perm_i_query_x)
        perm_area[perms] = (perm_areaDIAG - perm_areaWP) / perm_areaDIAG

    # permutation loop end
    permdist_areas = np.sort(perm_area)  # sort nonabsolutes for idxing thresholds
    sort_abs_permdist_areas = np.sort(abs(permdist_areas))
    thresh_abs = round(len(permdist_areas) * 0.95)
    # find 5% sig. thresh using sorted absolutes
    thresh_abs_val = sort_abs_permdist_areas[thresh_abs]
    thresh1 = -thresh_abs_val
    thresh2 = thresh_abs_val

    # p value (two tailed test)
    p = sum(abs(permdist_areas) >= abs(area)) / permutations
    if p < 0.0001:
        plot_p = "p < 0.0001"
    else:
        plot_p = "p = " + str(round(p, 3))
    print("\n *** %s *** \n" % (plot_p))

    # Plot figure of observed area-size & its permutation distribution
    f, ax = plt.subplots()
    plt.tight_layout()
    ax.hist(permdist_areas, bins=30, facecolor="#c7fdb5", edgecolor="#e4cbff")
    ax.axvline(thresh1, linewidth=1.5, linestyle="-", color="#ff796c")
    ax.axvline(thresh2, linewidth=1.5, linestyle="-", color="#ff796c")
    ax.axvline(area, linewidth=1.5, linestyle="-", color="#8ab8fe")
    yticks = ax.get_yticks()
    yticklabels_num = yticks / permutations
    yticklabels_str = [str(tick) for tick in yticklabels_num]
    ax.set_yticks(yticks, labels=yticklabels_str)
    ax.set_ylabel("Proportion", fontsize=12)
    ax.set_xlabel("DTW - Area", fontsize=12)
    ax.set_title(
        "Latency diff.: %i ms - %s" % (lat_in_ms, plot_p),
        fontweight="bold",
        fontsize=13,
    )
    f.savefig(
        os.path.join(filepath, "Permutation Distribution.png"),
        bbox_inches="tight",
        dpi=300,
    )
    f.savefig(
        os.path.join(filepath, "Permutation Distribution.svg"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()

    # Plot figure of observed area-size & its permutation distribution
    f, ax = plt.subplots()
    plt.tight_layout()
    ax.hist(abs(permdist_areas), bins=30, facecolor="#c7fdb5", edgecolor="#e4cbff")
    ax.axvline(thresh_abs_val, linewidth=1.5, linestyle="-", color="#ff796c")
    ax.axvline(abs(area), linewidth=1.5, linestyle="-", color="#8ab8fe")
    yticks = ax.get_yticks()
    yticklabels_num = yticks / permutations
    yticklabels_str = [str(tick) for tick in yticklabels_num]
    ax.set_yticks(yticks, labels=yticklabels_str)
    ax.set_ylabel("Proportion", fontsize=12)
    ax.set_xlabel("DTW - Area (absolute)", fontsize=12)
    ax.set_title(
        "Absolute Areas - Latency diff.: %i ms - %s" % (lat_in_ms, plot_p),
        fontweight="bold",
        fontsize=13,
    )
    f.savefig(
        os.path.join(filepath, "Permutation Distribution (absolutes).png"),
        bbox_inches="tight",
        dpi=300,
    )
    f.savefig(
        os.path.join(filepath, "Permutation Distribution (absolutes).svg"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()

    # store variables to a pickle file
    vars_to_save = {
        "area": area,
        "p": p,
        "thresh1": thresh1,
        "thresh2": thresh2,
        "perm_x": perm_x,
        "perm_y": perm_y,
        "permdist_areas": permdist_areas,
    }
    picklepath = os.path.join(
        filepath, ("DTW Permutation Test - %s subjects.pkl" % analysis_design)
    )
    with open(picklepath, "wb") as file:
        pickle.dump(vars_to_save, file)
    print(
        "\n\n****************************"
        + "\nDTW Latency Analysis is done"
        + "\n****************************\n\nFigures & vars were saved to filepath."
    )


# %% Local functions - plotting
def plot_standardised_and_aligned_series(
    z_avg_query, z_avg_reference, warped_query, warped_reference, units, filepath
):
    """Plots a figure of original/standardised and aligned time series"""
    f, ax = plt.subplots(4, 1)
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.9)
    ax[0].plot(z_avg_query, color="k")  # plot the data
    ax[1].plot(z_avg_reference, color="r")
    ax[2].plot(warped_query, color="k")
    ax[3].plot(warped_reference, color="r")
    # ylims of standardised time series
    stand_min = min([z_avg_query.min(), z_avg_reference.min(), -1])
    stand_max = max([z_avg_query.max(), z_avg_reference.max(), 1])
    ax[0].set_ylim([stand_min, stand_max])
    ax[1].set_ylim([stand_min, stand_max])
    # ylims of warped time series
    warped_min = min([warped_query.min(), warped_reference.min(), -1])
    warped_max = max([warped_query.max(), warped_reference.max(), -1])
    ax[2].set_ylim([warped_min, warped_max])
    ax[3].set_ylim([warped_min, warped_max])
    # labels, titles & save
    ax[-1].set_xlabel("Datapoints", fontsize=SUPLABEL_FONTSIZE)
    f.suptitle("Alignment", y=0.99, fontsize=TITLE_FONTSIZE)
    f.supylabel(units + " (z-scored)", x=0.01, fontsize=SUPLABEL_FONTSIZE)
    f.savefig(os.path.join(filepath, "Alignment.png"), bbox_inches="tight", dpi=300)
    f.savefig(os.path.join(filepath, "Alignment.svg"), bbox_inches="tight", dpi=300)
    plt.show()


def plotWP(
    i_query_x,
    i_reference_y,
    plot_avg_query,
    plot_avg_reference,
    name_query,
    name_reference,
    units,
    lat_in_ms,
    filepath,
):
    """Plot the warping path, the main diagonal and both original time series"""
    f, ax_global = plt.subplots()
    ax_global.set_xticks([])
    ax_global.set_yticks([])
    gs = f.add_gridspec(4, 4)  # prepare the grid
    plt.tight_layout()
    plt.subplots_adjust(left=0.1, bottom=0.1, top=0.9, hspace=0.5)

    # *** grid space 1: WP & main diagonal ***
    # main diagonal
    ax = f.add_subplot(gs[0:3, 1:4])  # set wanted gridcells as ax
    ax.set_xlim(right=i_query_x[-1])
    ax.set_ylim(top=i_reference_y[-1])
    xmax = i_query_x[-1] + 1  # bc. indexing starts at 0 & we use this var only w. range
    ax.plot(
        np.arange(xmax), np.arange(xmax), color="#ff796c", linewidth=1.5
    )  # c==salmon
    # warping path
    ax.plot(i_query_x, i_reference_y, color="#8ab8fe", linewidth=1.5)  # c==carolinablue
    # fill area between WP & main diagonal
    ax.fill_between(
        i_query_x, i_query_x, i_reference_y, color="#c7fdb5", alpha=0.5, edgecolor="k"
    )
    # remove ticklabels
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # *** grid space 2: query on x axis ***
    ax_2 = f.add_subplot(gs[-1, 1:])
    ax_2.plot(plot_avg_query, color="r")

    # *** grid space 3: reference on y axis ***
    ax_3 = f.add_subplot(gs[0:-1, 0])
    # IMPORTANT NOTE - plt does not have camroll as matlab does. Thus we:
    # 1. plot reference's y values (e.g. in mV) on "x axis"
    # 2. plot len of idxs (xmax) as "y values"
    # 3. we now HAVE TO invert the xaxis to get the same effect as with camroll
    ax_3.plot(plot_avg_reference, np.arange(xmax), color="k")
    ax_3.invert_xaxis()

    # axis stuff
    ylimmax = max([plot_avg_query.max(), plot_avg_reference.max()])
    ylimmax = ylimmax + (ylimmax * 0.1)
    ax_2.set_ylim(top=ylimmax)
    ax_3.set_ylim(top=ylimmax)
    xticks = ax.get_xticks()  # set xaxis of ax_2 and yaxis of ax_3 according to WP axis
    xlims = ax.get_xlim()
    ax_2.set_xticks(xticks)
    ax_2.set_xlim(xlims)
    ax_3.set_yticks(xticks)
    ax_3.set_ylim(xlims)
    # query & reference text & labels
    ax_2.text(
        0.77,
        1.1,
        name_query,
        transform=ax_2.transAxes,
        fontsize=12,
        fontweight="bold",
        fontstyle="italic",
    )
    ax_3.text(
        0.2,
        1.03,
        name_reference,
        transform=ax_3.transAxes,
        fontsize=12,
        fontweight="bold",
        fontstyle="italic",
    )
    ax_2.set_ylabel(units, labelpad=1)
    ax_3.set_xlabel(units, labelpad=1)
    for spine in ax_global.spines.values():  # despine outmost axis
        spine.set_visible(False)

    # suplabels & save the figure
    ax.set_title(
        "Latency difference: %i ms" % lat_in_ms, fontweight="bold", fontsize=13
    )
    ax_2.set_xlabel("Datapoints", fontsize=11)
    ax_3.set_ylabel("Datapoints", fontsize=11)
    f.savefig(os.path.join(filepath, "WarpPath.png"), bbox_inches="tight", dpi=300)
    f.savefig(os.path.join(filepath, "WarpPath.svg"), bbox_inches="tight", dpi=300)
    plt.show()
