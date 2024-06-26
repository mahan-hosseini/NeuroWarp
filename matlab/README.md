# Matlab NeuroWarp

Matlab NeuroWarp consists of two main scripts, **neurowarp_timeseries_correlation.m & neurowarp_latency_difference.m** - these scripts enable DTW analyses of the kind presented in *Transient Attention Gates Access Consciousness: Coupling N2pc and P3 Latencies using Dynamic Time Warping*. We additionally provide the Event Related Potentials presented in our paper in the **example_series_N2pcP3s.mat** file for users to familiarise themselves with our scripts before using them with their own data. Please see below for  further instructions.

The scripts and data used to generate our analyses and figures as presented in our paper can be found in the **replicating_paper** folder.  

## Dependencies
*Matlab NeuroWarp requires the following toolboxes to be installed and available in your Matlab path:*
- Matlab’s Statistics & Machine Learning Toolbox
- Matlab’s Signal Processing Toolbox
- [Hline & vline](https://de.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline )
- [RGB](https://de.mathworks.com/matlabcentral/fileexchange/46872-intuitive-rgb-color-values-from-xkcd)
- [Suplabel](https://de.mathworks.com/matlabcentral/fileexchange/7772-suplabel )

## DTW Temporal Correlation - neurowarp_timeseries_correlation.m
Perform a **DTW based bootstrap analysis** to assess the **temporal correlation between two time series** (Figure 5 of our paper).

### Important Notes
- Time series must be **2D matrices**
	- I.e., data points (e.g. time) x subjects (i.e., replications)
- We provide the ERPs of correct and intrusion trials for users to explore this function

### Running neurowarp_timeseries_correlation.m using our ERPs:
*Enter the following into MATLAB’s command window. Note that example_series_N2pcP3.mat has to be stored in your current directory (see step #1)*
1.	`load("example_series_N2pcP3s.mat")`
2.	`series_1 = P3_Correct`
3.	`series_2 = N2pc_Correct`
4.	`name_1 = “P3”`
5.	`name_2 = “N2pc”`
6.	`savepath = cd` (or wherever else you would like to store results to)
7.	`num_boots = 10000`
	- The number of bootstrap samples that you want to implement
8.	`outlier = 0`
	- Exclude outliers if their DTW area is +-5 standard deviations from the mean
9.	`try_to_fix_ylims = 1`
	- Attempt to standardise y-limits of marginals
10.	`neurowarp_timeseries_correlation(series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims)`

*Note that the figure will look slightly different to that of our paper due to different x/y limits. See the replicate_figures folder if you want to replicate our figure as it was printed.*

## DTW Latency Difference - neurowarp_latency_difference.m
Assess the **latency difference** between **two conditions** (i.e., within-subjects effect) or between **two groups** (i.e., across-subjects effect) of any signal of interest (in milliseconds).

*Figures 3 & 4 of our paper show a between-conditions (within-subjects) analysis*

### Important Notes
- Reference and query time series must be **2D matrices**
	- I.e., data points (e.g., time) x subjects (i.e., replications)
	- Time series have to be of **equal sizes**
- **analysis_design** determines whether you want to assess a within- or between-subjects latency effect (can only take “within” or “between” as input)
- We provide the ERPs of correct and intrusion trials for users to explore this function

### Running neurowarp_latency_difference.m using our ERPs:
*Enter the following into MATLAB’s command window. Note that example_series_N2pcP3.mat has to be stored in your current directory (see step #1)*
1.	`load("example_series_N2pcP3s.mat")`
2.	`analysis_design = “within”`
3.	`query = N2pc_Intrusion`
4.	`reference = N2pc_Correct`
5.	`name_query = "Intrusion"`
6.	`name_reference = "Correct"`
7.	`units = "\muV"`
8.	`sampling_rate = 500`
	- The number of data points per second in Hertz
9.	`savepath = cd` (or wherever else you would like to store results to)
10.	`permutations = 10000`
	- The number of permutations you would like to implement in statistical testing (we recommend >=10000)
11.	`neurowarp_latency_difference(analysis_design, query, reference, name_query, name_reference, units, sampling_rate, savepath, permutations)`

## Tests
Matlab NeuroWarp was tested with Matlab 2020b & 2023b on Windows and 2023b on Mac OS.