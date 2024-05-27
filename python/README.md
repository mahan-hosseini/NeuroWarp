# Python NeuroDTW

Python NeuroDTW is provided as a module that can be installed via PyPi. It consists of two functions that can be called after importing neurodtw: timeseries_correlation and latency_difference -  these are the general-purpose scripts enable the DTW analyses presented in *Transient Attention Gates Access Consciousness: Coupling N2pc and P3 Latencies using Dynamic Time Warping*.

## Installation

We recommend that you create a new virtual environment for our module via:
1. Open a terminal / cmd.exe window, navigate (via `cd`) to the directory you want to create the environment and enter:
`python -m venv neurodtw_env`
2. Activate the neurodtw_env via:

	**Windows:** `path\to\neurodtw_env\Scripts\activate`

	**MacOS:** `source neurodtw_env/bin/activate`

3. Install neurodtw via pip (enter while neurodtw_env is active):
`pip install neurodtw`
4. The functions can be used as explained below after importing the module

For more detail on virtual environments & pip [click here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)


## DTW Temporal Correlation - neurodtw.timeseries_correlation()
Perform a DTW based bootstrap analysis to assess the temporal correlation between two time series (Figure 5 of our paper).

Time series must be 2D matrices
I.e., data points (e.g. time) x subjects (i.e., replications)
We provide the ERPs of correct and intrusion trials for users to explore this function


**Running timeseries_correlation using our ERPs**
__Make sure to enter your actual paths!__
1. Run Python, import the neurodtw package and load our ERPs
`import neurodtw`	
`from scipy.io import loadmat`
`data = loadmat("your/path/to/example_series_N2pcP3s")`

**Configure Variables**
2. series_1 = data["P3_Correct"]
3. series_2 = data["P3_Intrusion"]
4. name_1 = "Correct"
5. name_2 = "Intrusion"
6. savepath = "where/to/store/results/to"
7. num_boots = 10000
	- The number of bootstrap samples that you want to implement
8. outlier = 0
	- Exclude outliers if their DTW area is +-5 standard deviations from the mean
9. try_to_fix_ylims = 1
	- Attempt to standardise y-limits of marginals

**Call Function**
10. `neurodtw.timeseries_correlation(series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims)`

## DTW Latency Difference - neurodtw.latency_difference()
Assess the latency difference between two conditions (i.e., within-subject effect) or between two groups (i.e., across-subject) effect of any signal of interest (in milliseconds).

Figures 3 & 4 of our paper show a two conditions analysis

Reference and query time series must be 2D matrices
I.e., data points (e.g., time) x subjects (i.e., replications)
Time series have to be of equal sizes
analysis_design determines whether you want to assess a within- or between-subjects latency effect (can only take “within” or “between” as input)
We provide the ERPs of correct and intrusion trials for users to explore this function


**Running timeseries_correlation using our ERPs**
__Make sure to enter your actual paths!__
1. Run Python, import the neurodtw package and load our ERPs
`import neurodtw`	
`from scipy.io import loadmat`
`data = loadmat("your/path/to/example_series_N2pcP3s")`

**Configure Variables**
2. analysis_design = "within"
3. query = data["N2pc_Intrusion"]
4. reference = data["N2pc_Correct"]
5. name_query = "Intrusion"
6. name_reference = "Correct"
7. units = "\u03BCV"
	- The units that your signals are measured in (in our case micro volts)
8. sampling_rate = 500
	- The number of data points per second in Hertz
9. filepath = "where/to/store/results/to"
10. permutations = 10000
	The number of permutations you would like to implement in statistical testing (we recommend >=10000)

**Call Function**
`neurodtw.latency_difference(analysis_design, query, reference, name_query, name_reference,units, sampling_rate, filepath, permutations)`

## Dependencies
*Python NeuroDTW requires the following toolboxes which are automatically installed via `pip install neurodtw`*
- Numpy
- Matplotlib
- Tslearn
- Scipy

## Tutorials & Walkthroughs (similar to those provided for Matlab NeuroDTW) will be provided after uploading Python NeuroDTW to pip