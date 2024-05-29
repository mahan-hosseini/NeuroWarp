# Python NeuroWarp

Python NeuroWarp is provided as a package that can be installed via PyPi. It consists of two functions that can be called after importing neurowarp: timeseries_correlation and latency_difference -  these are the general-purpose scripts enable the DTW analyses presented in *Transient Attention Gates Access Consciousness: Coupling N2pc and P3 Latencies using Dynamic Time Warping*.

## Installation

We recommend that you create a new virtual environment for our module via:
1. Open a terminal / cmd.exe window, navigate (via `cd`) to the directory you want to create the environment and enter:
`python -m venv neurowarp_env`
2. Activate the neurowarp_env via:

	**Windows:** `path\to\neurowarp_env\Scripts\activate`

	**MacOS:** `source neurowarp_env/bin/activate`

3. Install neurowarp via pip (enter while neurowarp_env is active):
`pip install neurowarp`
4. The functions can be used as explained below after importing the module

For more detail on virtual environments & pip [click here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)


## DTW Temporal Correlation - neurowarp.timeseries_correlation()
Perform a **DTW based bootstrap analysis** to assess the **temporal correlation between two time series** (Figure 5 of our paper).

### Important Notes
- Time series must be **2D matrices**
	- I.e., data points (e.g. time) x subjects (i.e., replications)
- We provide the ERPs of correct and intrusion trials for users to explore this function

### Running timeseries_correlation using our ERPs
*Enter the following into Python and make sure to enter your actual paths!*
1. `import neurowarp`	
2. `from scipy.io import loadmat`
3. `data = loadmat("your/path/to/example_series_N2pcP3s")`
4. `series_1 = data["P3_Correct"]`
5. `series_2 = data["N2pc_Correct"]`
6. `name_1 = "P3"`
7. `name_2 = "N2pc"`
8. `savepath = "where/to/store/results/to"`
9. `num_boots = 10000`
	- The number of bootstrap samples that you want to implement
10. `outlier = 0`
	- Exclude outliers if their DTW area is +-5 standard deviations from the mean
11. `try_to_fix_ylims = 1`
	- Attempt to standardise y-limits of marginals
12. `neurowarp.timeseries_correlation(series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims)`

*Note that the figure will look slightly different to that of our paper due to different x/y limits. See the replicate_figures folder if you want to replicate our figure as it was printed.*

## DTW Latency Difference - neurowarp.latency_difference()
Assess the **latency difference** between **two conditions** (i.e., within-subjects effect) or between **two groups** (i.e., across-subjects effect) of any signal of interest (in milliseconds).

*Figures 3 & 4 of our paper show a two conditions analysis*

### Important Notes
- Reference and query time series must be **2D matrices**
	- I.e., data points (e.g., time) x subjects (i.e., replications)
	- Time series have to be of **equal sizes**
- **analysis_design** determines whether you want to assess a within- or between-subjects latency effect (can only take “within” or “between” as input)
- We provide the ERPs of correct and intrusion trials for users to explore this function

### Running timeseries_correlation using our ERPs
*Enter the following into Python and make sure to enter your actual paths!*
1. `import neurowarp`	
2. `from scipy.io import loadmat`
3. `data = loadmat("your/path/to/example_series_N2pcP3s")`
4. `analysis_design = "within"`
5. `query = data["N2pc_Intrusion"]`
6. `reference = data["N2pc_Correct"]`
7. `name_query = "Intrusion"`
8. `name_reference = "Correct"`
9. `units = "\u03BCV"`
	- The units that your signals are measured in (in our case micro volts)
10. `sampling_rate = 500`
	- The number of data points per second in Hertz
11. `savepath = "where/to/store/results/to"`
12. `permutations = 10000`
	- The number of permutations you would like to implement in statistical testing (we recommend >=10000)
13. `neurowarp.latency_difference(analysis_design, query, reference, name_query, name_reference,units, sampling_rate, savepath, permutations)`

## Dependencies
*Python NeuroWarp requires the following toolboxes which are automatically installed via `pip install neurowarp`*
- Numpy
- Matplotlib
- Tslearn
- Scipy

## Tests
Python NeuroWarp was tested with Python 3.10.9.