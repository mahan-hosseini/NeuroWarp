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
4. The functions can be used as explained below after importing the module via (while neurodtw_env is active):
`python`
`import neurodtw`

For more detail on virtual environments & pip [click here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)

## Dependencies
*Python NeuroDTW requires the following toolboxes which are automatically installed via ˚pip install neurodtw˚
- Numpy
- Matplotlib
- Tslearn
- Scipy

## Tutorials & Walkthroughs (similar to those provided for Matlab NeuroDTW) will be provided after uploading Python NeuroDTW to pip