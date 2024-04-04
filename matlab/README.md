# Matlab NeuroDTW

Matlab NeuroDTW consists of two main scripts, dtw_timeseries_correlation.m and dtw_latency_difference.m - these are the general-purpose scripts enable the DTW analyses presented in *Transient Attention Gates Access Consciousness: Coupling N2pc and P3 Latencies using Dynamic Time Warping*.
The scripts used to replicate our analyses and figures as presented in our paper can be found in the *replicating_figures* folder.

## Dependencies
- Matlab’s Statistics & Machine Learning Toolbox
- Matlab’s Signal Processing Toolbox
- [Hline & vline](https://de.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline )
- [RGB](https://de.mathworks.com/matlabcentral/fileexchange/46872-intuitive-rgb-color-values-from-xkcd)
- [Suplabel](https://de.mathworks.com/matlabcentral/fileexchange/7772-suplabel )
- [nhist](https://de.mathworks.com/matlabcentral/fileexchange/27388-plot-and-compare-histograms-pretty-by-default)

## dtw_timeseries_correlation.m
- Peform a **DTW based bootstrap analysis** to assess the **temporal correlation between two time series** (Figure 5 of our paper)
- Time series must be 2D matrices
	- I.e., data points (e.g. time) x subjects (i.e., replications)


## Tests
Matlab NeuroDTW was tested with Matlab 2020b & 2023b on Windows. On a M1 Macbook the run_DTW local function caused crashes for my machine. If you have any insight on this please let me know. The rest worked as expected.