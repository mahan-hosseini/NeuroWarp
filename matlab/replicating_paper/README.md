# Matlab NeuroDTW - Replicating our Figures

## Dependencies
*Matlab NeuroDTW requires the following toolboxes to be installed and available in your Matlab path:*
- Matlab’s Statistics & Machine Learning Toolbox
- Matlab’s Signal Processing Toolbox
- [Hline & vline](https://de.mathworks.com/matlabcentral/fileexchange/1039-hline-and-vline )
- [RGB](https://de.mathworks.com/matlabcentral/fileexchange/46872-intuitive-rgb-color-values-from-xkcd)
- [Suplabel](https://de.mathworks.com/matlabcentral/fileexchange/7772-suplabel )
- [nhist](https://de.mathworks.com/matlabcentral/fileexchange/27388-plot-and-compare-histograms-pretty-by-default)

**Please download our ERPs and code-files to the same folder and configure your MATLAB current directory to that folder!**

## Figure 2 - ERPs 
1. load("LPF_ERP_structs.mat")
2. plotpath = "enter/the/path/to/save/plots/to/" (end with a "/"!)
3. ´plotERPs(LPF_P3s, LPF_N2pcs, plotpath, 1, 0)´

## Figure 3 & 4 - DTW for ERP Latency Differences 
1. experiment = "Experiment 1"
2. filepath = "enter/the/path/to/LPF_ERP_structs.mat/" (end with a "/"!)
3. dtwpath = filepath (or another path if you want to save the output somewhere else)
4. permutations = 10000
5. ´runDTW(experiment, filepath, dtwpath, permutations, 1)´

## Figure 5 - DTW Bootstrap Analysis of N2pc-P3 correlation
1. filepath = "enter/the/path/to/LPF_ERP_structs.mat/" (end with a "/"!)
2. P3timewin = [250 800]
3. N2pctimewin = [150 400]
4. num_boots = 10000
5. outlier = 0 (set to 1 to replicate the results mentioned in the "large P3 DTW areas" section
6. ´correlateP3andN2pc(filepath, P3timewin, N2pctimewin, num_boots, outlier)´

## Figure 6 - DTW Bootstrap N2pc & P3 Marginal Distribution Comparison Plots
1. filepath = "enter/the/path/to/Marginal Distributions.mat/" (end with a "/"!)
Marginal Distributions.mat is created by correlateP3andN2pc function!
2. ´compareMarginals(filepath, 0)´

## Figure 7 - DTW Simulations with (Human) Noise
1. make sure "meanpower.mat" is on filepath
2. filepath = "enter/the/path/to/LPF_ERP_structs.mat&meanpower.mat/" (end with a "/"!)
3. dtwpath = filepath
4. ´DTWSimulations(filepath, dtwpath)´
