import numpy as np
import timeseries_correlation
from scipy.io import loadmat


# # %% BETWEEN SUBJS DESIGN
# # used loadmat of scipy.io to load my example data mat file and:
# data = loadmat("/Users/mahan/sciebo/Matlab_N2pcP3/example_series_N2pcP3s")
# series_1 = data["P3_Correct"]
# series_2 = data["N2pc_Correct"]
# name_1 = "P3"
# name_2 = "N2pc"
# savepath = "/Users/mahan/sciebo/PythonCode/N2pcP3"
# num_boots = 10000
# outlier = 0
# try_to_fix_ylims = 1


# timeseries_correlation.main(
#     series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims
# )

# %% WITHIN SUBJS DESIGN
# used loadmat of scipy.io to load my example data mat file and:
data = loadmat("/Users/mahan/sciebo/Matlab_N2pcP3/example_series_N2pcP3s")
series_1 = data["P3_Correct"]
series_2 = data["P3_Intrusion"]
name_1 = "Correct"
name_2 = "Intrusion"
savepath = "/Users/mahan/sciebo/PythonCode/N2pcP3"
num_boots = 10000
outlier = 0
try_to_fix_ylims = 1


timeseries_correlation.main(
    series_1, series_2, name_1, name_2, savepath, num_boots, outlier, try_to_fix_ylims
)

