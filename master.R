#This is just for running the different "main" files in the right order,
#thereby producing all the results.
#
#Lei Zhao
#_________

#***DAN: either now or later, put in checkpoint functionality and also see if we can install reumanplatz as it existed on a particular date
#***DAN: we also need to review the wavelet tests

#***load the raw data and get the seasonal data ***
#1. choose right sites and get spatial averaged data
#2. get monthly data (averaged for each month)
#3. get seasonal data (averaged for each season based on monthly data)
#4. get linear interpolation for NAs
# finally get "Data_all_seasons_interpolate.RDS"

source("main_data.R") # <<<---becareful, time-consuming

#***detrend for each timeseries (output "Data_all_seasons_detrend.RDS")

source("main_data_detrend.R")

#*** average to get time series for shallow and deep ***
# output three versions: 
#1. original ("Data_TwoLayers.RDS")
#2. detrend ("Data_TwoLayers_detrend.RDS")
#3. clean: Box-Cox for each timeseries; detrend; standarizing variance ("Data_TwoLayers_clean.RDS")

source("main_data_TwoLayers.R")


#*** show clusters of correlation network ***
source("main_network.R")

#*** show timeseries for shallow and deep, and show correlations between shallow and deep ***
source("main_TwoLayers_ts_correlation.R")

#*** get spatial coherence of chla, T, and N between shallow and deep for near-shore and off-shore ***
source("main_TwoLayers_coherence.R")

#*** get spatial coherence of chla vs. T and chla vs. N for shallow or deep layers for near-shore and off-shore ***
source("main_TwoLayers_coherence_factors.R")

#*** test causal mechanisms ****
source("main_TwoLayers_mechanisms.R")

#*** application: show community stability ***
source("main_TwoLayers_stability.R")


