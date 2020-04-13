# SST-Rain-Response-Code
Code for analysis related to SST rain response paper, including updates to Deb's DYNAMO anlaysis and Falkor 2016 &amp; 2019 analysis

Deb's notes:
SST_rain_comparison_ID_rainpeaks.m  does the basic peak detection and collates met data from Edson.

SST_rain_comparison_bin_average_interim.m   plots some histograms and scatter
plots

SST_rain_comparison_with_anomaly.m    Collates variables into structure Var

SST_rain_comparison_peak_response.m   Adds times to max SST response, time to
recover, etc.

SST_rain_comparison_bin_average_more.m   bin averages to create composite
event - scales to beginning of event, and produces unscaled time as well (want
scaled!)

SST_rain_comparison_bin_average_long_short.m   bin averages to create
composite but selects out long events
