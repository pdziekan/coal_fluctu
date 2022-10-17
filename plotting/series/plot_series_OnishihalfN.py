from plot_coal_series import plot_coal_series, plot_coal_series_diff, labeldict
import matplotlib.pyplot as plt
import numpy as np

# font setup
plt.rcParams.update({'font.size': 8})

# dpi of raster output
plt.rcParams.update({'savefig.dpi': 300})

plots = ["tau", "nrain", "rmax"]

data = {}


# opisy symulacji sa w tekscie na overleaf
directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/ares/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/1/"]   = "nx300 ares float"
data[directory + "/2/"]   = "nx300 ares float"
data[directory + "/3/"]   = "nx300 ares float"
data[directory + "/4/"]   = "nx300 ares float"
data[directory + "/5/"]   = "nx300 ares float"
data[directory + "/6/"]   = "nx300 ares float"
data[directory + "/7/"]   = "nx300 ares float"
data[directory + "/8/"]   = "nx300 ares float"
data[directory + "/9/"]   =  "nx300 ares GPU float, reszta double"
data[directory + "/10/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/11/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/12/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/13/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/14/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/15/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/16/"]   = "nx300 ares GPU float, reszta double"
data[directory + "/17/"]   = "nx300 ares double"
data[directory + "/18/"]   = "nx300 ares double"
data[directory + "/19/"]   = "nx300 ares double"
data[directory + "/20/"]   = "nx300 ares double"
data[directory + "/21/"]   = "nx300 ares double"
data[directory + "/22/"]   = "nx300 ares double"
data[directory + "/23/"]   = "nx300 ares double"
data[directory + "/24/"]   = "nx300 ares double"

#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/cuda_k_4/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
#data[directory + "/"]   = "nx300"
#data[directory + "_2/"] = "nx300"
#data[directory + "_3/"] = "nx300"
#data[directory + "_4/"] = "nx300"
#data[directory + "_5/"] = "nx300"

#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/nx30/GA17_Np1e3_nx30_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
#data[directory + "_part/"]   = "nx30"
#data[directory + "_2_part/"] = "nx30"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/nx3/GA17_Np1e6_nx3_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
#data[directory + "/"]   = "nx3"
#data[directory + "_2/"] = "nx3"
#data[directory + "_3/"] = "nx3"
#data[directory + "_4/"] = "nx3"

#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/cuda_k_4/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
## same labels means that binned statistics will be calculated from all these simulations
##data[directory + "/"]   = "nx1" # this set has strange dip in std dev of nrain at around 600 s; it differs in that it has sedi=1, but that doesnt matter
#data[directory + "_2/"] = "nx1"
#data[directory + "_3/"] = "nx1"
#data[directory + "_4/"] = "nx1"
#data[directory + "_5/"] = "nx1"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/ares/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/1/"] = "nx1 ares float"
data[directory + "/2/"] = "nx1 ares GPU float, reszta double"
data[directory + "/3/"] = "nx1 ares double"

#fig, axs = plt.subplots(3, 2, figsize=(8, 8))
fig, axs = plt.subplots(3, 2, figsize=(8, 8))

legend_loc = {}
legend_loc["tau"] = "upper center"
legend_loc["rmax"] = "lower right"
legend_loc["nrain"] = "lower right"

for plot in plots:
  time, mean, mean_err, std_dev, std_dev_err = plot_coal_series(plot, data, fig, axs)

#  if plot == "tau":
#    # correction for the way tau was calculated on ares ?
#    m = mean["nx300 ares"]
#    m = np.where(m != 0, m - m[1], m)
#    mean["nx300 ares"] = m
#    # correction for a small difference in the way tau was calculated in different runs
#    m = mean["nx300"]
#    m = np.where(m != 0, m - m[1], m)
#    mean["nx300"] = m
#    # y scale of plots with relative differences starting with values from the moment there was 0.1% of rain
#    t_1prom_it = np.argmax(m > 0.001) 
#    t_1prom = time["nx300"][t_1prom_it]
#    print("t_1prom:",t_1prom)
##    print("t_1prom_it:",t_1prom_it)

  t_1prom_it=0

#  plot_coal_series_diff(plot, "nx300", "S", "nx1", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx3", "M", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)

#  plot_coal_series_diff(plot, "nx300", "S", "nx1 1", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx1 2", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx1 3", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx1 4", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx1 5", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)

#  plot_coal_series_diff(plot, "nx300", "S", "nx3 1", "M", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx3 2", "M", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx3 3", "M", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "S", "nx3 4", "M", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)

#  plot_coal_series_diff(plot, "nx300", "L", "nx300 2", "S", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "L", "nx300 3", "S", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "L", "nx300 4", "S", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)
#  plot_coal_series_diff(plot, "nx300", "L", "nx300 5", "S", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)

  axs[0,0].tick_params(labelbottom=False)    
  axs[0,1].tick_params(labelbottom=False)    
  axs[1,0].tick_params(labelbottom=False)    
  axs[1,1].tick_params(labelbottom=False)    
#  axs[0].set_xticks([])
  axs[2,0].set_xlabel('time [s]')
  axs[2,1].set_xlabel('time [s]')

  for i,ax in enumerate(axs.flatten()):
    ax.text(0.03, 0.9, labeldict[i], fontsize=10, transform=ax.transAxes)


  dirname = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/OnishihalfN"

  handles, labels = axs[0,0].get_legend_handles_labels()
  axs[0,0].legend(handles, labels, loc=legend_loc[plot])
  axs[0,1].legend(handles, labels, loc=legend_loc[plot])
 # fig.legend()

  fig.tight_layout(pad=2)
  fig.savefig(dirname+"_series_"+plot+"_stats.pdf")
#  fig.savefig(dirname+"_series_"+plot+"_stats.svg")
#  fig.savefig(dirname+"_series_"+plot+"_stats_large.pdf")
  for ax in axs.flatten():
    ax.clear()
