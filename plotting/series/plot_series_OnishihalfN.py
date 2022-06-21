from plot_coal_series import plot_coal_series, plot_coal_series_diff, labeldict
import matplotlib.pyplot as plt
import numpy as np

# font setup
plt.rcParams.update({'font.size': 8})

# dpi of raster output
plt.rcParams.update({'savefig.dpi': 300})

plots = ["tau", "nrain", "rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/nx300/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "nx300"
data[directory + "_2/"] = "nx300"
data[directory + "_3/"] = "nx300"
data[directory + "_4/"] = "nx300"
data[directory + "_5/"] = "nx300"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/nx30/GA17_Np1e3_nx30_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "_part/"]   = "nx30"
data[directory + "_2_part/"] = "nx30"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/nx3/GA17_Np1e6_nx3_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "nx3"
data[directory + "_2/"] = "nx3"
data[directory + "_3/"] = "nx3"
data[directory + "_4/"] = "nx3"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
# same labels means that binned statistics will be calculated from all these simulations
#data[directory + "/"]   = "nx1" # this set has strange dip in std dev of nrain at around 600 s; it differs in that it has sedi=1, but that doesnt matter
data[directory + "_2/"] = "nx1"
data[directory + "_3/"] = "nx1"
data[directory + "_4/"] = "nx1"
data[directory + "_5/"] = "nx1"

#fig, axs = plt.subplots(3, 2, figsize=(8, 8))
fig, axs = plt.subplots(3, 2, figsize=(8, 8))

legend_loc = {}
legend_loc["tau"] = "upper center"
legend_loc["rmax"] = "lower right"
legend_loc["nrain"] = "lower right"

for plot in plots:
  time, mean, mean_err, std_dev, std_dev_err = plot_coal_series(plot, data, fig, axs)

#  if plot == "tau":
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
