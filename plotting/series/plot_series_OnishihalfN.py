from plot_coal_series import plot_coal_series, plot_coal_series_diff, labeldict
import matplotlib.pyplot as plt
import numpy as np

# font setup
plt.rcParams.update({'font.size': 8})

plots = ["tau", "nrain", "rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
# same labels means that binned statistics will be calculated from all these simulations
data[directory + "/"]   = "large coalescence cells"
data[directory + "_2/"] = "large coalescence cells"
data[directory + "_3/"] = "large coalescence cells"
data[directory + "_4/"] = "large coalescence cells"
data[directory + "_5/"] = "large coalescence cells"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "small coalescence cells"
data[directory + "_2/"] = "small coalescence cells"
data[directory + "_3/"] = "small coalescence cells"
data[directory + "_4/"] = "small coalescence cells"
data[directory + "_5/"] = "small coalescence cells"

fig, axs = plt.subplots(3, 2, figsize=(6.3,11))

legend_loc = {}
legend_loc["tau"] = "upper center"
legend_loc["rmax"] = "lower right"
legend_loc["nrain"] = "lower right"

for plot in plots:
  time, mean, mean_err, std_dev, std_dev_err = plot_coal_series(plot, data, fig, axs)

  if plot == "tau":
    # correction for a small difference in the way tau was calculated in different runs
    m = mean["small coalescence cells"]
    m = np.where(m != 0, m - m[1], m)
    mean["small coalescence cells"] = m
    # y scale of plots with relative differences starting with values from the moment there was 0.1% of rain
    t_1prom_it = np.argmax(m > 0.001) 
    t_1prom = time["small coalescence cells"][t_1prom_it]
#    print("t_1prom:",t_1prom)
#    print("t_1prom_it:",t_1prom_it)

  plot_coal_series_diff(plot, "small coalescence cells", "S", "large coalescence cells", "L", time, mean, mean_err, std_dev, std_dev_err, t_1prom_it, fig, axs)

#  axs[0].tick_params(labelbottom=False)    
#  axs[1].tick_params(labelbottom=False)    
#  axs[2].tick_params(labelbottom=False)    
#  axs[3].tick_params(labelbottom=False)    
#  axs[0].set_xticks([])
#  axs[4].set_xlabel('time [s]')

#  for i,ax in enumerate(axs):
#    ax.text(0.01, 0.9, labeldict[i], fontsize=10, transform=ax.transAxes)


  dirname = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/OnishihalfN"

  handles, labels = axs[0].get_legend_handles_labels()
  axs[0].legend(handles, labels, loc=legend_loc[plot])
  axs[1].legend(handles, labels, loc=legend_loc[plot])
 # fig.legend()

  fig.tight_layout(pad=2)
  fig.savefig(dirname+"_series_"+plot+"_stats.pdf")
  fig.savefig(dirname+"_series_"+plot+"_stats.svg")
  for ax in axs:
    ax.clear()
