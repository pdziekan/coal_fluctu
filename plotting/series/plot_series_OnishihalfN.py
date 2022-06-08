from plot_coal_series import plot_coal_series, plot_coal_series_diff
import matplotlib.pyplot as plt
import numpy as np

plots = ["tau", "nrain", "rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
# same labels means that binned statistics will be calculated from all these simulations
data[directory + "/"]   = "single coalescence cell"
data[directory + "_2/"] = "single coalescence cell"
data[directory + "_3/"] = "single coalescence cell"
data[directory + "_4/"] = "single coalescence cell"
data[directory + "_5/"] = "single coalescence cell"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "multiple coalescence cells"
data[directory + "_2/"] = "multiple coalescence cells"
data[directory + "_3/"] = "multiple coalescence cells"
data[directory + "_4/"] = "multiple coalescence cells"
data[directory + "_5/"] = "multiple coalescence cells"

fig, axs = plt.subplots(4, 1, figsize=(12,9))

for plot in plots:
  time, mean, mean_err, std_dev, std_dev_err = plot_coal_series(plot, data, fig, axs)

  # correction for a small difference in the way tau was calculated in different runs
  if plot == "tau":
    m = mean["multiple coalescence cells"]
    m = np.where(m != 0, m - m[1], m)
    mean["multiple coalescence cells"] = m

  plot_coal_series_diff("single coalescence cell", "multiple coalescence cells", time, mean, mean_err, std_dev, std_dev_err, fig, axs)

  axs[0].set_xticks([])
  axs[1].set_xticks([])
  axs[2].set_xticks([])
  axs[3].set_xlabel('time [s]')

  dirname = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/OnishihalfN"
  fig.tight_layout()
  fig.savefig(dirname+"_series_"+plot+"_stats.svg")
  for ax in axs:
    ax.clear()
