from plot_coal_series import plot_coal_series

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

plot_coal_series(plots, data, "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/OnishihalfN")
