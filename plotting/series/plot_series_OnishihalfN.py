from plot_coal_series import plot_coal_series

plots = ["tau", "nrain", "rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
# same labels means that binned statistics will be calculated from all these simulations
data[directory + "/"]   = "LCM Np27e6 nx1 stopr300, var dt HallDavis"
data[directory + "_2/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis"
data[directory + "_3/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis"
data[directory + "_4/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis"
data[directory + "_5/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "LCM Np1 nx300 stopr300, var dt HallDavis"
data[directory + "_2/"] = "LCM Np1 nx300 stopr300, var dt HallDavis"
data[directory + "_3/"] = "LCM Np1 nx300 stopr300, var dt HallDavis"
data[directory + "_4/"] = "LCM Np1 nx300 stopr300, var dt HallDavis"
data[directory + "_5/"] = "LCM Np1 nx300 stopr300, var dt HallDavis"

plot_coal_series(plots, data, "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/OnishihalfN")
