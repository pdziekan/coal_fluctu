from plot_coal_series import plot_coal_series

plots = ["rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/MarshallPalmer_R10/LCM/zero_dim/GA17_Np64e6_nx1_dt3var_MP_R10_HallDavis_cutoff2.5mm_stopr3mm_0"
data[directory + "/"]   = "LCM Np64e6 nx1 stopr3000, var dt HallDavis"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/MarshallPalmer_R10/LCM/GA17_from_libcloudphxx/GA17_Np1_nx400_eps0.1_dt0.01var_MP_R10_HallDavis_cutoff2.5mm_stopr3mm"
data[directory + "/"]   = "LCM Np1 nx400 stopr3000, var dt HallDavis"

plot_coal_series(plots, data, "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/MarshallPalmer_R10")
