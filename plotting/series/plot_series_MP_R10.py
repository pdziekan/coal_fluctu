import numpy as np
import matplotlib.pyplot as plt

plots = ["tau", "nrain", "rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/MarshallPalmer_R10/LCM/zero_dim/GA17_Np64e6_nx1_dt3var_MP_R10_HallDavis_cutoff2.5mm_stopr3mm_0"
data[directory + "/"]   = "LCM Np64e6 nx1 stopr3000, var dt HallDavis (1)"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/MarshallPalmer_R10/LCM/GA17_from_libcloudphxx/GA17_Np1_nx400_eps0.1_dt0.01var_MP_R10_HallDavis_cutoff2.5mm_stopr3mm_part"
data[directory + "/"]   = "LCM Np1 nx400 stopr3000, var dt HallDavis (1)"

for plot in plots:
  fig, ax = plt.subplots(figsize=(12,9))

  for pre in data:
    print(plot)
    print(pre)
    ft = open(pre+"time.dat","r")
    fd = open(pre+plot+".dat","r")
    all_time=[x.split() for x in ft.readlines()]
    all_data=[x.split() for x in fd.readlines()]
    for idx, (row_t, row_d) in enumerate(zip(all_time, all_data)):
      time = np.array(row_t, dtype=np.float32)
      series = np.array(row_d, dtype=np.float32)
      ax.plot(time, series[0:len(time)], label=data[pre]+'('+str(idx)+')')

  #plt.xlim([0, 800])

  
  #plt.legend()
  fig.savefig("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/MarshallPalmer_R10_series_"+plot+".png")
