import numpy as np
import matplotlib.pyplot as plt

plots = ["tau", "nrain", "rmax"]

data = {}

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/zero_dim/GA17_Np27e6_nx1_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "LCM Np27e6 nx1 stopr300, var dt HallDavis (1)"
data[directory + "_2/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis (2)"
data[directory + "_3/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis (3)"
data[directory + "_4/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis (4)"
data[directory + "_5/"] = "LCM Np27e6 nx1 stopr300, var dt HallDavis (5)"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/OnishihalfN/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps0.1_dt1var_OnishihalfN_HallDavis_stopr300"
data[directory + "/"]   = "LCM Np1 nx300 stopr300, var dt HallDavis (1)"
data[directory + "_2/"] = "LCM Np1 nx300 stopr300, var dt HallDavis (2)"
data[directory + "_3/"] = "LCM Np1 nx300 stopr300, var dt HallDavis (3)"
data[directory + "_4/"] = "LCM Np1 nx300 stopr300, var dt HallDavis (4)"
data[directory + "_5/"] = "LCM Np1 nx300 stopr300, var dt HallDavis (5)"

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
  fig.savefig("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/OnishihalfN_series_"+plot+".png")
