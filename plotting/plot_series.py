import numpy as np
import matplotlib.pyplot as plt

plots = ["tau", "nrain", "rmax"]

#data = {"/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/ST_Np1e4_nx1_eps1_dt0.1var_Wang/1/" : "LCM Np1e4 nx1 normr Onishi"}

directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/ST_Np1e3_nx20_eps1_dt0.1var_Wang/"
data = {
         directory+ "/1/" : "LCM Np1e3 nx20 rmr250 (1)"
        ,directory+"/2/" : "LCM Np1e3 nx20 rmr250 (2)"
        ,directory+"/3/" : "LCM Np1e3 nx20 rmr250 (3)"
       }

directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/ST_Np8e6_nx1_eps1_dt0.1var_Wang/"
data[directory+ "/1/"] = "LCM Np8e6 nx1 rmr250 (1)"

directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/ST_Np1e6_nx1_eps1_dt0.1var_Wang_normr/"
data[directory+ "/1/"] = "LCM Np1e6 nx1 normr (1)"

directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/ST_Np1e7_nx1_eps1_dt0.1var_Wang_normr/"
data[directory+ "/1/"] = "LCM Np1e7 nx1 normr (1)"
data[directory+ "/2/"] = "LCM Np1e7 nx1 normr (2)"

for plot in plots:
  fig, ax = plt.subplots()

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
      ax.plot(time, series, label=data[pre]+'('+str(idx)+')')

  # plot tau from Bott EFM
  if plot=="tau":
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/rain_ratio_normr.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_normr")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/rain_ratio_rmr250.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr250")

  
  plt.xlim([2000, 5000])
  plt.legend()
  fig.savefig("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/img/series_"+plot+".png")


