import numpy as np
import matplotlib.pyplot as plt

plots = ["tau", "nrain", "rmax"]

data = {}

# --- Wang data ---
#directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/synth_turb_from_separate_code/ST_Np1e3_nx20_eps1_dt0.1var_Wang/"
#data = {
#         directory+ "/1/" : "LCM Np1e3 nx20 rmr250 (1)"
#        ,directory+"/2/" : "LCM Np1e3 nx20 rmr250 (2)"
#        ,directory+"/3/" : "LCM Np1e3 nx20 rmr250 (3)"
#       }
#
#directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/synth_turb_from_separate_code/ST_Np8e6_nx1_eps1_dt0.1var_Wang/"
#data[directory+ "/1/"] = "LCM Np8e6 nx1 rmr250 (1)"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/synth_turb_from_separate_code/ST_Np1e6_nx1_eps1_dt0.1var_Wang_normr/"
#data[directory+ "/1/"] = "LCM Np1e6 nx1 normr (1)"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/synth_turb_from_separate_code/ST_Np1e7_nx1_eps1_dt0.1var_Wang_normr/"
#data[directory+ "/1/"] = "LCM Np1e7 nx1 normr (1)"
#data[directory+ "/2/"] = "LCM Np1e7 nx1 normr (2)"
#
directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx454_eps1_dt0.1var_Wang_HallDavis_rmr100"
data[directory+ "/"] = "LCM Np1e0 nx454 rmr100, var dt 0.1, GA17 libcloudphxx, HallDavis"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx100_eps1_dt0.1var_Wang_HallDavis_rmr100"
#data[directory+ "/"] = "LCM Np1e0 nx100 rmr100, var dt 0.1, GA17 libcloudphxx, HallDavis"
###
directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/GA17_from_libcloudphxx/GA17_Np1e6_nx1_eps1_dt0.1var_Wang_HallDavis_rmr100"
data[directory+ "/"] = "LCM Np1e6 nx1 rmr100, var dt 0.1, GA17 libcloudphxx, HallDavis"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/LCM/GA17_from_libcloudphxx/GA17_Np1e6_nx1_eps1_dt0.1var_Wang_HallDavis_normr"
#data[directory+ "/"] = "LCM Np1e6 nx1 normr, var dt 0.1, GA17 libcloudphxx, HallDavis"

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

  # plot tau from Bott EFM

## -- Wang ---
  if plot=="tau":
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/Wang_HallDavis/dt0.10_rq0 9.3_xmw 1.0_scal10.0_isw3_cut 20._rmr****_rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_normr dt0.1")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/Wang_HallDavis/dt 0.1_rq0 9.3_xmw 1.0_scal10.0_isw3_cut 20._rmr250._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr250 dt0.1")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/Wang_HallDavis/dt 0.1_rq0 9.3_xmw 1.0_scal10.0_isw3_cut 20._rmr100._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr100 dt0.1")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/Wang_HallDavis/dt0.10_rq0 9.3_xmw 1.0_scal15.0_isw3_cut 20._rmr100._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr100 dt0.1 scal15")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/data/EFM/Wang_HallDavis/dt0.10_rq0 9.3_xmw 1.0_scal40.0_isw3_cut 20._rmr100._1000bins_rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr100 dt0.1 scal40 1000bins")
  plt.xlim([2000, 5000])
  
  plt.legend()
  fig.savefig("/home/piotr/praca/coal_fluctu_dim/Smolu_min_cell/img/Wang_HallDavis/series_"+plot+"_Wang.png")


