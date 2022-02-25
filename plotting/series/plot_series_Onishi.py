import numpy as np
import matplotlib.pyplot as plt

plots = ["tau", "nrain", "rmax"]

data = {}

# --- Onishi data ---
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/zero_dim/ST_Np1e4_nx1_eps1_dt0.1var_Onishi_HallDavis_normr/"
#data = {
#         directory+ "/1/" : "LCM Np1e4 nx1 normr (1)",
#         directory+ "/2/" : "LCM Np1e4 nx1 normr (2)"
#       }
#
directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/zero_dim/ST_Np1e5_nx1_eps1_dt0.1var_Onishi_HallDavis_normr/"
data[directory+ "/1/"] = "LCM ST Np1e5 nx1 normr, var dt HallDavis (1)"
#data[directory+ "/2/"] = "LCM Np1e5 nx1 normr, var dt HallDavis (2)"


directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/zero_dim/ST_Np1e5_nx1_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
data[directory+ "/1/"] = "LCM ST Np1e5 nx1 rmr250, var dt HallDavis (1)"
#data[directory+ "/2/"] = "LCM ST Np1e5 nx1 rmr250, var dt HallDavis (2)"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/zero_dim/GA17_Np1e5_nx1_eps1_dt0.1var_Onishi_HallDavis_rmr100/"
data[directory+ "/"] = "LCM GA Np1e5 nx1 rmr100, var dt HallDavis"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/zero_dim/GA17_Np1e6_nx1_eps1_dt0.1var_Onishi_HallDavis_rmr100/"
data[directory+ "/"] = "LCM GA Np1e6 nx1 rmr100, var dt HallDavis"

#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/synth_turb_from_separate_code/ST_Np1e2_nx10_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
#data[directory+ "/1/"] = "LCM Np1e2 nx10 rmr250 (1)"
#data[directory+ "/2/"] = "LCM Np1e2 nx10 rmr250 (2)"

#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/synth_turb_from_separate_code/ST_Np1e0_nx46_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
#data[directory+ "/1/"] = "LCM Np1e0 nx46 rmr250 (1)"
#data[directory+ "/2/"] = "LCM Np1e0 nx46 rmr250 (2)"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/synth_turb_from_separate_code/ST_Np1e0_nx100_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
#data[directory+ "/1/"] = "LCM Np1e0 nx100 rmr250 (1)"
#
directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/synth_turb_from_separate_code/ST_Np1e0_nx200_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
data[directory+ "/1/"] = "LCM ST Np1e0 nx200 rmr250 (1)"

directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/synth_turb_from_libcloudphxx/ST_Np1e0_nx300_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
data[directory+ "/"] = "LCM ST_lib Np1e0 nx300 rmr250, var dt, HallDavis"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/synth_turb_from_libcloudphxx/ST_Np1e0_nx300_eps1_dt0.1_Onishi_HallDavis_rmr250/"
#data[directory+ "/"] = "LCM Np1e0 nx300 rmr250, const dt, synth turb libcloudphxx, HallDavis"
#
directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps1_dt0.1var_Onishi_HallDavis_rmr250/"
data[directory+ "/"] = "LCM GA Np1e0 nx300 rmr250, vat dt, GA17 libcloudphxx, HallDavis"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps1_dt0.1_Onishi_HallDavis_rmr250/"
#data[directory+ "/"] = "LCM Np1e0 nx300 rmr250, const dt, GA17 libcloudphxx, HallDavis"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps1_dt0.1_Onishi_Hall_rmr250/"
#data[directory+ "/"] = "LCM Np1e0 nx300 rmr250, const dt, GA17 libcloudphxx, Hall"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx300_eps1_dt0.1_Onishi_Hall_rmr250_kappa0.1/"
#data[directory+ "/"] = "LCM Np1e0 nx300 rmr250, const dt, GA17 libcloudphxx, Hall, kappa 0.1"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx454_eps1_dt0.1_Onishi_Hall_rmr250/"
#data[directory+ "/"] = "LCM Np1e0 nx454 rmr250, const dt, GA17 libcloudphxx, Hall"
#
#directory = "/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/LCM/GA17_from_libcloudphxx/GA17_Np1e0_nx454_eps1_dt0.1_Onishi_Hall_normr/"
##data[directory+ "/"] = "LCM Np1e0 nx454 no rmr, const dt, GA17 libcloudphxx, Hall"

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

# -- Onishi ---
  if plot=="tau":
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/EFM/dt0.10_rq015.0_xmw 2.0_scal10.0_isw3_cut 40._rmr****_rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_normr dt0.1")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/EFM/dt0.10_rq015.0_xmw 2.0_scal10.0_isw3_cut 40._rmr250._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr250 dt0.1")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/EFM/dt0.10_rq015.0_xmw 2.0_scal10.0_isw3_cut 40._rmr100._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr100 dt0.1")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/EFM/dt0.01_rq015.0_xmw 2.0_scal10.0_isw3_cut 40._rmr100._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr100 dt0.01")
    efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/data/Onishi/EFM/dt0.01_rq015.0_xmw 2.0_scal15.0_isw3_cut 40._rmr100._rain_ratio.out", unpack=True)
    ax.plot(efm_data[0], efm_data[1], label="EFM_rmr100 dt0.01 scal15")
  plt.xlim([0, 800])

  
  plt.legend()
  fig.savefig("/home/piotr/praca/coal_fluctu_dim/well_mixed_cell_size/img/series/Onishi_series_"+plot+".png")


