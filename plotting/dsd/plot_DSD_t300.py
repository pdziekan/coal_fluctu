import numpy as np
import matplotlib.pyplot as plt

data_labels = {}
data_colors = {}
data_nsd = {}

directory_base = "/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/LCM/Onishi/"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = "DT 0.1 ConstMulti1"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300/"] = 64e6
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_ConstMulti1_Ens10_T300_RngSeedInit44/"] = "DT 0.1 ConstMulti1 RngSeedInit"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = "DT 0.1 SstpCoal1  Sd1e5Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e5Tail_Ens100_T300/"] = 1e5

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] = "DT 0.1 SstpCoal1  Sd1e4Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e4Tail_Ens100_T300/"] = 1e4

# nie ma roznic miedzy tymi trzema, dt i sstp_coal nie maja wplywu
data_labels[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.01 SstpCoal1  Sd1e3Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.01_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"] = 1e3
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e3Tail_Ens100_T300/"]  = "DT 0.1  SstpCoal1  Sd1e3Tail"
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal10_Onishi_SdConc1e3Tail_Ens100_T300/"] = "DT 0.1  SstpCoal10 Sd1e3Tail"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e2Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens1e4_T300/"] = 1e2
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e2Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e2Tail RngSeedInit"

data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] =               "DT 0.1  SstpCoal1  Sd1e1Tail"
data_nsd[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens1e4_T300/"] = 1e1
#data_labels[directory_base + "GA17_Np64e6_nx1_dt0.1_SstpCoal1_Onishi_SdConc1e1Tail_Ens100_T300_RngSeedInit44/"] = "DT 0.1  SstpCoal10 Sd1e1Tail RngSeedInit"

four_over_three_pi_rhow = 4./3. * np.pi * 1e3 # [kg/m3]
volume = 0.451 # cell volume [m3]

fig, ax = plt.subplots(1, 2, figsize=(24,9))

for pre in data_labels:
  LCM_r         = np.zeros(0)
  LCM_m_r_end    = np.zeros(0)
  LCM_m_r_std_dev_end  = np.zeros(0)

  fs = open(pre+"size_spectr_mean.dat","r")
  rows = [x.split() for x in fs.readlines()]
  for row in rows[1:]: # first row is data description
    LCM_r                = np.append(LCM_r,                float(row[0])) # [um]
    LCM_m_r_end          = np.append(LCM_m_r_end,          float(row[3])) # same
    LCM_m_r_std_dev_end  = np.append(LCM_m_r_std_dev_end,  float(row[4])) # same

  LCM_dlogr = np.log(LCM_r[1]) - np.log(LCM_r[0]) # [log(um)]
  LCM_r *= 1e-6 # [m]

  LCM_m_logr_end = LCM_m_r_end * LCM_r # [g/m3 / unit(log[um])] 
  LCM_m_logr_std_dev_end = LCM_m_r_std_dev_end * LCM_r # [g/m3 / unit(log[um])] 
  LCM_M_end = LCM_m_logr_end * LCM_dlogr * volume * 1e-3 # mass of droplets of this size in the cell [kg]
  LCM_N_end = LCM_M_end / four_over_three_pi_rhow / np.power(LCM_r,3) # number of droplets of this size in the cell [1]
  print("total number of droplets in the cell in LCM at the end: ", np.sum(LCM_N_end))
  print("total mass of droplets in the cell in LCM at the end: ", np.sum(LCM_M_end))

# plot m(log r)
  ax[0].plot(LCM_r * 1e6, LCM_m_logr_end, label= data_labels[pre])
  ax[1].plot(LCM_r * 1e6, LCM_m_logr_std_dev_end, label= data_labels[pre])
# std dev scaled as sqrt(N_SD)
#  ax[1].plot(LCM_r * 1e6, LCM_m_logr_std_dev_end  * np.sqrt(data_nsd[pre]) / np.sqrt(64e6), label='end ' + data_labels[pre])

# plot m(r)
#  ax[0].plot(LCM_r * 1e6, LCM_m_r_end, label='end ' + data_labels[pre])
#  ax[1].plot(LCM_r * 1e6, LCM_m_r_std_dev_end, label='end ' + data_labels[pre])

# plot N (number of droplets in the bin)
#  ax[0].plot(LCM_r * 1e6, LCM_N_end, label='end ' + data_labels[pre])
  

#  ax[0].plot(LCM_r, LCM_m_r_init, label='init ' + data_labels[pre]) # LCM data is mass density m(r) with LCM_rius in meters
#  ax[1].plot(LCM_r, LCM_m_r_std_dev_init, label='init ' + data_labels[pre])
 # ax[1].plot(LCM_r, LCM_m_r_std_dev_init * np.sqrt(data_nsd[pre]) / np.sqrt(64e6), label='init ' + data_labels[pre])
 # ax[1].plot(LCM_r, LCM_m_r_std_dev_end  * np.sqrt(data_nsd[pre]) / np.sqrt(64e6), label='end ' + data_labels[pre])

# EFM results (Smoluchowski)
efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/EFM/dt 0.1_rq015.0_xmw 2.0_scal10.0_isw3_cut****_rmr****_tmax301._boplot00.out", unpack=True)
#efm_data = np.genfromtxt("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/data/EFM/dt 0.1_rq015.0_xmw 2.0_scal15.0_isw3_cut****_rmr****_tmax301._boplot00.out", unpack=True)
# efm_data[0] - radius[um]
# efm_data[1] - mass density m(ln r) [kg/m3] with radius in microns
EFM_r = efm_data[0] * 1e-6 # [m]
EFM_m_logr = efm_data[1] # [kg/m3 / unit(log[um])]
EFM_dlogr = np.log(efm_data[0][1]) - np.log(efm_data[0][0]) # [log(um)]
EFM_m_r = EFM_m_logr / EFM_r # mass density m(r) [kg/m3 / m]
EFM_M = EFM_m_logr * EFM_dlogr * volume # mass of droplets of this size in the cell [kg]
EFM_N = EFM_M / four_over_three_pi_rhow / np.power(EFM_r,3) # number of droplets of this size in the cell [1]
EFM_M_std_dev = four_over_three_pi_rhow * np.power(EFM_r,3) * np.sqrt(EFM_N) # std_dev of total mass of droplets of this size in the whole cell [kg]
# both give the same
EFM_m_logr_std_dev = EFM_M_std_dev / EFM_dlogr / volume # std_dev of mass density m(log r) [kg/m3]
#EFM_m_logr_std_dev = np.sqrt(four_over_three_pi_rhow * np.power(EFM_r,3) * EFM_m_logr / EFM_dlogr / volume) # std_dev of mass density m(log r) [kg/m3]

# plot m(log r)
EFM_m_logr_init = EFM_m_logr[:int(EFM_N.size/2)]
EFM_m_logr_end = EFM_m_logr[int(EFM_N.size/2):]
EFM_m_logr_std_dev_init = EFM_m_logr_std_dev[:int(EFM_N.size/2)]
EFM_m_logr_std_dev_end = EFM_m_logr_std_dev[int(EFM_N.size/2):]
ax[0].plot(EFM_r[:int(EFM_N.size/2)] * 1e6, EFM_m_logr_end * 1e3, label="EFM") # radius in [um], mass density in [g/m3 / unit(log[um])]
ax[1].plot(EFM_r[:int(EFM_N.size/2)] * 1e6, EFM_m_logr_std_dev_end * 1e3, label="EFM") # radius in [um], mass density in [g/m3 / unit(log[um])]

# plot N (number of droplets in the bin)
#ax[0].plot(EFM_r * 1e6, EFM_N, label='end ' + data_labels[pre])

#  ax[0].plot(EFM_r * 1e6, efm_data[1] / efm_data[0] * 1e6 * 1e3, label="EFM") 
#  ax[1].plot(EFM_r * 1e6, np.sqrt(four_over_three_pi_rhow * np.power(efm_data[0] * 1e-6,3) * EFM_M), label="EFM") # std_dev of total mass of droplets of this size in the whole cell [kg]
#  ax[0].plot(efm_data[0], efm_data[1] * EFM_dlogr / four_over_three_pi_rhow / np.power(efm_data[0] * 1e-6,3), label="EFM") # EFM data is mass density m(ln r) [kg/m3] with radius in microns; LCM data is m(r) with radius in meters

print("total number of droplets in the cell in EFM at the end: ", np.sum(EFM_N[int(EFM_N.size/2):])) # EFM_N contains init spectrum and final spectrum
print("total mass of droplets in the cell in EFM at the end: ", np.sum(EFM_M[int(EFM_M.size/2):]))

#plt.title('Time at which largest droplet exceeds 3 mm radius and size of the droplet')
ax[0].set_xlabel('radius [um]')
ax[1].set_xlabel('radius [um]')
#ax[0].set_ylabel('mean mass density m(r) [g/m^3 / m]')
ax[0].set_ylabel('mean mass density m(log r) [g/m^3 / unit(log[um])')
ax[1].set_ylabel('standard deviation of mass density m(log r) [g/m^3 / unit(log[um])')

ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
#ax[0].set_ylim(1e-7,)
ax[1].set_ylim(1e-8,)

#ax[0].set_xlim(40,)
#ax[0].set_ylim(0,0.30)

plt.grid(axis='y', color='0.5')
plt.legend()
fig.savefig("/home/piotr/praca/coal_fluctu_dim/LCM_DSD_fluctuations/img/DSD_t300.png")
plt.show()
