from plot_DSD_common import mean_DSD

datasets = [
"/home/piotr/tmp/04_07_2023_GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/",
"/home/piotr/tmp/05_07_2023_GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/2/",
"/home/piotr/tmp/05_07_2023_GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/3/",
]

r, m_pre, m_pre_std, m_post, m_post_std = mean_DSD(datasets)

fo = open("/home/piotr/tmp/GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/size_spectr_mean.dat", "w")
fo.write('# radius[um]  mass_density@init_ensemble_mean_of_mean_from_cells[g/m3 / m] mass_density@init_ensemble_std_dev_of_mean_from_cells[g/m3 / m] mass_density@end_ensemble_mean_of_mean_from_cells[g/m3 / m] mass_density@end_ensemble_std_dev_of_mean_from_cells[g/m3 / m]\n')
for l in zip(r, m_pre, m_pre_std, m_post, m_post_std):
  fo.write(f'{l[0]:.6}' + " ")
  fo.write(f'{l[1]:.6}' + " ")
  fo.write(f'{l[2]:.6}' + " ")
  fo.write(f'{l[3]:.6}' + " ")
  fo.write(f'{l[4]:.6}' + "\n")
# 

print(r, m_pre, m_pre_std, m_post, m_post_std)


#"/home/piotr/tmp/03_07_2023_GA17_Np64e6_nx1_eps0.1_dt0.1var_sd1e2_Onishi_HallDavis/",
#"/home/piotr/tmp/04_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_sd1_Onishi_HallDavis/",
#"/home/piotr/tmp/04_07_2023_GA17_Np1_nx400_eps0.1_dt0.01var_ConstMulti1_Onishi_HallDavis/",
#"/home/piotr/tmp/06_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_Spinup30_Onishi_HallDavis/",
#"/home/piotr/tmp/07_07_2023_GA17_Np64e3_nx1e1_eps0.1_dt0.1var_SdConc1Tail_DomainInit_Spinup30_Onishi_HallDavis/",
