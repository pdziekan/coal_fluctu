/*
 * calculate coalescence fluctuations and statistics
 * ensamble of nx cells * n_rep repetitions
 * statistics of max_rw can be calcul;ated using larger cells, conataining more than one cimulation cell
*/

#include <iostream>
#include <fstream>
#include <libcloudph++/lgrngn/factory.hpp>
#include <boost/assign/ptr_map_inserter.hpp>
#include <stdio.h>
#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/vterm.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>
#include <time.h>
#include <libcloudph++/common/earth.hpp>
#include <numeric>
#include <future>
#include <chrono>


#define Onishi
//#define Wang

//#define sgs_ST
#define sgs_GA17

//#define variable_dt

#define cutoff 40e-6
#define HallDavis
#define HIST_BINS 21
#define BACKEND CUDA
#define N_SD_MAX 27e6
#define NXNYNZ 300 //300 //720 // number of cells in each direction
#define SEDI 1
#define RCYC 0
#define N_REP 1e0
#define SIMTIME 800 // [s]
#define NP 1e0 // init number of droplets per cell
#define DT 0.1 // [s]
#define DISS_RATE 1 // [cm^2 / s^3]
#define LKOL 1e-3 // Kolmogorov length scale[m]. Smallest synthetic eddies are od this size 
#define NModes 20 // number of synthethic turbulence modes.
#define NWaves 50 // (max)number of wave vectors for each synthethic turbulence mode.
#define MaxCourant 1 // dt will be adjusted to keep courants less than this
#define OUTFREQ 800 // output done every SIMTIME / OUTFREQ seconds
#define MAXRINTERVAL 0.1 // maximum r is diagnosed every MAXRINTERVAL seconds
#define REMOVE_R 250 // [um] droplets larger than this will be removed 


#if defined sgs_ST && defined sgs_GA17
  #error Both sgs_ST and sgs_GA17 defined
#endif

#if defined Onishi && defined Wang
  #error Both Wang and Onishi defined
#endif

using namespace std;
using namespace libcloudphxx::lgrngn;

using real_t = float;

namespace hydrostatic = libcloudphxx::common::hydrostatic;
namespace theta_std = libcloudphxx::common::theta_std;
namespace theta_dry = libcloudphxx::common::theta_dry;
namespace lognormal = libcloudphxx::common::lognormal;

#ifdef Onishi
const quantity<si::length, real_t>
  mean_rd1 = real_t(15e-6) * si::metres;  
const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
  n1_stp = real_t(142e6) / si::cubic_metres; 
#endif

#ifdef Wang
const quantity<si::length, real_t>
  mean_rd1 = real_t(9.3e-6) * si::metres;  // WANG 2007 (and Unterstrasser 2017)
const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
  n1_stp = real_t(297e6) / si::cubic_metres; // WANG 2007 (and Unter 2017)
#endif

//  mean_rd1 = real_t(0.02e-6) * si::metres;  // api_lgrngn
//  mean_rd1 = real_t(10.177e-6) * si::metres;  // Shima small
const quantity<si::dimensionless, real_t>
  sdev_rd1 = real_t(1.4);
//  n1_stp = real_t(60e6) / si::cubic_metres; // api_lgrngn
//  n1_stp = real_t(226.49e6) / si::cubic_metres; // Shima small


//globals
std::array<real_t, HIST_BINS> rad_bins;
const int n_rep = N_REP; // number of repetitions of simulation
const int nx = NXNYNZ; // total number of collision cells
const int ny = NXNYNZ;
const int nz = NXNYNZ;
constexpr int n_cell = NXNYNZ * NXNYNZ * NXNYNZ;
constexpr int n_courant = (NXNYNZ+1) * NXNYNZ * NXNYNZ;

constexpr real_t Np = NP; // number of droplets per simulation (collision cell)
constexpr real_t Np_in_avg_r_max_cell = Np; // number of droplets per large cells in which we look for r_max
//#ifdef Onishi
//  const int n_cells_per_avg_r_max_cell = Np_in_avg_r_max_cell / Np;
  constexpr real_t cell_vol = Np /  (n1_stp * si::cubic_metres); // for Onishi comparison
  constexpr real_t dx = pow(cell_vol, real_t(1./3.));
  constexpr real_t dy = pow(cell_vol, real_t(1./3.));
  constexpr real_t dz = pow(cell_vol, real_t(1./3.));
//#else
//  const int n_cells_per_avg_r_max_cell = 1; // r_max in each small cell separately
//  const real_t dx = 1e-6; // for bi-disperse (alfonso) comparison
//  const real_t dx = 1e6; // for Shima comparison
//#endif
//const int n_large_cells = (nx * ny * nz) / n_cells_per_avg_r_max_cell;
const int sstp_coal = 1;

// sizes of smallest and largest eddies from the synthetic turbulence scheme
const real_t Lmin = std::max(dx, real_t(LKOL));
const real_t Lmax = nx * dx;

const real_t theta_val = 300;
const real_t rv_val = 0.01;

const auto rho_stp_f = (libcloudphxx::common::earth::rho_stp<real_t>() / si::kilograms * si::cubic_metres);
const auto temperature = libcloudphxx::common::theta_dry::T<real_t>(theta_val * si::kelvins, rho_stp_f * si::kilograms / si::cubic_meters);
const auto pressure = libcloudphxx::common::theta_dry::p<real_t>(rho_stp_f * si::kilograms / si::cubic_meters, rv_val, temperature);
const auto visc = libcloudphxx::common::vterm::visc(temperature);

const real_t outinterval = SIMTIME / OUTFREQ;

const int sd_const_multi = 1; const real_t sd_conc = 0; const bool tail = 0;
//  const int sd_const_multi = 0; const real_t sd_conc = 1e3; const bool tail = 1;

// timing stuff
std::chrono::system_clock::time_point tbeg, tend;
std::chrono::milliseconds tloop, tinit, tadjust, tsync, tasync, trmr, tdiagmax, tdiag;

// lognormal aerosol distribution
template <typename T>
struct log_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    return T(( 
        lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd))
      // +  lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
      ) * si::cubic_metres
    );  
  }   

  log_dry_radii *do_clone() const 
  { return new log_dry_radii( *this ); }
};  

// aerosol distribution exponential in droplet volume as a function of ln(r)
template <typename T>
struct exp_dry_radii : public libcloudphxx::common::unary_function<T>
{
  T funval(const T lnrd) const
  {   
    T r = exp(lnrd);
#ifdef cutoff
    if(r>= cutoff) return 0.; else 
#endif
return (n1_stp * si::cubic_metres) * 3. * pow(r,3) / pow(mean_rd1 / si::metres, 3) * exp( - pow(r/(mean_rd1 / si::metres), 3));
  }   

  exp_dry_radii *do_clone() const 
  { return new exp_dry_radii( *this ); }
};  

void diag(particles_proto_t<real_t> *prtcls, std::array<real_t, HIST_BINS> &res_bins, std::array<real_t, HIST_BINS> &res_stddev_bins)
{
  prtcls->diag_all();
  prtcls->diag_sd_conc();
  std::cout << "sd conc: " << prtcls->outbuf()[0] << std::endl;

/*
  prtcls->diag_all();
  prtcls->diag_dry_mom(0);
  drop_no = prtcls->outbuf()[0];
  std::cout << "number of droplets: " << drop_no << std::endl;

  prtcls->diag_all();
  prtcls->diag_wet_mom(1);
  std::cout << "mean wet radius: " << prtcls->outbuf()[0] / drop_no << std::endl;

  prtcls->diag_all();
  prtcls->diag_dry_mom(1);
  std::cout << "mean dry radius: " << prtcls->outbuf()[0] / drop_no << std::endl;

*/
  prtcls->diag_all();
  prtcls->diag_wet_mom(3);
  real_t sum = 0;
  auto out = prtcls->outbuf();
  #pragma omp parallel for reduction(+ : sum)
  for(int c=0; c < n_cell; ++c)
    sum += out[c];
  std::cout << "3rd wet mom mean: " << sum / n_cell << std::endl;

  real_t r_max_poss = pow(real_t(sum * rho_stp_f * cell_vol), real_t(1./3.));
  real_t vt_max_poss = libcloudphxx::common::vterm::vt_beard76(r_max_poss * si::meters, temperature, pressure, rho_stp_f * si::kilograms / si::cubic_meters, visc) * si::seconds / si::meters;

  std::cout << "max possible rad (based on mean 3rd wet mom): " << r_max_poss * 1e6 << " [um]" << std::endl;
  std::cout << "max possible sedimentation velocity [cm/s] (based on mean 3rd wet mom): " << vt_max_poss*1e2 << std::endl;

  // get spectrum
//  #pragma omp parallel for
  for (int i=0; i <rad_bins.size() -1; ++i)
  {
//    prtcls->diag_all();
  //  prtcls->diag_wet_size_spectr( rad_bins[i], 0.62 ); //sigma0 = 0.62 like in Shima (2009)
    prtcls->diag_wet_rng(rad_bins[i], rad_bins[i+1]);
    prtcls->diag_wet_mom(0);
    real_t rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    real_t mean = 0;
    auto buf = prtcls->outbuf();
    for(int c=0; c < n_cell; ++c)
    {
//      std::cout << buf[c] << " ";
      mean += buf[c];
    }
    mean /= n_cell;

    real_t std_dev = 0;
    for(int c=0; c < n_cell; ++c)
    {
      std_dev += pow(buf[c] - mean, 2.);
    }
    std_dev = sqrt(std_dev / n_cell);

    mean = mean * rho_stp_f; // mean number of droplets of radius rad [1/m^3]
    std_dev *= rho_stp_f;
    
    // to get mass in bins in [g/cm^3]
    
    res_bins[i]= mean / 1e6 // now its number per cm^3
                     * 3.14 *4. / 3. *rad * rad * rad * 1e3 * 1e3;  // to get mass in grams
    res_stddev_bins[i]= std_dev / 1e6 // now its number per cm^3
                     * 3.14 *4. / 3. *rad * rad * rad * 1e3 * 1e3;  // to get mass in grams
    

    // to get mass density function (not through libcloudphxx estimator)
    /*
    res_bins[i] = mean / (rad_bins[i+1] - rad_bins[i]) // number density 
                    * 4./3.*3.14*pow(rad,4)*1e3        // * vol * rad * density
                    * 1e3;                             // to get grams
    */
  }
    std::cout << "res_bins sum (LWC?): " << std::accumulate(res_bins.begin(), res_bins.end(), 0.) << std::endl;
}


int main(int argc, char *argv[]){
  if(argc != 2) throw std::runtime_error("Please specify output prefix");
//  std::cerr << "main start" << std::endl;

  // sanity check
//#ifdef Onishi
//  if(n_cells_per_avg_r_max_cell * Np != Np_in_avg_r_max_cell)
//    throw std::runtime_error("Np_in_avg_r_max_cell nie jest wilokrotnoscia Np");
//  if(n_large_cells * n_cells_per_avg_r_max_cell != n_cell)
//    throw std::runtime_error("n_cell nie jest wilokrotnoscia n_cells_per_avg_r_max_cell");
//#endif

  std::string outprefix(argv[1]);

  std::ofstream of_size_spectr(outprefix+"size_spectr.dat");
  std::ofstream of_series(outprefix+"series.dat");
  std::ofstream of_tau(outprefix+"tau.dat");
  std::ofstream of_rmax(outprefix+"rmax.dat");
  std::ofstream of_nrain(outprefix+"nrain.dat");
  std::ofstream of_time(outprefix+"time.dat");
  std::ofstream of_setup(outprefix+"setup.dat");
  std::ofstream of_t10_tot(outprefix+"t10_tot.dat");

#ifdef cutoff
  of_setup << "init distr cutoff at " << cutoff << " microns!" << std::endl;
#endif

#ifdef variable_dt
  of_setup << "variable dt run" << std::endl;
#endif

#ifdef Onishi
  of_setup << "Onishi (expvolume) run!" << std::endl;
#elif defined Wang
  of_setup << "Wang (expvolume) run!" << std::endl;
#endif

#ifdef sgs_ST
  of_setup << "ST_periodic sgs_adve" << std::endl;
#elif defined sgs_GA17
  of_setup << "GA17 sgs_adve" << std::endl;
#endif

#if defined Onishi || defined Wang
  of_setup << "Np = " << Np << std::endl;
  of_setup << "Np per avg cell = " << Np_in_avg_r_max_cell << std::endl;
#else
  of_setup << "Alfonso (bi-disperse) run!" << std::endl;
#endif
  of_setup << "dx = " << dx * 1e2  << "cm (cell vol = " << cell_vol * 1e6 << " cm^3)"<< std::endl;
  of_setup << "x1 = " << dx * nx * 1e2  << "cm (domain vol = "<< dx * nx * dy * ny * dz * nz  << " m^3)" << std::endl;

  of_setup << "n_rep = " << n_rep 
//            << " n_large_cells = " << n_large_cells
            << " n_cell = " << n_cell
            << " SIMTIME = " << SIMTIME
            << " DT = " << DT
            << " sstp_coal = " << sstp_coal
            << " const_multi = " << sd_const_multi
            << " sd_conc = " << sd_conc
            << " tail = " << tail
            << " mean_rd1 = " << mean_rd1
            << " n1_stp = " << n1_stp
            << " sedi = " << SEDI
            << " rcyc = " << RCYC
            << " backend = " << BACKEND
            << " NModes = " << NModes
            << " NWaves = " << NWaves
            << " Lmin = " << Lmin
            << " Lmax = " << Lmax
            << " MaxCourant = " << MaxCourant
            << " OUTFREQ = " << OUTFREQ
            << " REMOVE_R = " << REMOVE_R
#ifdef HallDavis
            << " kernel: Hall & Davis"
#else
            << " kernel: Hall"
#endif
            << std::endl;

  of_setup << std::flush;


  std::vector<real_t> init_cloud_mass(n_cell);
  std::vector<real_t> init_rain_mass(n_cell);
  real_t init_tot_cloud_mass;
  real_t init_tot_rain_mass;

  std::vector<real_t> t10_tot(n_rep, 0);

  std::vector<std::array<real_t, HIST_BINS>> res_bins_pre(n_rep);
  std::vector<std::array<real_t, HIST_BINS>> res_stddev_bins_pre(n_rep);
  std::vector<std::array<real_t, HIST_BINS>> res_bins_post(n_rep);
  std::vector<std::array<real_t, HIST_BINS>> res_stddev_bins_post(n_rep);
  std::iota(rad_bins.begin(), rad_bins.end(), 0);
  #pragma omp parallel for
  for (auto &rad_bin : rad_bins)
  {
    rad_bin = rad_bin * 1e-6;// + 10e-6; 
  }

  // repetitions loop
  for(int rep = 0; rep < n_rep; ++rep)
  {
    opts_init_t<real_t> opts_init;
  
    opts_init.dt=DT;
    opts_init.sstp_coal = sstp_coal; 
    opts_init.sstp_cond = 1; 
#ifdef variable_dt
    opts_init.variable_dt_switch = true;
#endif
//    opts_init.kernel = kernel_t::hall_pinsky_1000mb_grav;
#ifdef HallDavis
    opts_init.kernel = kernel_t::hall_davis_no_waals;
#else
    opts_init.kernel = kernel_t::hall;
#endif

    opts_init.terminal_velocity = vt_t::beard76;
    opts_init.dx = dx;
    opts_init.dy = dy;
    opts_init.dz = dz; 
  //  cell_vol = opts_init.dx * opts_init.dy * opts_init.dz;
  
    opts_init.sedi_switch=1;
    opts_init.src_switch=0;
    opts_init.chem_switch=0;

#ifdef sgs_ST
    opts_init.sgs_adve=sgs_adve_t::ST_periodic;
    opts_init.ST_Lmax = Lmax;
    opts_init.ST_Lmin = Lmin;
    opts_init.ST_Nmodes = NModes;
    opts_init.ST_Nwaves_max = NWaves;
    opts_init.ST_eps = DISS_RATE * 1e-4;
#endif

#ifdef sgs_GA17
    opts_init.sgs_adve=sgs_adve_t::GA17;
    opts_init.SGS_mix_len = std::vector<real_t>(nz, Lmax);
#endif

//    opts_init.adve_scheme = as_t::pred_corr;

    opts_init.periodic_topbot_walls = 1;
  
    opts_init.nx = nx; 
    opts_init.ny = ny; 
    opts_init.nz = nz; 
    opts_init.x0 = 0;
    opts_init.y0 = 0;
    opts_init.z0 = 0;
    opts_init.x1 = opts_init.nx * opts_init.dx;
    opts_init.y1 = opts_init.ny * opts_init.dy;
    opts_init.z1 = opts_init.nz * opts_init.dz;
    opts_init.rng_seed = time(NULL);
  
    opts_init.sd_conc = sd_conc;//int(1024);
    opts_init.sd_conc_large_tail = tail;
    opts_init.sd_const_multi = sd_const_multi;
  //  opts_init.n_sd_max = 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
    opts_init.n_sd_max = N_SD_MAX;// 20e6 * opts_init.x1 * opts_init.y1 * opts_init.z1 + 1;
//  std::cout << "opts_init.n_sd_max: " << opts_init.n_sd_max << std::endl; 
  
#if defined Onishi || defined Wang
    opts_init.dry_distros.emplace(
      0, // key (kappa)
      std::make_shared<exp_dry_radii<real_t>> () // value
    );
#else
    opts_init.dry_sizes.emplace(
      0, //0 // key (kappa)
      std::map<real_t, std::pair<real_t, int> > {
        {17e-6  , {20e6, 20e6 * cell_vol}}, // radius, STP concentration, number of SD
        {21.4e-6, {10e6, 10e6 * cell_vol}}, // radius, STP concentration, number of SD
      }
    );
#endif

    std::unique_ptr<particles_proto_t<real_t>> prtcls(
      factory<real_t>(
        (backend_t)BACKEND, 
        opts_init
      )
    );
  
  
    std::vector<real_t> pth(n_cell, theta_val); // 300
    std::vector<real_t> prhod(n_cell, rho_stp_f);
    std::vector<real_t> prv(n_cell, rv_val); // .01
  
    const long int strides[] = {1 * NXNYNZ * NXNYNZ, 1 * NXNYNZ, 1};

    arrinfo_t<real_t> th(pth.data(), strides);
    arrinfo_t<real_t> rhod(prhod.data(), strides);
    arrinfo_t<real_t> rv(prv.data(), strides);

#ifdef sgs_GA17
    std::vector<real_t> pdiss_rate(n_cell, DISS_RATE * 1e-4); // 1e-4 to turn cm^2/s^3 to m^2/s^3
    arrinfo_t<real_t> diss_rate(pdiss_rate.data(), strides);
#endif

    prtcls->init(th,rv,rhod, arrinfo_t<real_t>());
//    std::cerr << "post init" << std::endl;
  
    opts_t<real_t> opts;
    opts.adve = 1;
    opts.sedi = SEDI;
    opts.cond = 0;
    opts.coal = 1;
    opts.rcyc = RCYC;
    opts.sgs_adve = 1;
  
    std::fill(res_bins_pre[rep].begin(), res_bins_pre[rep].end(), 0.);
    std::fill(res_bins_post[rep].begin(), res_bins_post[rep].end(), 0.);
    std::fill(res_stddev_bins_pre[rep].begin(), res_stddev_bins_pre[rep].end(), 0.);
    std::fill(res_stddev_bins_post[rep].begin(), res_stddev_bins_post[rep].end(), 0.);
    diag(prtcls.get(), res_bins_pre[rep], res_stddev_bins_pre[rep]);
  
  //  prtcls->step_sync(opts,th,rv);//,rhod);
  //  cout << prtcls->step_async(opts) << endl;

    init_tot_cloud_mass = 0;
    init_tot_rain_mass = 0;
  
    prtcls->diag_wet_rng(0, 40e-6); // cloud water (like in Onishi)
    prtcls->diag_wet_mom(3);
    auto arr = prtcls->outbuf();
    #pragma omp parallel for
    for(int j=0; j<n_cell; ++j)
    {
      init_cloud_mass[j] = arr[j];
    }
    #pragma omp parallel for reduction(+ : init_tot_cloud_mass)
    for(int j=0; j<n_cell; ++j)
    {
      init_tot_cloud_mass += arr[j];
    }
  
    prtcls->diag_wet_rng(40e-6, 1); // rain water (like in Onishi)
    prtcls->diag_wet_mom(3);
    arr = prtcls->outbuf();
    #pragma omp parallel for
    for(int j=0; j<n_cell; ++j)
    {
      init_rain_mass[j] = arr[j];
    }
    #pragma omp parallel for reduction(+ : init_tot_rain_mass)
    for(int j=0; j<n_cell; ++j)
    {
      init_tot_rain_mass += arr[j];
    }
  
    real_t rep_max_rw = 0.;
    // get max rw
    prtcls->diag_max_rw();
    arr = prtcls->outbuf();
    int large_cell_idx = -1;

    // max sgs velocity
    real_t max_sgs_vel = 0;
    real_t Cmax_sgs = 0; // it's courant
    // terminal velocity of the largest droplet
    real_t Cmax_vt = 0;
    real_t outtime = outinterval;
    real_t maxrtime = MAXRINTERVAL;
  
    // output init state
    of_time << "0 ";
    of_tau << "0 ";
    of_nrain << "0 ";
    prtcls->diag_max_rw();
    arr = prtcls->outbuf();
    of_rmax << *(std::max_element(arr, arr+n_cell)) << " ";

    real_t time = 0;
#ifdef variable_dt
    opts.dt = DT;
#endif

    std::cerr << "init cloud_mass_tot: " << init_tot_cloud_mass << std::endl;
    std::cerr << "init rain_mass_tot: " << init_tot_rain_mass << std::endl;
    
    // simulation loop
    while(time <= SIMTIME)
    {
      tbeg = std::chrono::system_clock::now();

#ifdef variable_dt
      if(time>0) // adjust dt
      {
        auto up_minmax = prtcls->diag_up_minmax();
        auto vp_minmax = prtcls->diag_vp_minmax();
        auto wp_minmax = prtcls->diag_wp_minmax();
        max_sgs_vel = std::max(
          max_sgs_vel,
          std::max(std::abs(up_minmax.first), std::abs(up_minmax.second))
        );
        max_sgs_vel = std::max(
          max_sgs_vel,
          std::max(std::abs(vp_minmax.first), std::abs(vp_minmax.second))
        );
        max_sgs_vel = std::max(
          max_sgs_vel,
          std::max(std::abs(wp_minmax.first), std::abs(wp_minmax.second))
        );
        Cmax_sgs = max_sgs_vel * opts.dt / opts_init.dx;

        // Courant number ofr the terminal velocity of the largest droplet
        Cmax_vt = (libcloudphxx::common::vterm::vt_beard76(rep_max_rw * si::meters, temperature, pressure, *(rhod.data) * si::kilograms / si::cubic_meters, visc) * si::seconds / si::meters) * opts.dt / opts_init.dx;

        const real_t max_Cmax = Cmax_sgs + Cmax_vt;// std::max(Cmax, Cmax_vt);
        if(max_Cmax > MaxCourant)
          opts.dt *= MaxCourant / max_Cmax;
      }
#endif
      tend = std::chrono::system_clock::now();
      tadjust += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      tbeg = tend;


#ifdef sgs_ST
      prtcls->step_sync(opts,arrinfo_t<real_t>{},arrinfo_t<real_t>{});
#endif
#ifdef sgs_GA17
      prtcls->step_sync(opts,arrinfo_t<real_t>{},arrinfo_t<real_t>{},arrinfo_t<real_t>{},arrinfo_t<real_t>{},arrinfo_t<real_t>{},arrinfo_t<real_t>{},diss_rate);
#endif

      tend = std::chrono::system_clock::now();
      tsync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      tbeg = tend;

      prtcls->step_async(opts);
      tend = std::chrono::system_clock::now();
      tasync += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      tbeg = tend;

      prtcls->remove_wet_rng(REMOVE_R*1e-6, 1);
      tend = std::chrono::system_clock::now();
      trmr += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );

#ifdef variable_dt
      time += opts.dt;
#else
      time += opts_init.dt;
#endif

      // --- diagnostics ---

      tbeg = std::chrono::system_clock::now();
      if(time > maxrtime)
      {
        maxrtime += MAXRINTERVAL;
        // get max rw
        prtcls->diag_max_rw();
        arr = prtcls->outbuf();
  //      int large_cell_idx = -1;
        #pragma omp parallel for reduction(max : rep_max_rw)
        for(int j=0; j<n_cell; ++j)
          rep_max_rw = rep_max_rw > arr[j] ? rep_max_rw : arr[j];
      }
      tend = std::chrono::system_clock::now();
      tdiagmax += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
      tbeg = tend;

      if(time > outtime)
      {
        outtime += outinterval;

        // sd conc
        prtcls->diag_all();
        prtcls->diag_sd_conc();
        arr = prtcls->outbuf();
        real_t mean_sd_conc = 0;
        #pragma omp parallel for reduction(+ : mean_sd_conc)
        for(int j=0; j<n_cell; ++j)
        {
          mean_sd_conc += arr[j]; 
        }
        mean_sd_conc /= real_t(n_cell);

        // cloud mass
        prtcls->diag_wet_rng(0, 40e-6); // cloud water (like in Onishi)
        prtcls->diag_wet_mom(3);
        arr = prtcls->outbuf();
        real_t cloud_mass_tot = 0;
        #pragma omp parallel for reduction(+ : cloud_mass_tot)
        for(int j=0; j<n_cell; ++j)
          cloud_mass_tot += arr[j];

        // rain mass
        prtcls->diag_wet_rng(40e-6, 1); // rain water (like in Onishi)
        prtcls->diag_wet_mom(3);
        arr = prtcls->outbuf();
        real_t rain_mass_tot = 0;
        #pragma omp parallel for reduction(+ : rain_mass_tot)
        for(int j=0; j<n_cell; ++j)
          rain_mass_tot += arr[j];
        // add removed particles (due to large r) to rain water mass
        rain_mass_tot += init_tot_cloud_mass - (cloud_mass_tot + rain_mass_tot);

        // rain conc
        prtcls->diag_wet_mom(0);
        arr = prtcls->outbuf();
        real_t nrain_mean = 0; // concentration of rain droplets, averaged over all cells (whole domain)
        #pragma omp parallel for reduction(+ : nrain_mean)
        for(int j=0; j<n_cell; ++j)
        {
          nrain_mean += arr[j];
        }
        nrain_mean = nrain_mean * rho_stp_f / n_cell; // [1/m^3]

        // t 10%
        if(t10_tot[rep] == 0. && rain_mass_tot >= init_tot_cloud_mass * .1)
        {
          t10_tot[rep] = time;
          of_t10_tot << t10_tot[rep] << std::endl;
        }
    
        // output
        printf("rep no: %3d progress: %3d%%: dt: %lf sstp_coal: %3d rw_max %lf [um] mean_sd_conc %lf max_sedi_courant %lf max_sgs_courant %lf t10_tot %lf\n", rep, int(time / SIMTIME * 100), opts.dt, prtcls->diag_sstp_coal(), rep_max_rw * 1e6, mean_sd_conc, Cmax_vt, Cmax_sgs, t10_tot[rep]);
        of_rmax << rep_max_rw << " ";
        of_tau << rain_mass_tot / init_tot_cloud_mass << " ";
        of_time << time << " ";
        of_nrain << nrain_mean << " ";

//        std::cerr << "cloud_mass_tot: " << cloud_mass_tot << std::endl;
//        std::cerr << "rain_mass_tot: " << rain_mass_tot << std::endl;

        std::cout << std::flush;
        of_rmax << std::flush;
        of_tau << std::flush;
        of_nrain << std::flush;
        of_time << std::flush;
      }
      tend = std::chrono::system_clock::now();
      tdiag += std::chrono::duration_cast<std::chrono::milliseconds>( tend - tbeg );
    }
    of_rmax << std::endl;
    of_tau << std::endl;
    of_nrain << std::endl;
    of_time << std::endl;
  
    diag(prtcls.get(), res_bins_post[rep], res_stddev_bins_post[rep]);
    std::cout << std::endl;

    auto raw_ptr = prtcls.release();
    delete raw_ptr;

  }

  std::cout <<  "wall time in milliseconds: " << std::endl
    << "adjust dt:      " << tadjust.count() << std::endl
    << "sync:           " << tsync.count() << std::endl
    << "async:          " << tasync.count() << std::endl
    << "remove large r: " << trmr.count() << std::endl
    << "diagmax rw:     " << tdiagmax.count() << std::endl
    << "diag:           " << tdiag.count() << std::endl;


  real_t mean_t10_tot=0, std_dev_t10_tot=0;
  for(int j=0; j<n_rep; ++j)
    mean_t10_tot += t10_tot[j];
  mean_t10_tot /= n_rep;
  for(int j=0; j<n_rep; ++j)
    std_dev_t10_tot += pow(t10_tot[j] - mean_t10_tot,2);
  std_dev_t10_tot = sqrt(std_dev_t10_tot / (n_rep-1));
  std::cout << "mean(t10% in the domain) = " << mean_t10_tot << std::endl;
  std::cout << "std_dev(t10% in the domain) = " << std_dev_t10_tot << std::endl;

// output
  for (int i=0; i <rad_bins.size() -1; ++i)
  {
    real_t rad = (rad_bins[i] + rad_bins[i+1]) / 2.;
    real_t pre = 0, post = 0, stddev_post = 0;
    for(int j=0; j< n_rep; ++j)
    {
      pre += res_bins_pre[j][i];
      post += res_bins_post[j][i];
      stddev_post += res_stddev_bins_post[j][i];
    }
    pre /= n_rep;
    post /= n_rep;
    stddev_post /= n_rep;
    of_size_spectr << rad * 1e6 << " " << pre << " " << post << " " << stddev_post << std::endl; 
  }
}
