using namespace std;

//Random numbers
#include "../lez00/ParRandGen/random.h"
Random rnd;
int seed[4];
int p1, p2;

//parameters, observables
const int m_props=1000;
int n_props, iH;
double walker[m_props];
const int nbin=100;
double histo[nbin];
const double xd=3.0, xs=-3.0, binsize=(xd-xs)/(double)nbin;
int norm=0;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima, err;

// simulation
int nstep, nblk;
double passo, mu, sigma, x0, x1;
double ratio;

int att, succ;
int iterations, stepforT, n_moves;
double temp, fract, delta_mu, delta_sigma;
double new_mu, new_sigma, p, xstart;
double E, E_new, old_mu, old_sigma;

//functions
void input(void);
void settingRnd(void);
void SigmaMu(int);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(double, double );
void PrintFinal(void);
void Measure(double, double );
double psi(double, double, double );
double psi2(double, double, double );
double psi_second(double, double, double);
double potential(double );
double Error(double,double,int);


