/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <string>
using namespace std;
//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie, iw, igofr;
double vtail,ptail,bin_size;
const int nbins = 100;
double walker[nbins];

double pot, kin, tempe, etot;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblk;
double delta;
double d_V;

//Blocks
//int NTot, TrBlk, NBlk;
//double* ave_epot; 
//double* ave_ekin; 
//double* ave_etot; 
//double* ave_temp; 

// averages
double blk_av[nbins],blk_norm,accepted,attempted;
double glob_av[nbins],glob_av2[nbins];
double stima_pot,stima_pres,err_pot,err_press,err_gdir, stima_gofr, err_gofr;

// restart
int restart;
double change, fs, temp_att;
//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(int);
double Force(int, int);
double Pbc(double);
//void Blocks(double *, string );
//void PrintBlocks();
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
