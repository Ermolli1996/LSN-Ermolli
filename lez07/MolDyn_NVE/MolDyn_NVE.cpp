/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; iblk++){
    Reset(iblk);
    for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure(istep);     //Properties measurement
        Accumulate();
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
    }
    Averages(iblk);
  }
  ConfFinal();         //Write final configuration to restart
	//PrintBlocks();
	
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> restart;
  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> nblk;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps in one block = " << nstep << endl;
  cout << "Blocks for Data-Blocking = " << nblk << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
	
	
	//g(r)
  bin_size = (box/2.0)/(double)nbins;
  
  
	/*NTot=nstep/10.;
  ave_epot=new double[NTot]; 
  ave_ekin=new double[NTot]; 
  ave_etot=new double[NTot]; 
  ave_temp=new double[NTot]; */

//Read initial configuration
	if(restart==0){
  	cout << "Read initial configuration from file config.0 " << endl << endl;
  	ReadConf.open("config.0");
  	for (int i=0; i<npart; ++i){
    	ReadConf >> x[i] >> y[i] >> z[i];
    	x[i] = x[i] * box;
    	y[i] = y[i] * box;
    	z[i] = z[i] * box;
  	}
  	ReadConf.close();
  	
//Prepare initial velocities
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
     	vx[i] = rand()/double(RAND_MAX) - 0.5;
     	vy[i] = rand()/double(RAND_MAX) - 0.5;
     	vz[i] = rand()/double(RAND_MAX) - 0.5;

     	sumv[0] += vx[i];
     	sumv[1] += vy[i];
     	sumv[2] += vz[i];
   	}
   	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   	double sumv2 = 0.0;
   	for (int i=0; i<npart; ++i){
     	vx[i] = vx[i] - sumv[0];
     	vy[i] = vy[i] - sumv[1];
     	vz[i] = vz[i] - sumv[2];

     	sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	sumv2 /= (double)npart;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   	for (int i=0; i<npart; ++i){
     	vx[i] *= fs;
     	vy[i] *= fs;
     	vz[i] *= fs;

     	xold[i] = Pbc(x[i] - vx[i] * delta);
     	yold[i] = Pbc(y[i] - vy[i] * delta);
     	zold[i] = Pbc(z[i] - vz[i] * delta);
   		}
  }
   //Restart!!
	else if(restart==1){
   	//read final and old
   	cout << "Read initial configuration from file config.final " << endl << endl;
    ReadConf.open("config.final");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
   
    cout << "Read initial old configuration from file old.final " << endl << endl;
    ReadConf.open("old.final");
    for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();
   	//actual T
   	Move();
   	temp_att=0;
		for(int i=0; i<npart; ++i) temp_att+= 0.5* (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  	temp_att=(2.0 / 3.0) * temp_att/(double)npart;
  	//scaling factor
  	fs = sqrt(temp / temp_att); 	
    for (int i=0; i<npart; ++i){	
      vx[i] *= fs;			
      vy[i] *= fs;
      vz[i] *= fs;
    }
    //xold 
    for (int i=0; i<npart; ++i){	
      xold[i] = Pbc( x[i] - delta * vx[i] );
      yold[i] = Pbc( y[i] - delta * vy[i] );
      zold[i] = Pbc( z[i] - delta * vz[i] );
    }
	}	    
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(int istep){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
  
//reset the hystogram of g(r)
  for (int k=0; k<nbins; k++) walker[k]=0.0;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
		 if(dr<box/2.0){
		 	bin=int(dr/bin_size);
     	walker[bin] += 2.0;
     }
               
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
//Potential energy
       v += vij;
     }
    }
         
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    pot = v/(double)npart; //Potential energy per particle
    kin = t/(double)npart; //Kinetic energy per particle
    tempe = (2.0 / 3.0) * t/(double)npart; //Temperature
    etot = (t+v)/(double)npart; //Total energy per particle

    Epot << pot << endl;
    Ekin << kin  << endl;
    Temp << tempe << endl;
    Etot << etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
	
		//ave_epot[istep/10 -1]=stima_pot;
    //ave_ekin[istep/10 -1]=stima_kin;
    //ave_temp[istep/10 -1]=stima_temp;
    //ave_etot[istep/10 -1]=stima_etot;
    return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<nbins; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<nbins; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<nbins; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
    
   double r, gdir;
   ofstream Gave, Epot;
   const int wd=12;
    
   cout << "Block number " << iblk << endl << endl;
   
   Epot.open("output.epot.0",ios::app);
   //Gofr.open("output.gofr.0",ios::app);
   Gave.open("output.gave.0",ios::app);
    
   stima_pot = pot/blk_norm; //Potential energy
   glob_av[iv] += stima_pot;
   glob_av2[iv] += stima_pot*stima_pot;
   err_pot=Error(glob_av[iv],glob_av2[iv],iblk);

//Potential energy per particle
   Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;

//g(r)
		for(int i=0; i<nbins; i++){
      r=(i)*bin_size;
      d_V=4./3.*M_PI*(pow(r+bin_size,3)-pow(r,3));
      stima_gofr = blk_av[i]/blk_norm/(rho*npart*d_V);
      glob_av[i] += stima_gofr;
      glob_av2[i] += stima_gofr*stima_gofr;
      err_gofr = Error(glob_av[i], glob_av2[i], iblk);
      //Gofr << setw(wd) << iblk << setw(wd) << r << setw(wd) << stima_gofr << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gofr << endl;
      //output.gofr.0: (iblk), (r), (stimag_gofr), (<gofr>), (sigma_gofr)
      if(iblk==nblk)
        Gave << setw(wd) << r << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gofr << endl;
    }
      

    cout << "----------------------------" << endl << endl;

    Epot.close();
}



void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  
  //Saving also the n-1 th cofiguration
  ofstream WriteConf_Old;

  cout << "Print final old configuration to file old.final " << endl << endl;
  WriteConf_Old.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf_Old << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf_Old.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/*void Blocks(double *meanval, string outfile){
   TrBlk=int(NTot/NBlk);
   ofstream out;
   out.open(outfile,ios::app);
   double ave;
   double sigma;
   double true_sum=0, true_sum2=0;
   double true_ave, true_ave2;
   int i;
   for(i=0; i<NBlk; i++){
     double sum=0;
     for(int j=0; j<TrBlk; j++)
        sum+=meanval[i*TrBlk+j];		
     ave=(double)sum/(double)TrBlk;
     true_sum+=ave;											
     true_sum2+=ave*ave;
     true_ave=true_sum/(double)(i+1);
     true_ave2=true_sum2/(double)(i+1);
     if(i!=0)
       sigma=sqrt((double)(true_ave2-pow(true_ave,2))/i);
     else
       sigma=0;
     out<<i<<setw(15)<<true_ave<<setw(15)<<sigma<<endl;
   }
   out.close();
   return;
}

void PrintBlocks(){

  string output_ave_epot="output_ave_epot.dat";
  string output_ave_etot="output_ave_etot.dat";
  string output_ave_ekin="output_ave_ekin.dat";
  string output_ave_temp="output_ave_temp.dat";

  cout<<"Printing average potential energy to file output_ave_epot.dat"<<endl;
  Blocks(ave_epot, output_ave_epot);

  cout<<"Printing average total energy to file output_ave_etot.dat"<<endl;
  Blocks(ave_etot, output_ave_etot);

  cout<<"Printing average kinetic energy to file output_ave_ekin.dat"<<endl;
  Blocks(ave_ekin, output_ave_ekin);

  cout<<"Printing average temperature to file output_ave_temp.dat"<<endl;
  Blocks(ave_temp, output_ave_temp);

  delete[] ave_epot;
  delete[] ave_ekin;
  delete[] ave_etot;
  delete[] ave_temp;
  
  return;
}*/
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
