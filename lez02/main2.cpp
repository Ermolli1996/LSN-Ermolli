/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../lez00/ParRandGen/random.h"
#include "sigma.h"

using namespace std;
 
int main (int argc, char *argv[]){
//---------RANDOM GENERETOR SETTING--------------
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("../lez00/ParRandGen/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../lez00/ParRandGen/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
//-------------PART 1 (RW 3D lattice)------------------
	int M=1E4;
	int N=100;
	int L=int(M/N);
	int nsteps=100;	//Steps of RW
   ofstream out("Es02.2.1.res");
   double ave[N][nsteps];
   double av2[N][nsteps];
   double a=1.;
   double sum[nsteps];
	
	for(int i=0; i<N; i++){
   	for(int k=0; k<nsteps; k++){
     		sum[k]=0;
      	ave[i][k]=0;
      	av2[i][k]=0;
      }
     	for(int j=0; j<L; j++){
     		double r[3]={0,0,0};
			for(int k=0; k<nsteps; k++){
				int z=rnd.Rannyu(0.,3.);
				double q=rnd.Rannyu(0.,1.);
				int sign=1;
				if(q<0.5){
					sign=-1;
				}
				r[z]+=a*sign;
				sum[k]+=pow(r[0],2.)+pow(r[1],2.)+pow(r[2],2.);
			}
		}
		for(int k=0; k<nsteps; k++){
      	ave[i][k]=(double)sum[k]/(double)(L);
      	av2[i][k]=pow(ave[i][k],2);
      }
   }
	
	double sum_prog[N][nsteps];
	double su2_prog[N][nsteps];
	double err_prog[N][nsteps];
	  
	for(int k=0; k<nsteps; k++){
		for (int i=0; i<N; i++){
			sum_prog[i][k]=0;
			su2_prog[i][k]=0;	
			err_prog[i][k]=0;
		for (int j=0; j<i+1; j++){
			sum_prog[i][k]+= ave[j][k];
			su2_prog[i][k]+= av2[j][k];
		}
		sum_prog[i][k]/=(double)(i+1);
		su2_prog[i][k]/=(double)(i+1);
		//err_prog[i][k] = error(sum_prog[k],su2_prog[k],i);
		if(i!=0){
      	err_prog[i][k]=sqrt((double)(su2_prog[i][k]-pow(sum_prog[i][k],2))/i);
      }
     	else{
     		err_prog[0][k]=0;
     	}
     	if(i==N-1){  
			out<<k+1<<setw(20)<<sqrt(sum_prog[i][k])<<setw(20)<<err_prog[i][k]/(2*sqrt(sum_prog[i][k]))<<endl;  
//stampo la radice! E per l'errore propago l'errore
		}
		}	
	}	
			
	out.close();

//-------------PART 2 (RW in continuum)------------------
	ofstream out2("Es02.2.2.res");
	for(int i=0; i<N; i++){
   	for(int k=0; k<nsteps; k++){
     		sum[k]=0;
      	ave[i][k]=0;
      	av2[i][k]=0;
      }
     	for(int j=0; j<L; j++){
     		double r[3]={0,0,0};
			for(int k=0; k<nsteps; k++){
				double phi=rnd.Rannyu(0.,2*M_PI);
				double theta=acos(1-2*rnd.Rannyu(0.,1.)); //p(theta)=1/2 * sin(theta)
				r[0]+=a*sin(theta)*cos(phi);
				r[1]+=a*sin(theta)*sin(phi);
				r[2]+=a*cos(theta);
				sum[k]+=pow(r[0],2.)+pow(r[1],2.)+pow(r[2],2.);
			}
		}
		for(int k=0; k<nsteps; k++){
      	ave[i][k]=(double)sum[k]/(double)(L);
      	av2[i][k]=pow(ave[i][k],2);
      }
   }
	  
	for(int k=0; k<nsteps; k++){
		for (int i=0; i<N; i++){
			sum_prog[i][k]=0;
			su2_prog[i][k]=0;	
			err_prog[i][k]=0;
		for (int j=0; j<i+1; j++){
			sum_prog[i][k]+= ave[j][k];
			su2_prog[i][k]+= av2[j][k];
		}
		sum_prog[i][k]/=(double)(i+1);
		su2_prog[i][k]/=(double)(i+1);
		//err_prog[i][k] = error(sum_prog[k],su2_prog[k],i);
		if(i!=0){
      	err_prog[i][k]=sqrt((double)(su2_prog[i][k]-pow(sum_prog[i][k],2))/i);
      }
     	else{
     		err_prog[0][k]=0;
     	}
     	if(i==N-1){  
			out2<<k+1<<setw(20)<<sqrt(sum_prog[i][k])<<setw(20)<<err_prog[i][k]/(2*sqrt(sum_prog[i][k]))<<endl;  
//stampo la radice! E per l'errore propago l'errore
		}
		}	
	}	
	
	out2.close();	

   rnd.SaveSeed();

   return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
