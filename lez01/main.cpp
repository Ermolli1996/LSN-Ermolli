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
   
//-------------PART 1 (<r> as integral)------------------
	int M=100000;
	int N=100;
	int L=int(M/N);

	ofstream out("Es01.1.1.res"); //3 columns: step,result,sigma
	double r[M];
	double ave[N];
	double av2[N];
	
   for(int i=0; i<M; i++){
      r[i]=rnd.Rannyu();
	}
	
	for (int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			int k=j+i*L;
			sum +=r[k];
		}
		ave[i]=(double)sum/(double)L;			//r_i
		av2[i]=pow((ave[i]),2);					//(r_i)^2
	}
	
	double sum_prog[N];
	double su2_prog[N];
	double err_prog[N];
	for (int i=0; i<N; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;	
		for (int j=0; j<i+1; j++){
			sum_prog[i]+= ave[j];
			su2_prog[i]+= av2[j];
		}
		sum_prog[i]/=(double)(i+1);
		su2_prog[i]/=(double)(i+1);
		err_prog[i] = error(sum_prog,su2_prog,i);
		out<<i<<setw(20)<<sum_prog[i]<<setw(20)<<err_prog[i]<<endl;
	}
	out.close();
	
//-------------PART 2 (sigma2 as integral)------------------
	ofstream out2("Es01.1.2.res"); //3 columns: step,result,sigma
	
	for (int i=0; i<N; i++){
		double sum=0;
		ave[i]=0;
     	av2[i]=0;
		for (int j=0; j<L; j++){
			int k=j+i*L;
			sum +=pow((r[k]-0.5),2);
		}
		ave[i]=(double)sum/(double)L;			//r_i
		av2[i]=pow((ave[i]),2);					//(r_i)^2
	}
	
	for (int i=0; i<N; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;	
		err_prog[i]=0;
		for (int j=0; j<i+1; j++){
			sum_prog[i]+= ave[j];
			su2_prog[i]+= av2[j];
		}
		sum_prog[i]/=(double)(i+1);
		su2_prog[i]/=(double)(i+1);
		err_prog[i] = error(sum_prog,su2_prog,i);
		out2<<i<<setw(20)<<sum_prog[i]<<setw(20)<<err_prog[i]<<endl;
	}
	out2.close();
	
//-------------PART 3 (chi2)------------------

	int S=100; 	//number of sub-intervals
	int n=1E4;	//number of pseudo numbers each time
	int t=100;	//times
	double a[n];
	double chi2=0;
	int occ=0;
	ofstream out3("Es01.1.3.res");
	
	for(int i=0; i<t; i++){
		chi2=0;
		for(int j=0; j<n; j++){
			a[j]=rnd.Rannyu();
		}
		for(int w=0; w<S; w++){
			occ=0;
			for(int j=0; j<n; j++){
				if(a[j]>=((double)w/(double)(S))&&a[j]<((double)(w+1)/(double)(S))) occ++;}
			chi2+=(double)pow((occ-n/S),2.)/(double)(n/S);
		}
		out3<<i<<setw(20)<<chi2<<endl;
		}	
		
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
