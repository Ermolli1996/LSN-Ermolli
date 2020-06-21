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
   
//-------------PART 1 (integral sampling a uniform distribution)------------------
	int M=1E6;
	int N=100;
	int L=int(M/N);

	ofstream out("Es02.1.1.res"); //3 columns: step,result,sigma
	double ave[N];
	double av2[N];
	
   for(int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			sum +=(double)(M_PI/2)*(double)cos(M_PI*rnd.Rannyu()/2);
			ave[i]=(double)sum/(double)L;
			av2[i]=pow((ave[i]),2);	
		}
	}
	double sum_prog[N];
	double su2_prog[N];
	double err_prog[N];
	
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
		out<<i<<setw(20)<<sum_prog[i]<<setw(20)<<err_prog[i]<<endl;
	}
	out.close();

//-------------PART 2 (integral with importance sampling)------------------
	ofstream out2("Es02.1.2.res"); //3 columns: step,result,sigma
	//p(x)=2(1-x), z=1-sqrt(1-y), g=g_prec/p 
	
	double z[M];
	
   for(int i=0; i<M; i++){
      z[i]=(double)(1-sqrt(1-rnd.Rannyu()));
	}
	
	for (int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			int k=j+i*L;
			sum +=(double)(M_PI/2)*(double)cos(M_PI*z[k]/2)/(2*(1-z[k]));
		}
		ave[i]=(double)sum/(double)L;
		av2[i]=pow((ave[i]),2);
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
