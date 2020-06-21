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
   
//------------- European call-option and put-option prices------------------
	int M=1E5;
	int N=100;
	int L=int(M/N);
	
	double S0=100;
	double T=1;
	double K=100;
	double r=0.1;
	double sigma=0.25;
	double S;
	
//sampling directly the final asset price  𝑆(𝑇) 

//CALL
	ofstream out("Es03.1.1.res"); //3 columns: step,C,sigma_C,
	double ave[N];
	double av2[N];

   for(int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			S=S0*exp((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0, T));
			sum+=exp(-r*T)*(double)max(0.,S-K);
		}
		ave[i]=(double)sum/(double)L;
		av2[i]=pow((ave[i]),2);	
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

//PUT
	ofstream out2("Es03.1.2.res"); //3 columns: step,P,sigma_P


   for(int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			S=S0*exp((r-0.5*sigma*sigma)*T+sigma*rnd.Gauss(0, T));
			sum+=exp(-r*T)*(double)max(0.,K-S);
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
 
//sampling the discretized  𝐺𝐵𝑀(𝑟,𝜎2)

	int n=100;
   double t=T/(double)(n);
//CALL
	ofstream out3("Es03.1.3.res"); //3 columns: step,C,sigma_C,
	for(int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			S=S0;
     		for(int z=1; z<(n+1); z++){
       		S=S*exp((r-0.5*pow(sigma,2.))*(t*z-t*(z-1))+sigma*rnd.Gauss(0.,1.)*sqrt(t*z-t*(z-1)));
      	}
			sum+=exp(-r*T)*(double)max(0.,S-K);
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
		out3<<i<<setw(20)<<sum_prog[i]<<setw(20)<<err_prog[i]<<endl;
	}
	out3.close();
	
//PUT
	ofstream out4("Es03.1.4.res"); //3 columns: step,P,sigma_P
	for(int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){
			S=S0;
     		for(int z=1; z<(n+1); z++){
       		S=S*exp((r-0.5*pow(sigma,2.))*(t*z-t*(z-1))+sigma*rnd.Gauss(0.,1.)*sqrt(t*z-t*(z-1)));
      	}
			sum+=exp(-r*T)*(double)max(0.,K-S);
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
		out4<<i<<setw(20)<<sum_prog[i]<<setw(20)<<err_prog[i]<<endl;
	}
	out4.close();
	
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
