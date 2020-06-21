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

//--------------Uniform-exp-lorentz dices (Central Limit Theorem)-----------------

	ofstream out("Es01.2.1.res"); //uniform dice
	ofstream out2("Es01.2.2.res"); //exponential dice
	ofstream out3("Es01.2.3.res"); //lorentzian dice
	int R=1E4;
	int N[4]={1, 2, 10, 100};
	double S[4][R]; //uniform dice
	double E[4][R]; //exponential dice
	double L[4][R]; //lorentzian dice	
	
	double l=1; //lambda (exp)
	double g=1; //gamma (lorentz)
	double mu=0; //mu (lorentz)
	
	
	for(int j=0; j<R; j++) {
		for(int k=0; k<4; k++){
			S[k][j]=0;
			E[k][j]=0;
			L[k][j]=0;
		}
	}
	
	for(int k=0; k<4; k++){
     for(int i=0; i<N[k]; i++){
       for(int j=0; j<R; j++){
       	S[k][j]+=rnd.Rannyu()/(double)(N[k]);
       	//E[k][j]+=(double)(-log(1-rnd.Rannyu())/l)/(double)(N[k]);
       	//L[k][j]+=(double)(g*tan(M_PI*(rnd.Rannyu()-0.5)))/(double)(N[k]);
       	E[k][j]+=rnd.Exp(l)/(double)(N[k]);
       	L[k][j]+=rnd.Lorentz(g,mu)/(double)(N[k]);
       	}
     }
   }
   out<<setw(20);
   out2<<setw(20);
   out3<<setw(20);
	for(int j=0; j<R; j++){
		for(int k=0; k<4; k++){
			out<<S[k][j]<<setw(20);
			out2<<E[k][j]<<setw(20);
			out3<<L[k][j]<<setw(20);
		}
		out<<endl;
		out2<<endl;
		out3<<endl;
	}
	
	out.close();
	out2.close();
	out3.close();

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
