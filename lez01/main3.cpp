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

//--------------Buffon's experiment: estimation of pi-----------------

	ofstream out("Es01.3.res"); 
	double d=15.; //distance
	double L=8.; //needle's lenght
	int Nthr=1E7; //trials (lanci del bastonino)
	int Nblk=100;	//blocks
	int Lb=int(Nthr/Nblk); //throws per block
	double pi[Nblk];
	double ave[Nblk];
	double av2[Nblk];
	
		for(int s=0; s<Nblk; s++){
		double Nhit=0;
		for(int i=0; i<Lb; i++){
			double c=rnd.Rannyu(0., d); // needle's center
			double x,y,a;
			do{
				x=rnd.Rannyu(-1., 1.);
				y=rnd.Rannyu(0., 1.);
				a=acos(x/(pow((x*x+y*y),0.5))); //a U[0,pi] no using pi (cos pari)
				} while(x*x+y*y>=1); //(x,y) nella semi-cir unitaria
		
			double y_up=c+(L/2)*sin(a);
			double y_down=c-(L/2)*sin(a);
			if(y_down<=0. || y_up>=d){
				Nhit++;
				}
		}
		pi[s]=(2*L*Lb)/(Nhit*d);
	}

	for (int i=0; i<Nblk; i++){
		ave[i]=pi[i];			//r_i
		av2[i]=pow((ave[i]),2);					//(r_i)^2
	}
	
	double sum_prog[Nblk];
	double su2_prog[Nblk];
	double err_prog[Nblk];
	for (int i=0; i<Nblk; i++){
		sum_prog[i]=0;
		su2_prog[i]=0;	
		for (int j=0; j<i+1; j++){
			sum_prog[i]+= pi[j];
			su2_prog[i]+= av2[j];
		}
		sum_prog[i]/=(double)(i+1);
		su2_prog[i]/=(double)(i+1);
		err_prog[i] = error(sum_prog,su2_prog,i);
		out<<i<<setw(20)<<sum_prog[i]<<setw(20)<<err_prog[i]<<endl;
	}

	out.close();
	
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
