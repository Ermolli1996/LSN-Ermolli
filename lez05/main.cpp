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
#include "main.h"

using namespace std;
 
int main (int argc, char *argv[]){
	
  settingRnd();
	input();
	
	int L=int(M/N);
	double ave[N];
	double av2[N];
	
	string f1,f2,f3;
	if(distribution=='U') f1="_Uni";
	else f1="_Gauss";
	if(state==1) f2="_1s";
	else f2="_2p";
	f3=""; //per simulzioni particolari (scentrata,...)
	ofstream out("../lez05/Risultati/coordinates"+f1+f2+f3+".res");

	int cont=0;
	x0=x_start;
	b0=b_start;
	z0=z_start;
	//Metropolis
	for(int i=0; i<N; i++){
		
		double sum=0;
		for (int j=0; j<L; j++){
			if(distribution=='U'){
				x1=x0+rnd.Rannyu(-0.5,0.5)*passo;
				b1=b0+rnd.Rannyu(-0.5,0.5)*passo;
				z1=z0+rnd.Rannyu(-0.5,0.5)*passo;
			}
			else{
				x1=x0+rnd.Gauss(0,0.5)*passo;
				b1=b0+rnd.Gauss(0,0.5)*passo;
				z1=z0+rnd.Gauss(0,0.5)*passo;
			}
			if(state==1){
				ratio=psi_1s(x1,b1,z1)/psi_1s(x0,b0,z0);
			}
			else{
				ratio=psi_2p(x1,b1,z1)/psi_2p(x0,b0,z0);
			}
			if(rnd.Rannyu() <= min(1.,ratio)){
       	x0=x1;
       	b0=b1;
       	z0=z1;
       	cont++;
     	}
    	sum+=sqrt(x0*x0+b0*b0+z0*z0);
			out<<j+i*L<<setw(20)<<x0<<setw(20)<<b0<<setw(20)<<z0<<endl;
		}
		ave[i]=(double)sum/(double)L;
		av2[i]=pow((ave[i]),2);	
	}
	
	out.close();
	ofstream out2("../lez05/Risultati/rmean"+f1+f2+f3+".res");
	double sum_prog[N];
	double su2_prog[N];
	double err_prog[N];
	
	//Data Blocking
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
	cout<<"Rapporto di accettazione:  "<<(double)cont/M<<endl;
	
	rnd.SaveSeed();
  return 0;
}


void settingRnd(){
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
}   

void input(){
	ifstream in("input.dat");
  if(!in.good()){
     cerr<<"Error reading file"<<endl;
   }
   in>>M;
   in>>N;
   in>>passo;
   in>>x_start;
   in>>b_start;
   in>>z_start;
   in>>state;
   in>>distribution;
   if(state!=1 && state!=2){
     cerr<<"Error in selected state: 1s (state=1) or 2p (state=2)"<<endl;
   }
   if(distribution!='U' && distribution!='G'){
     cerr<<"Error in selected distribution: U (Uniform) and G (Gaussian)"<<endl;
   }
   in.close();
}

double error(double * AV, double * AV2, int n){
	if (n==0){
		return 0;
	}
	else return sqrt((double)(AV2[n] - pow(AV[n],2))/n);
}

double psi_1s(double x,double y,double z){
 	double r=sqrt(x*x+y*y+z*z);
  double psi=exp(-r);
  return psi*psi;
}

double psi_2p(double x,double y,double z){
 	double r=sqrt(x*x+y*y+z*z);
  double psi=exp(-r/2)*z;
  return psi*psi;
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
