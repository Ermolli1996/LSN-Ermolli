#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../lez00/ParRandGen/random.h"
#include "main.h"

using namespace std;

int main (){
	
  settingRnd();
	input();
	
	for(int i=0; i<iterations; i++){
		SigmaMu(i);
		temp=temp/fract;
	}
	cout <<endl <<endl;
	cout<<"Data-Blocking"<<endl;
	cout<< nblk << " blocks"<<endl;
	cout<< nstep << " MC steps in each block"<<endl <<endl;
	
	for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(sigma,mu);
      Measure(sigma,mu);
      Accumulate();
    }
    Averages(iblk);   //Print results for current block
  }
  PrintFinal(); //Write final configuration

  return 0;
}

void settingRnd(){
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
  in >> nstep;
  in >> nblk;
  in >> passo;
  in >> mu;
  in >> sigma;
  in >> x0;
  in >> iterations;
  in >> temp;
  in >> fract;
  in >> stepforT;
  in >> delta_mu;
  in >> delta_sigma;
  in >> n_moves;
  in.close();

//Prepare arrays for measurements
  iH = 0; //energy
  n_props = 1; //Number of observables

  for(int i=0; i<nbin; ++i) histo[i]=0.;
	cout << "Varational Monte Carlo: " << endl;
	cout << n_moves <<" moves to evaluate energy at every T " << endl;
 	cout << iterations << " temperatures" << endl;
 	cout << stepforT << " steps for every temperature" << endl << endl;
  return;
}
	
void SigmaMu(int i)
{	
	att=0;
	succ=0;
	for(int istep=0; istep<stepforT; istep++){
   	new_mu=mu + rnd.Rannyu(-0.5,0.5)*delta_mu;
   	new_sigma=sigma + rnd.Rannyu(-0.5,0.5)*delta_sigma;
   	//Uso come funzione costo l'energia del sistema mediata su n_moves passi montecarlo con stesso xstart!
   	E=0; 
   	E_new=0;
   	xstart=x0;
   	for(int k=0; k<n_moves; k++){
   		Move(sigma,mu);	// mu e sigma old
   		E+=(-0.5*psi_second(x0,sigma,mu)+potential(x0)*psi(x0,sigma,mu))/psi(x0,sigma,mu);
   	}
   	E=E/n_moves;
   	x0=xstart;
   	for(int k=0; k<n_moves; k++){
   		Move(new_sigma, new_mu); // mu e sigma new
   		E_new+=(-0.5*psi_second(x0,new_sigma,new_mu)+potential(x0)*psi(x0,new_sigma,new_mu))/psi(x0,new_sigma,new_mu);
   	}
   	E_new=E_new/n_moves;
   	p=min(1., exp(-( E_new - E )/(double)temp));
  	if(rnd.Rannyu()<=p){
  		succ++;
  		mu=new_mu;
  		sigma=new_sigma;
  	}
  	att++;
  }
  if(i%10==0){
		cout << "T = " << temp << endl;
  	cout << "Acceptance rate : " << (double)succ/att << endl;
  }
  return;
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void) //Update block averages
{
   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk)
{
   ofstream out;
  
   cout << "Block number " << iblk << endl;
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
   out.open("output.ene.dat",ios::app);
    
   stima = blk_av[iH]/blk_norm;
   glob_av[iH] += stima;
   glob_av2[iH] += stima*stima;
   err=Error(glob_av[iH],glob_av2[iH],iblk);

//Potential energy per particle
    out << setw(15) << iblk <<  setw(15) << stima << setw(15) << glob_av[iH]/(double)iblk << setw(15) << err << setw(15) << mu << setw(15) << sigma << endl;

    cout << "----------------------------" << endl << endl;

    out.close();
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void Move(double sig, double m)
{
	x1=x0+rnd.Rannyu(-0.5,0.5)*passo;
	ratio=psi2(x1,sig,m)/psi2(x0,sig,m);
	if(rnd.Rannyu() <= min(1.,ratio)){
  	x0=x1;
    accepted++;
  }
  attempted++;
  return;
}

double psi(double x, double sig, double m)
{	 
	return exp( -((x-m)*(x-m))/(2.0*sig*sig) ) + exp( -((x+m)*(x+m))/(2.0*sig*sig) );
}

double psi2(double x, double sig, double m)
{
  return (psi(x,sig,m))*(psi(x,sig,m));
}

void Measure(double sig, double m)
{
  walker[iH]=(-0.5*psi_second(x0, sig, m)+potential(x0)*psi(x0, sig, m))/psi(x0, sig, m);
  norm++;
  int k = int((x0-xs)/double(binsize));
  histo[k]=histo[k]+1;
  return;
}

double psi_second(double x, double sig, double m)
{
	return (1/(sig*sig))*( -psi(x,sig,m)+(((x-m)*(x-m))/(sig*sig))*exp(-((x-m)*(x-m))/(2*sig*sig)) + (((x+m)*(x+m))/(sig*sig))*exp(-((x+m)*(x+m))/(2.0*sig*sig)) );
	/*double e1 = exp(-(x-m)*(x-m)/(double)(2*sig*sig) ) / (sig*sig);
	double e2 = exp(-(x+m)*(x+m)/(double)(2*sig*sig) ) / (sig*sig);
	double appo1 = (1.-(x-m)*(x-m)/(sigma*sig));
	double appo2 = (1.-(x+m)*(x+m)/(sigma*sig));
	return - e1*appo1 - e2*appo2;*/

}

double potential(double x)
{
  return pow(x,4)-5.0*x*x/2.0;
}

void PrintFinal()
{
	ofstream histogr;
  histogr.open("output.histo.dat");
  for(int i=0; i<nbin; ++i){
    histogr << xs+(i+0.5)*binsize << setw(15) << (double)histo[i]/((double)norm*binsize) << endl;
  }
  histogr.close();
  cout << endl;
	cout <<"*****************************"<<endl;
	cout<<"Estimated mu: " << mu << endl;
	cout<<"Estimated sigma: " << sigma << endl;
	cout <<"*****************************"<<endl;
	cout << endl;
  return;
}

