#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "../lez00/ParRandGen/random.h"
#include "main.h"
#include <iterator>		//advance
#include <algorithm>	//rotate, iter_swap, swap_ranges, reverse, sort

using namespace std;
 
int main (){
	
  settingRnd();
	input();
	temp=T_in;
	for(int i=0; i<niterations; i++){
		if(!Check()){
    	cerr << "Check error ..." << endl;
    	return -1;
  	}
		Mutations();
		Print_evolution(i);
		temp=temp/fract;
	}
	
	if(!Check()){
    cerr << "Check error ..." << endl;
    return -1;
  }
 	
	Print_best(pop[0]);
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

//Input: Leggo input, Genero le posizioni delle città, Genero la popolazione iniziale casualmente e la ordino
void input(){
  ifstream ReadInput;
  if(!ReadInput.good()){
     cerr<<"Error reading file"<<endl;
  }
  ReadInput.open("input.dat");
  ReadInput >> ncities;
  ReadInput >> niterations;
  ReadInput >> type; 
  ReadInput >> R; 
  ReadInput >> T_in; 
  ReadInput >> fract;
  ReadInput >> steps_forT;
  ReadInput.close();
  
  cout << "Traveling Salesman Problem             " << endl <<endl;
  cout << "Simulated annealing             " << endl << endl;
	
	popdim=1;
	
  GeneratePosition(type);
  
  pop.resize(popdim);
  for(int i=0; i<popdim; i++){
    for(int j=0; j<ncities; j++){
      pop[i].push_back(int(j+1));
    }
    random_shuffle(pop[i].begin()+1, pop[i].end());
  }
}

//MUTATIONS
void Permutation_1(vector<int>& traj){				// traiettoria
  old_traj=traj;
  auto it1 = traj.begin();	//puntatore al primo elemento della traiettoria
  auto it2 = traj.begin();
  advance(it1,int(rnd.Rannyu(1,traj.size())));		//avanzo il puntatore di almeno 1
  advance(it2,int(rnd.Rannyu(1,traj.size())));
  iter_swap(it1, it2);		//cambio gli elementi puntati
  p=min(1., exp(-( L2(traj) - L2(old_traj) )/temp));
  r=rnd.Rannyu();
  if(r<=p) succ1++;		//successo! Lascio la nuova traj
  else traj=old_traj;	//Mutazione fallita! Riscelgo la old_traj
  att1++;
}

void Shift(vector<int>& traj){
	old_traj=traj;
  double n = int( rnd.Rannyu(1, traj.size()) );
  rotate(traj.begin()+1, traj.begin() + n, traj.end()); 		//shift tale che middle diventa il secondo elemento
  p=min(1., exp(-( L2(traj) - L2(old_traj) )/temp));
  r=rnd.Rannyu();
  if(r<=p) succ2++;		//successo! Lascio la nuova traj
  else traj=old_traj;	//Mutazione fallita! Riscelgo la old_traj
  att2++;
}

void Inversion(vector<int>& traj){
  old_traj=traj;
  int n = int( rnd.Rannyu(1, traj.size()) );		//numero di elementi da scambiare
  int st = int( rnd.Rannyu(1, traj.size()-n) );	//posizione del primo da scambiare
  reverse(traj.begin()+st, traj.begin()+st+n);
  p=min(1., exp(-( L2(traj) - L2(old_traj) )/temp));
  r=rnd.Rannyu();
  if(r<=p) succ3++;		//successo! Lascio la nuova traj
  else traj=old_traj;	//Mutazione fallita! Riscelgo la old_traj
  att3++;
}

void Permutation_m(vector<int>& traj){
  old_traj=traj;
  int n = int( rnd.Rannyu(1, traj.size()/2) );				//numero scambi
  int m = int( rnd.Rannyu(0, traj.size()-2*n) );			//distanza range (b-c)
  int st = int( rnd.Rannyu(1, traj.size()-2*n-m) );		//posizione primo scambiato
  swap_ranges(traj.begin()+st, traj.begin()+st+n, traj.begin()+st+n+m); //Cambio nel range [a,b) con quelli nel range che parte con c 
	p=min(1., exp(-( L2(traj) - L2(old_traj) )/temp));
  r=rnd.Rannyu();
  if(r<=p) succ4++;		//successo! Lascio la nuova traj
  else traj=old_traj;	//Mutazione fallita! Riscelgo la old_traj
  att4++;
}

//Mutations: applico la mutazione ai figli con la probabilità associata
void Mutations(){
	att1=0, att2=0, att3=0, att4=0;
	succ1=0, succ2=0, succ3=0, succ4=0;
	for(int i=0; i<steps_forT; i++){
   	Permutation_1(pop[0]);
   	Shift(pop[0]);
   	Inversion(pop[0]);
   	Permutation_m(pop[0]);
	}
	cout << "T = " << temp << endl;
  cout << "Permutation_1 acceptance ratio : " << (double)succ1/att1 << endl;
  cout << "Shift acceptance ratio : " << (double)succ2/att2 << endl;
  cout << "Inversion acceptance ratio : " << (double)succ3/att3 << endl;
  cout << "Permutation_m acceptance ratio : " << (double)succ4/att4 << endl << endl;
}

//Coordinates of the cities
void GeneratePosition(int type){

	coord.resize(ncities);
  if (type==1){
    for(int i=0; i<ncities; i++){
    	double theta=rnd.Rannyu(0.,2*M_PI);
      coord[i].push_back(R*cos(theta));
      coord[i].push_back(R*sin(theta));
    }
  }

  else if (type==2) {
    for(int i=0; i<ncities; i++){
      coord[i].push_back(rnd.Rannyu(-R,R));
      coord[i].push_back(rnd.Rannyu(-R,R));
    }
  }
}

//Cost function
double L2(vector<int> traj){
  double sum=0;
  for(unsigned int i=0; i<traj.size(); i++){
    double x1=coord[traj[i]-1][0];
    double y1=coord[traj[i]-1][1];
    double x2=coord[traj[Pbc(i+1)]-1][0];
    double y2=coord[traj[Pbc(i+1)]-1][1];
    sum+=fabs((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  }
  return sum;
}

//Check
bool Check(){

  for(unsigned int i=0; i<pop.size(); ++i){
  	if(pop[i][0] != 1) return false;
    for(unsigned int j=0; j<pop[i].size()-1; ++j){
      for(unsigned int z=j+1; z<pop[i].size(); ++z){
        if(pop[i][j]==pop[i][z]) return false;
      }
    }
  }

  return true;
}

//Print evolution: stampa L2 mediato sulla miglior metà ad ogni iterazione
void Print_evolution(int num){
	ofstream out;
	out.open("Evolution.res", ios::app);
	out << num << setw(15) << L2(pop.front()) <<endl; //Niter, L2best
	out.close();
}
	
//Print best: stampa le coordinate del percorso migliore
void Print_best(vector<int> traj_best){
	ofstream outbest;
	outbest.open("Best_coordinates.res");
	for(unsigned int i=0; i<traj_best.size(); i++){
		outbest << coord[traj_best[i]-1][0] << setw(15) <<  coord[traj_best[i]-1][1] << endl;
	}
	outbest.close();
	cout<<"Finish! Take a look at output files"<<endl;
}	

int Pbc(int i){
    if(i >= ncities) i = i - ncities;
    else if(i < 0) i = i + ncities;
    return i;
}

