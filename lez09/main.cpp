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
	
	for(int i=0; i<ngenerations; i++){
		if(!Check()){
    	cerr << "Check error ..." << endl;
    	return -1;
  	}
		Selection();
		Crossover();
		Mutations(son1, son2);
		Update();
		if(i%100==0){
		 	cout << "Iteration number  " << i << " of " << ngenerations << endl;
			Print_evolution(i);
		}
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
	cout<<p1<<"  "<<p2<<endl;
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
  ReadInput >> popdim;
  ReadInput >> ngenerations;
  ReadInput >> type;
  ReadInput >> R;
  ReadInput >> pr1;		
  ReadInput >> pr2;
  ReadInput >> pr3;
  ReadInput >> pr4;
  ReadInput >> prcr;
  ReadInput >> p_elite;
  ReadInput.close();
  
  cout << "Traveling Salesman Problem             " << endl;
  cout << "Genetic algorithm             " << endl << endl;

  GeneratePosition(type);
  
  pop.resize(popdim);
  for(int i=0; i<popdim; i++){
    for(int j=0; j<ncities; j++){
      pop[i].push_back(int(j+1));
    }
    random_shuffle(pop[i].begin()+1, pop[i].end());
  }

  Sort();
}

//Selection: scelta dei percorsi padre e madre con probabilità proporzionale all'inverso della lunghezza
void Selection(){
	prob.resize(popdim);
	probcum.resize(popdim+1);
	probcum[0]=0;
	norm=0;
	cum=0;
	for(int i=0; i<popdim; i++){
		prob[i]=1/(pow(L2(pop[i]),3));
		norm+=prob[i];
	}
	for(int i=0; i<popdim; i++){
		prob[i]=prob[i]/norm;
		cum+=prob[i];
		probcum[i+1]=cum;
	}
	pmom=rnd.Rannyu(0,1);
	pdad=rnd.Rannyu(0,1);
	for(int k=0; k<popdim; k++){
		if(probcum[k]<pmom  && pmom<probcum[k+1]){
			mom=k;
			}
		if(probcum[k]<pdad && pdad<probcum[k+1]){
			dad=k;
			}
	}
}
	
		
		
//CROSSOVER: a partire da mom e dad creo son1 e son2 con una certa p
void Crossover(){
	if(rnd.Rannyu()<=prcr){
		n = int( rnd.Rannyu(1, ncities) );
		son1 = pop[mom];
  	son2 = pop[dad];
  	cont_mom=0;
  	cont_dad=0;

  	for(int j=0; j<ncities; j++){
    	pos_dad=true;
    	pos_mom=true;
    	for(int i=0; i<n; i++){
      	if(pop[dad][j]==pop[mom][i]) pos_dad=false;		//true se non trovo l'elemento di dad nella testa di mom
      	if(pop[mom][j]==pop[dad][i]) pos_mom=false;
    	}
    	if(pos_dad){
      	son1[n+cont_mom]=pop[dad][j];	//inserisco nella coda di mom gli elementi che gli mancano nell'ordine di dad
      	cont_mom++;
    	}
    	if(pos_mom){
      	son2[n+cont_dad]=pop[mom][j];
      	cont_dad++;
    	}
  	}
  }
  //Se non si verifica il crossover non posso fare le mutazioni ai figi! Restituisco dei figli che in realtà sono una l'estrazione di una traiettoria nella popolazione in modo gaussiano 
  else{
  	/*do{
      k1 = int(rnd.Gauss(pop.size()/3, pop.size()/4) );
    } while(k1 < 0 || k1 >= pop.size());
    do{
      k2 = int(rnd.Gauss(pop.size()/3, pop.size()/2));
    } while(k2 < 0 || k2 >= pop.size() || k2 == k1);
  	son1=pop[k1];
  	son2=pop[k2];*/
  	son1 = pop[mom];
  	son2 = pop[dad];
	}
}

//MUTATIONS
void Permutation_1(vector<int>& traj){				// traiettoria
  auto it1 = traj.begin();	//puntatore al primo elemento della traiettoria
  auto it2 = traj.begin();
  advance(it1,int(rnd.Rannyu(1,traj.size())));		//avanzo il puntatore di almeno 1
  advance(it2,int(rnd.Rannyu(1,traj.size())));
  iter_swap(it1, it2);		//cambio gli elementi puntati
}

void Shift(vector<int>& traj){
  n = int( rnd.Rannyu(1, traj.size()) );
  rotate(traj.begin()+1, traj.begin() + n, traj.end()); 		//shift tale che middle diventa il secondo elemento
}

void Inversion(vector<int>& traj){
  int n = int( rnd.Rannyu(1, traj.size()) );		//numero di elementi da scambiare
  int st = int( rnd.Rannyu(1, traj.size()-n) );	//posizione del primo da scambiare
  reverse(traj.begin()+st, traj.begin()+st+n);
}

void Permutation_m(vector<int>& traj){
  int n = int( rnd.Rannyu(1, traj.size()/2) );				//numero scambi
  int m = int( rnd.Rannyu(0, traj.size()-2*n) );			//distanza range (b-c)
  int st = int( rnd.Rannyu(1, traj.size()-2*n-m) );		//posizione primo scambiato
  swap_ranges(traj.begin()+st, traj.begin()+st+n, traj.begin()+st+n+m); //Cambio nel range [a,b) con quelli nel range che parte con c 
}

void elitary(){
	Sort();
	elite=pop[0];
	pop.push_back(elite);		//aggiungo elite
  Sort();
  pop.erase(pop.end());
}

//Mutations: applico la mutazione ai figli con la probabilità associata
void Mutations(vector<int>& traj1, vector<int>& traj2){
   r1=rnd.Rannyu();
   r2=rnd.Rannyu();
   if(r1<pr1) {Permutation_1(traj1); }
   if(r2<pr1) {Permutation_1(traj2); }
   r1=rnd.Rannyu();
   r2=rnd.Rannyu();
   if(r1<pr2) {Shift(traj1); }
   if(r2<pr2) {Shift(traj2); }
   r1=rnd.Rannyu();
   r2=rnd.Rannyu();
   if(r1<pr3) {Inversion(traj1); }
   if(r2<pr3) {Inversion(traj2); }
   r1=rnd.Rannyu();
   r2=rnd.Rannyu();
   if(r1<pr4) {Permutation_m(traj1); }
   if(r2<pr4) {Permutation_m(traj2); }
   r1=rnd.Rannyu();
   if(r1<p_elite) {elitary(); }
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

//Sort
void Sort(){
	struct {
  	bool operator() (vector<int> traj1, vector<int> traj2) {
  		return L2(traj1) < L2(traj2);
  	}
	}lenght;
  sort(pop.begin(), pop.end(), lenght);
}

//Update
void Update(){
	pop.push_back(son1);		//aggiungo i figli alla popolazione
  pop.push_back(son2);
  Sort();
  pop.erase(pop.end());		//elimino gli ultimi due elementi dopo aver ordinato
  pop.erase(pop.end()-1);
	
}

//Print evolution: stampa L2 mediato sulla miglior metà ad ogni iterazione
void Print_evolution(int num){
	t=0;
	for(int i=0; i<int(pop.size()/2); i++){
	 t+=L2(pop[i]);
	}
	L=t/(pop.size()/2);
	ofstream out;
	out.open("Evolution.res", ios::app);
	out << num << setw(15) << L << setw(15) << L2(pop.front()) <<endl; //Ngen, <L2>best_half, L2best
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

