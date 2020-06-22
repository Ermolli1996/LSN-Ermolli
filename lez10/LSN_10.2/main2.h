#include "./ParRandGen/random.h"
using namespace std;
#include <vector>

Random rnd;

vector<vector<int> > pop; //popolazione: vector di traiettorie (ciascuna vector di città identificate da un numero int) (popdim x ncities)
vector<vector<double> > coord; 	//Coordinate delle città in 2D (ncities x 2)
vector<double> prob; //associo ad ogni percorso nella popolazione una probabilità
vector<double> probcum; 	//la probabilità che il percorso k-esimo sia scelto è quella di estrarre un numero fra probcum[k] e probcum[k+1]
int ncities, popdim, ngenerations, type, n_migr; 		
double R, pr1, pr2, pr3, pr4, prcr, p_elite;
string f;
int mom, dad, n;
unsigned int k1, k2;
double cum, pmom, pdad, norm, r1, r2, t, L;
int z, y;
double local_best, global_best;

vector<int> son1;
vector<int> son2;
vector<int> elite;

int cont_mom, cont_dad;
bool pos_dad, pos_mom;


void input();
void settingRnd(int );
void Selection();
void Crossover();
void Permutation(vector<int>& );
void Shift(vector<int>& );
void Inversion(vector<int>& );
void Permutation_m(vector<int>& );
void elitary();
void Mutations(vector<int>&, vector<int>& );
void GeneratePosition(int );
double L2(vector<int> );
bool Check();
bool myfunction (vector<int> , vector<int> );
void Sort();
void Update();
void Print_evolution(int, int );
void Print_best(vector<int>, int );
int Pbc(int );
