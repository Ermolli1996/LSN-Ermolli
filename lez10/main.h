#include "../lez00/ParRandGen/random.h"
using namespace std;
#include <vector>

Random rnd;

vector<vector<int>> pop; //popolazione: vector di traiettorie (ciascuna vector di città identificate da un numero int) (popdim x ncities)
vector<vector<double>> coord; 	//Coordinate delle città in 2D (ncities x 2)
int ncities, popdim, niterations, type, steps_forT; 		
double R, T_in, p, r, temp, fract;

double succ1, succ2, succ3, succ4;
double att1, att2, att3, att4;

vector<int> old_traj;

void input();
void settingRnd();
void Permutation(vector<int>& );
void Shift(vector<int>& );
void Inversion(vector<int>& );
void Permutation_m(vector<int>& );
void Mutations();
void GeneratePosition(int );
double L2(vector<int> );
bool Check();
void Print_evolution(int );
void Print_best(vector<int> );
int Pbc(int );
