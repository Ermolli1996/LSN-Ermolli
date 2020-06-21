#include "../lez00/ParRandGen/random.h"
Random rnd;

int M;
int N;
double passo;
double x_start, x0, x1;
double b_start, b0, b1;
double z_start, z0, z1;
double ratio;
int state;
char distribution;

void input();
void settingRnd();
double error(double *, double *, int);
double psi_1s(double , double , double);
double psi_2p(double , double , double);
