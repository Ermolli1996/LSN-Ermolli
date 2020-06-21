#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "sigma.h"

double error(double * AV, double * AV2, int n){
	if (n==0){
		return 0;
	}
	else return sqrt((double)(AV2[n] - pow(AV[n],2))/n);
}
