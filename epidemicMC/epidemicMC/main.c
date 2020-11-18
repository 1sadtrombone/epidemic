#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "vector_mtx.h"
#include "write.h"

int main(int argc, char **argv)
{
	srand(time(NULL));

	int N = 200;
	int Lx = 1000;
	int Ly = 1000;
	int tmax = 5000;

	double p = 0.55;

	char *outfile = "epidemic_above_pc.dat";

	double **Z = mtx_malloc(tmax, Lx);

	Z = MC(p, N, tmax, Lx, Ly);

	WriteResults(outfile, Z, tmax, Lx);

	return 0;
}