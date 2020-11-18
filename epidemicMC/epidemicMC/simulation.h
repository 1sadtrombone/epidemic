#pragma once

double **MC(double p, int N, int tmax, int Lx, int Ly);
double drand(void);
double zero_to_one(void);
int is_already_infected(int x, int y, int *Yx, int* Yy, int nmax);