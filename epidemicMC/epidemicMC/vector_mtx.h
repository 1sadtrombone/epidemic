#ifndef VECTOR_MTX_H
#define VECTOR_MTX_H

double *vector_malloc(int nmax);

double **mtx_malloc(int mmax, int nmax);

void mtx_free(double **mtx, int mmax);

void CopyVector(double *a, double *b, int nmax);

int *int_vector_malloc(int nmax);

#endif