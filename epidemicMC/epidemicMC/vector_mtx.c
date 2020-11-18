double *vector_malloc(int nmax) {
	return (double *) malloc(sizeof(double)*nmax);
}

int *int_vector_malloc(int nmax) {
	return (int *)malloc(sizeof(int)*nmax);
}

double **mtx_malloc(int mmax, int nmax) {
	double **mtx;
	mtx = (double **) malloc(sizeof(double *)*mmax);
	for (int m = 0; m < mmax; m++) {
		mtx[m] = (double *)malloc(sizeof(double)*nmax);
	} // m-loop
	return mtx;
}

void mtx_free(double **mtx, int mmax) {
	for (int m = 0; m < mmax; m++) {
		free(mtx[m]);
	} // m-loop
	free(mtx);
}

void CopyVector(double *a, double *b, int nmax) {
	for (int n = 0; n < nmax; n++) {
		b[n] = a[n];
	}
	// we're working directly with memory, so no need to return anything
}