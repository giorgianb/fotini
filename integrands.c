#include <stdio.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265358979323846

struct data {
	int n_terms;
	double *norms;
	double *coeffs_real;
	double *coeffs_imag;
	int n_permutation;
	int *permutations; 
	double *t_0;
	double *ω_0;
	double *Δ;
	int printed;

};

static inline double complex ξ_gaussian(int n, double *t, double *t_0, double *ω_0, double *Δ, int *permutation) {
	double complex sum = 0;
	for (int i = 0; i < n; ++i) {
		double δ = t_0[i] - t[permutation[i]];
		sum += -I*ω_0[i]*t[permutation[i]] - Δ[i]*Δ[i]*δ*δ;
	}

	sum = cexp(sum);

	return sum;
}

double β_gaussian(int n, double *xx, void *user_data) {
	struct data *data = user_data;

	double complex total = 0;
	double ξ_coeff = 1;
	for (int i = 0; i < n; ++i)
		ξ_coeff *= (2*data->Δ[i]*data->Δ[i]/PI);

	ξ_coeff = sqrt(sqrt(ξ_coeff));


	for (int i = 0; i < data->n_terms; ++i) {
		double complex subterm = 0;
		for (int j = 0; j < data->n_permutation; ++j)
			subterm += ξ_gaussian(n, xx, data->t_0, data->ω_0, data->Δ, &data->permutations[i*data->n_permutation*n + n*j]);

		double complex coeff = ξ_coeff*(data->coeffs_real[i] + data->coeffs_imag[i]*I);
		total += subterm*data->norms[i]*coeff/data->n_permutation;
	}

	data->printed = 1;

	return creal(total * conj(total));
}

// For testing
double ξ_r(int n, double *t, double *t_0, double *ω_0, double *Δ, int *permutation) {
	double complex sum = 0;
	for (int i = 0; i < n; ++i) {
		double δ = t_0[i] - t[permutation[i]];
		sum += -I*ω_0[i]*t[permutation[i]] - Δ[i]*Δ[i]*δ*δ;
	}

	double complex product = 1;
	for (int i = 0; i < n; ++i)
		product *= (2*Δ[i]*Δ[i]/PI);

	product = sqrt(sqrt(product));
	sum = cexp(sum);

	return creal(product * sum);
}

double ξ_i(int n, double *t, double *t_0, double *ω_0, double *Δ, int *permutation) {
	double complex sum = 0;
	for (int i = 0; i < n; ++i) {
		double δ = t_0[i] - t[permutation[i]];
		sum += -I*ω_0[i]*t[permutation[i]] - Δ[i]*Δ[i]*δ*δ;
	}

	double complex product = 1;
	for (int i = 0; i < n; ++i)
		product *= (2*Δ[i]*Δ[i]/PI);

	product = sqrt(sqrt(product));
	sum = cexp(sum);

	return cimag(product * sum);
}
