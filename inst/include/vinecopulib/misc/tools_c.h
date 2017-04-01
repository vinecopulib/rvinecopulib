#ifndef TOOLS_C_H
#define TOOLS_C_H

#ifdef __cplusplus
extern "C" {
#endif

void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V);
void ktau_matrix(double *data, int *d, int *N, double *out);
double debyen(const double x, const int n);

#ifdef __cplusplus
}
#endif

#endif
