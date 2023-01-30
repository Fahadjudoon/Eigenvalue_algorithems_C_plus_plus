#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

#include "vectors.h"
#include "matrix.h"

vectors mat_vec_prod(const matrix &, const vectors &);
matrix vec_vec_trans_prod(const vectors &v, const matrix &m);
double vec_trans_vec_prod(const matrix &lhs, const vectors &rhs);
double max_val(const vectors &v);
double signum(double x);
double vec_norm(const vectors &v);
vectors compute_eig_vec(const matrix &B, vectors &x_start);
matrix vec_trans(const vectors &v);
double vec_dot_prod(const vectors &v1, const vectors &v2);
//void eigen_v(const matrix &A, vectors &x);
void compute_all_eigenvectors(matrix &B);
//matrix calculate_Q(const matrix &A);
matrix calc_Q(const matrix &A);
matrix calculate_R(const matrix &A, const matrix &Q);
void eig(matrix &A);
matrix qr_givens_rotation(const matrix &H);			//This functon apply alg 4.2 on Hessenberg for qr factorization
#endif // LINEAR_ALGEBRA_H