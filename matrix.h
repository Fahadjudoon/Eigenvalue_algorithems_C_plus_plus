#ifndef MATRIX_H
#define MATRIX_H

#include<vector>

class matrix{
public:
	std::vector<double> data;
	unsigned int size_x_val, size_y_val;
public:
	matrix(unsigned int, unsigned int);

	double& operator()(unsigned int, unsigned int);
	const double& operator()(unsigned int, unsigned int)const;

	unsigned int size_x(void)const;
	unsigned int size_y(void)const;

	void print(void)const;
	
};

matrix mat_five(void);
matrix mat_ten(void);
matrix mat_fifteen(void);
matrix mat_twenty(void);
matrix mat_thirty(void);
matrix mat_fifty(void);
matrix mat_seventy_five(void);
matrix mat_hundred(void);
matrix mat_hundred_fifty(void);
matrix mat_two_hundred(void);


matrix operator*(const matrix &lhs, double rhs);
void operator*=(matrix &lhs, double rhs);

matrix operator/(const matrix &lhs, double rhs);
void operator/=(matrix &lhs, double rhs);

matrix operator*(double lhs, const matrix &rhs);

matrix operator+(const matrix &lhs, const matrix &rhs);
void operator+=(matrix &lhs, const matrix &rhs);

matrix operator-(const matrix &lhs, const matrix &rhs);
void operator-=(matrix &lhs, const matrix &rhs);

void operator==(matrix &lhs, const matrix &rhs);

matrix mat_mat_prod(const matrix & lhs, const matrix & rhs);

matrix transpose(const matrix &m);

matrix random_mat(unsigned int r, unsigned int c);
matrix identity_mat(unsigned int r, unsigned int c);
#endif // MATRIX_H