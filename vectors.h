#ifndef VECTORS_H
#define VECTORS_H

#include<vector>

class vectors {
private:
	std::vector<double> data;
	unsigned int size_x_val;
public:
	vectors(unsigned int);

	double& operator ()(unsigned int);
	const double& operator ()(unsigned int)const;


	unsigned int size(void)const;
	void print(void)const;
};
vectors vec_one(void);
vectors vec_two(void);
vectors random_vec(unsigned int r);

vectors operator*(const vectors &lhs, double rhs);
void operator*=(vectors &lhs, double rhs);

vectors operator/(const vectors &lhs, double rhs);
void operator/=(vectors &lhs, double rhs);

vectors operator*(double lhs, const vectors &rhs);

vectors operator*(const vectors &lhs, const vectors &rhs);
void operator*=(vectors &lhs, const vectors &rhs);


vectors operator+(const vectors &lhs, const vectors &rhs);
void operator+=(vectors &lhs, const vectors &rhs);

vectors operator-(const vectors &lhs, const vectors &rhs);
void operator-=(vectors &lhs, const vectors &rhs);

//bool operator==(vectors &lhs, const vectors &rhs);


#endif // VECTORS_H