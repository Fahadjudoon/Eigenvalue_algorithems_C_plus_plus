#include<iostream>
#include"vectors.h"
#include<assert.h>
#include<time.h>
#include <random>

	vectors::vectors(unsigned int size_x) : data(size_x, 0.), size_x_val(size_x) {}

	double& vectors::operator ()(unsigned int x) { assert(x < size_x_val && "you try to access vector-element that does not exist"); return data[x]; }
	const double& vectors::operator ()(unsigned int x)const { assert(x < size_x_val && "you try to access matrix-element that does not exist"); return data[x]; }


	unsigned int vectors::size(void)const { return size_x_val; }
	void vectors::print(void)const { 
		std::cout << std::endl;
		for (unsigned int i = 0; i < size(); i++) { std::cout << data[i] << std::endl; }}

	vectors random_vec(unsigned int r) {
		std::default_random_engine generator;
		std::normal_distribution<double> distribution(0.0, 1.0);
		vectors random(r);
		for (unsigned int i = 0; i < random.size(); i++)
		{	
				random(i) = distribution(generator);
		}
		return random;
	}

	vectors vec_one(void) {
		vectors vec(3);
		vec(0) = 1.;
		vec(1) = 1.;
		vec(2) = 1.;
		return vec;
	}
	vectors vec_two(void) {
		vectors vec(2);
		vec(0) = 3.;
		vec(1) = 4.;
		return vec;
	}

	vectors operator*(const vectors &lhs, double rhs) {
		vectors output = lhs;
		for (unsigned int i = 0; i < lhs.size(); i++) { output(i) *= rhs; }
		return output;
	}
	vectors operator/(const vectors &lhs, double rhs) {
		vectors output = lhs;
		for (unsigned int i = 0; i < lhs.size(); i++) { output(i) /= rhs; }
		return output;
	}
	void operator/=(vectors &lhs, double rhs) {
		for (unsigned int i = 0; i < lhs.size(); i++) { lhs (i) = lhs(i) / rhs; }
	}
	vectors operator*(double lhs, const vectors &rhs) {
		vectors output = rhs;
		for (unsigned int i = 0; i < rhs.size(); i++) { output(i) *= lhs; }
		return output;
	}
	
	void operator*=(vectors &lhs, double rhs) {
		for (unsigned int i = 0; i < lhs.size(); i++) { lhs(i) *= rhs; }
	}

	vectors operator*(const vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		vectors output(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++) { output(i) = lhs(i) * rhs(i); }
		return output;
	}

	void operator*=(vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++) { lhs(i) *= rhs(i); }
	}

	vectors operator+(const vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		vectors output(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++) { output(i) = lhs(i) + rhs(i); }
		return output;
	}

	void operator+=(vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++) { lhs(i) += rhs(i); }
	}

	vectors operator-(const vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		vectors output(lhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++) { output(i) = lhs(i) - rhs(i); }
		return output;
	}

	void operator-=(vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		for (unsigned int i = 0; i < lhs.size(); i++) { lhs(i) -= rhs(i); }
	}

	/*bool operator==(vectors &lhs, const vectors &rhs) {
		assert(lhs.size() == rhs.size());
		bool result;
		
		if (lhs == rhs)
			return (1);
		else
			return(0);
	}*/

	

	