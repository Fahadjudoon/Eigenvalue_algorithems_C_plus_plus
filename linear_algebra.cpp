#include "linear_algebra.h"
#include <assert.h>
#include <iostream>
#include <vector>


vectors mat_vec_prod(const matrix & lhs, const vectors & rhs) { 
	assert (lhs.size_y() == rhs.size());
		vectors output(rhs.size());
		for (unsigned int i = 0; i < lhs.size_x(); i++)
		{
			for (unsigned int j = 0; j < rhs.size(); j++)
			{
				output(i) += (lhs(i, j) * rhs(j));
			}
		}
		return output;
	}

//To take max value as common
double max_val(const vectors &v)
{
	double max_val = 0;
	for (unsigned int i = 0; i < v.size(); i++)
	{
		if (v(i) > max_val)
		{
			max_val = v(i);
		}
	}
	return max_val;
}
//Sign returning function
double signum(double x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

double vec_norm(const vectors &v) 
{
	double temp = 0.;
	for (unsigned int i = 0; i < v.size(); i++)
	{
		temp += v(i) * v(i);
	}
	double length = sqrt(temp);
	return length;
}
vectors compute_eig_vec(const matrix &B, vectors & x_start)
{
	for ( unsigned int i = 0; i < B.size_y(); i++)
	{
		x_start = mat_vec_prod(B, x_start);
		x_start /= vec_norm(x_start);
	}
	return x_start;
}

matrix vec_trans(const vectors &v) {
	matrix v_trans(1, v.size());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		v_trans(0, i) = v(i);
	}
	return v_trans;
}
/*
matrix vec_vec_trans_prod(const vectors &v, const matrix &m)
{
	assert(m.size_y() == v.size());
	matrix output(v.size(), v.size());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		for (unsigned int j = 0; j < m.size_y(); j++)
		{
			output(i,j) = (v(j) * m(0,j));
			std::cout << i << j << output(i, j) << std::endl;
		}
	}
	return output;
}
*/
matrix vec_vec_trans_prod(const vectors &v, const matrix &m)
{
	assert(m.size_y() == v.size());
	matrix output(v.size(), v.size());
	for (unsigned int i = 0; i < v.size(); i++)
	{
		for (unsigned int j = 0; j < m.size_y(); j++)
		{
			output(i, j) = (v(i) * m(0, j));
		}
	}
	return output;
}
double vec_trans_vec_prod(const matrix &lhs, const vectors &rhs)
{
	double output = 0.;
	for (unsigned int i = 0; i < rhs.size(); i++) { output += lhs(0,i) * rhs(i); }
	return output;
}

/*void eigen_v(const matrix &A, vectors &x) {
	do {
		vectors v = mat_vec_prod(A, x);
		std::cout << "Matrix vector Product" << std::endl;
		v.print();

		double z = max_val(v);
		std::cout << std::endl;
		std::cout << "max common value is: " << z << std::endl;

		v /= vec_norm(v);
		v.print();
		std::cout << std::endl;

		for (unsigned int j = 0; j < 5; j++) {
			
			vectors v1 = mat_vec_prod(A, v);
			std::cout << "Matrix vector Product" << std::endl;
			v1.print();
			
			double z = max_val(v1);
			std::cout << std::endl;
			std::cout << "max common value is: " << z << std::endl;

				v1 /= vec_norm(v1);
				v1.print();
				std::cout << std::endl;
			}
	} while (0);
}*/

void compute_all_eigenvectors(matrix &A) {
	matrix B = mat_mat_prod(transpose(A), A);
	B.print();
	double eigen = 0.0;
	std::vector <vectors> xi;
	xi.reserve(B.size_y());
	for (unsigned int i = 0; i < B.size_y(); i++)
	{
		xi.push_back(vectors(B.size_x()));
	}
	for (unsigned int i = 0; i < B.size_x(); i++)
	{
		vectors x_start = random_vec(B.size_x());
		xi[i] = compute_eig_vec(B, x_start);
		//std::cout << "first xi" << std::endl;
		//xi[i].print();
		eigen = vec_norm(mat_vec_prod(B, xi[i])) / vec_norm(xi[i]);
		std::cout << std::endl;
		std::cout << eigen << std::endl;
		xi[i] = xi[i] / vec_norm(xi[i]);

		xi[i].print();
		B = B - (eigen * vec_vec_trans_prod(xi[i], vec_trans(xi[i])));
	}
}
/*matrix calculate_Q (const matrix &A)
{
	vectors v1(A.size_x());
	vectors v2(A.size_x());
	vectors v3(A.size_x());
	matrix Q(A.size_x(), A.size_y());
	//Determining u1
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		v1(i) = A(i, 0);
	}

	vectors u1 = v1 / vec_norm(v1);
		//std::cout << "u1 is: " << std::endl;
		//u1.print();
		//Determining u2
		for (unsigned int i = 0; i < v1.size(); i++)
		{
			v2(i) = A(i, 1);
		}
		vectors u2 = (v2 - (u1 * v2) * u1) / vec_norm((v2 - (u1 * v2) * u1));
		//std::cout << "u2 is: " << std::endl;
		//u2.print();
		//Determining u3
		for (unsigned int i = 0; i < v1.size(); i++)
		{
			v3(i) = A(i, 2);
		}
		
		vectors upper_term_u3 = (v3 - (((u2 * v3) * u2) - ((u1 * v3) * u1)));
		double u3_length = vec_norm(upper_term_u3);
		vectors u3 = (v3 - (((u2 * v3) * u2) - ((u1 * v3) * u1))) / vec_norm((v3 - (((u2 * v3) * u2) - ((u1 * v3) * u1))));
		//std::cout << "u3 is: " << std::endl;
		//u3.print();
		
		for (unsigned int i = 0; i < u1.size(); i++)
		{
			Q(i, 0) = u1(i);
		}
		for (unsigned int i = 0; i < v2.size(); i++)
		{
			Q(i, 1) = u2(i);
		}
		for (unsigned int i = 0; i < v3.size(); i++)
		{
			Q(i, 2) = u3(i);
		}	
	return Q;
}
*/
matrix calc_Q(const matrix &A)
{
	matrix Q(A.size_x(), A.size_y());
	std::vector <vectors> v;
	v.reserve(A.size_y());
	for (unsigned int i = 0; i < A.size_y(); i++)
	{
		v.push_back(vectors(A.size_x()));
		for (unsigned int j = 0; j < A.size_x(); j++)
		{
			v[i](j) = A(j,i);
		}
		//v[i].print();
	}
	std::vector <vectors> u;
	u.reserve(A.size_y());
	for (unsigned int i = 0; i < A.size_y(); i++)
	{   u.push_back(vectors(A.size_x()));  }
	u[0] = v[0] / vec_norm(v[0]);
	//std::cout << std::endl;
	//std::cout <<"u[0] is" << std::endl;
	//u[0].print();

	std::vector <vectors> uu;
	uu.reserve((A.size_y()-1));
	for (unsigned int i = 0; i < A.size_y()-1; i++)
	{   uu.push_back(vectors(A.size_x()));}

	std::vector <vectors> vec_proj;
	vec_proj.reserve((A.size_y() - 1));
	for (unsigned int i = 0; i < A.size_y() -1; i++)
	{	vec_proj.push_back(vectors(A.size_x())); }
	for (unsigned int i = 0; i < A.size_y() -1  ; i++)
	{
		unsigned int k = A.size_y()-1;		
			if (k == A.size_y()-1)
			{
				vec_proj[i] = (vec_dot_prod(u[i], v[i + 1])) * u[i];
				//std::cout << std::endl;
				//std::cout << "vec_proj inner" << i << "is" << std::endl;
				//vec_proj[i].print();

				uu[i] = v[i + 1] - vec_proj[i];
				//std::cout << std::endl;
				//std::cout << "uu inner is" << std::endl;
				//uu[i].print();
			}
			for (unsigned int j = i ; j > 0; j--)
			{
				vec_proj[i] = (vec_dot_prod(u[j-1], v[i + 1])) * u[j-1];
				//std::cout <<"sasdasdasd" << std::endl;
				//std::cout << "vec_proj middel " << i << "is" << std::endl;
				//vec_proj[i].print();

				uu[i] -= vec_proj[i];
				//std::cout << "uu middel is" << std::endl;
				//uu[i].print();
			}
			/*uu[i+1] = v[i + 1] - vec_proj[i];
			std::cout << std::endl;
			std::cout << "uu is" << std::endl;
			uu[i+1].print(); */

			u[i+1] = uu[i] / vec_norm(uu[i]); 
			//std::cout << "u outer is" << std::endl;
		   // u[i+1].print();

		/*std::cout << std::endl;
		std::cout << "uu is" << std::endl;
		uu[i + 1].print();
		u[i + 1] = uu[i + 1] / vec_norm(uu[i + 1]);
		std::cout << "u is" << std::endl;
		u[i + 1].print();*/
	}
	for (unsigned int i = 0; i < A.size_y(); i++)
	{
		for (unsigned int j = 0; j < A.size_x(); j++)
		{
			Q(j, i) = u[i](j);
		}
	}
	//Q.print();
	return Q;
}
matrix calculate_R (const matrix &A, const matrix &Q)
{
	matrix Q_trans = transpose(Q);
	//std::cout << "Q_trans is" << std::endl;
	//Q_trans.print();
	matrix R = mat_mat_prod(Q_trans, A);
	return R;
}

double vec_dot_prod(const vectors &v1, const vectors &v2) {
	double temp = 0.;
	for (unsigned int i = 0; i < v1.size(); i++)
	{
		temp += v1(i) * v2(i);
	}
	//std::cout <<"dot product is " << temp << std::endl;
	return temp;
}

void eig(matrix &A)
{
	matrix Q = matrix(A.size_x(), A.size_y());
	matrix R = matrix(A.size_x(), A.size_y());
	matrix T = matrix(A.size_x(), A.size_y());
	matrix U = matrix(A.size_x(), A.size_y());
	unsigned int n;

	std::cout << "Type in number of iterations: " << std::endl;          
	std::cin >> n;

	std::vector <matrix> A_;
	A_.reserve(n);
	for (unsigned int i = 0; i < n; i++)
	{
		A_.push_back(matrix(A.size_x(), A.size_y()));
	}

	std::vector <matrix> Q_;
	Q_.reserve(n);
	for (unsigned int i = 0; i < n; i++)
	{
		Q_.push_back(matrix(A.size_x(), A.size_y()));
	}

	std::vector <matrix> R_;
	R_.reserve(n);
	for (unsigned int i = 0; i < n; i++)
	{
		R_.push_back(matrix(A.size_x(), A.size_y()));
	}

	std::vector <matrix> U_;
	U_.reserve(n);
	for (unsigned int i = 0; i <= n; i++)
	{
		U_.push_back(matrix(A.size_x(), A.size_y()));
	}

	A_[0] = A; U_[0] = identity_mat(A.size_x(), A.size_y());

	unsigned int k;
	for (k = 1; k < n; k++)
	{
		Q_[k] = calc_Q(A_[k - 1]);
		R_[k] = calculate_R(A_[k - 1], Q_[k]);
		//R.print();
		//A_[k - 1] = mat_mat_prod(Q_[k], R_[k]);
		//std::cout << "new Matrix A_[k-1] with" << "iteration number " << k - 1 << " is " << std::endl;
		//std::cout << std::endl;
		//A_[k - 1].print();

		A_[k] = mat_mat_prod(R_[k], Q_[k]);
		//std::cout << "new Matrix A_[k] with" << "iteration number " << k << " is " << std::endl;
		//A_[k].print();
		

		U_[k] = mat_mat_prod(U_[k - 1], Q_[k]);
		//std::cout << "Matrix Q_[k] with" << "iteration number " << k << " is " << std::endl;
		//U_[k].print();
		//Q_[k].print();

		T = A_[k];
		U = U_[k];
		//T.print();
		//U.print();
		if (k == n - 1)
		{
			for (unsigned int j = 0; j < 1; j++)
			{
				//std::cout << std::endl;
				matrix G = mat_mat_prod(U, T);
				matrix H = transpose(U);
				A = mat_mat_prod(G, H);
				//std::cout << "The last Matrix A (R * Q)is " << std::endl;
				//A.print();
			}
		}
	}
	std::cout << "Matrix R_Q after Basic QR is" << std::endl;
	T.print();		//extra just to check, ater should be deleted
}



matrix qr_givens_rotation(const matrix &H) {
	matrix R = matrix(H.size_x(), H.size_y());
	R = H;

	matrix G = identity_mat(H.size_x(), H.size_y());
	//matrix K = matrix(H.size_x(), H.size_y());



	std::vector <matrix> Q_;
	Q_.reserve(H.size_x());
	for (unsigned int i = 0; i < H.size_x(); i++)
	{
		Q_.push_back(matrix(H.size_x(), H.size_y()));
	}
	Q_[0] = identity_mat(2, 2);

	std::vector <matrix> G_;
	G_.reserve(H.size_x() - 1);
	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		G_.push_back(matrix(2, 2));
	}

	std::vector <matrix> H_rows;
	H_rows.reserve(H.size_x() - 1);
	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		H_rows.push_back(matrix(2, H.size_y()));
	}

	std::vector <matrix> G_T;
	G_T.reserve(H.size_x() - 1);
	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		G_T.push_back(matrix(2, 2));
	}

	double c = 0.;
	double s = 0.;

	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		//Wird zu 2*2 Matrix

			c = R(i, i) / (sqrt((R(i, i) * R(i, i)) + (R(i + 1, i) * R(i + 1, i))));
			s = -1 * (R(i + 1, i)) / (sqrt((R(i, i) * R(i, i)) + (R(i + 1, i) * R(i + 1, i))));

		//std::cout << std::endl;
		//std::cout << c << "   " << s << std::endl;


		G_[i](0, 0) = c;		G_[i](0, 1) = (-1 * s);	  G_[i](1, 0) = s;		G_[i](1, 1) = c;
		//std::cout << "G_ matrix in above loop is: " << std::endl;
		//G_[i].print();

		//std::cout << std::endl;


		for (unsigned int k = i; k < 2 + i; k++)
		{
			for (unsigned int l = 0; l < H.size_y(); l++)
			{
				H_rows[i](k - i, l) = R(k, l);
			}
		}

		//std::cout << "Matrix H_rows is" << std::endl;
		//H_rows[i].print();
		//Anwendung von 2*2 Matrix nur auf die beiden (2) entsprechenden Zeilen von H

		H_rows[i] = mat_mat_prod(G_[i], H_rows[i]);
		//std::cout << "Matrix H_rows after multiplication with G_ is" << std::endl;
		//H_rows[i].print();

		for (unsigned int k = i; k < 2 + i; k++)
		{
			for (unsigned int l = 0; l < H.size_y(); l++)
			{
				R(k, l) = H_rows[i](k - i, l);
			}
		}
		//std::cout << "Matrix R is" << std::endl;
		//R.print();

	}
	//std::cout << "R_Q matrix calculated inside qr_givens_rotation func is: " << std::endl;
	//R.print();
	//std::cout << std::endl;

	matrix Hh = matrix(H.size_x(), H.size_y());
	Hh = R;

	std::vector <matrix> Hh_sub;
	Hh_sub.reserve(H.size_x() - 1);
	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		Hh_sub.push_back(matrix(i + 2, 2));
	}



	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		//std::cout << "G_ matrix in under loop is: " << std::endl;
		//G_[i].print();

		G_T[i] = transpose(G_[i]);
		//std::cout << "G_ T in under loop is: " << std::endl;
		//G_T[i].print();

		for (unsigned int k = 0; k < i + 2; k++)
		{
			for (unsigned int l = i; l < i + 2; l++)
			{
				Hh_sub[i](k, l - i) = Hh(k, l);
			}
		}
		//std::cout << "Hh_sub matrix is: " << std::endl;
		//Hh_sub[i].print();

		Hh_sub[i] = mat_mat_prod(Hh_sub[i], G_T[i]);
		//std::cout << "Hh_sub matrix after multiplication is: " << std::endl;
		//Hh_sub[i].print();

		for (unsigned int k = 0; k < i + 2; k++)
		{
			for (unsigned int l = i; l < i + 2; l++)
			{
				Hh(k, l) = Hh_sub[i](k, l - i);
			}
		}
		//std::cout << "Matrix at the end: (inside second loop) of Givens Rotation is" << std::endl;
		//Hh.print();
	}

	return (Hh);
}

/*
matrix qr_givens_rotation(const matrix &H) {
	matrix R = matrix(H.size_x(), H.size_y());
	R = H;

	matrix G = identity_mat(H.size_x(), H.size_y());
	//matrix K = matrix(H.size_x(), H.size_y());



	

	matrix G_ = matrix(2, 2);

	matrix H_rows = matrix (2,2);
	
		//std::vector <double> c(H.size_x() - 1);
	//std::vector <double> s(H.size_x() - 1);
	double c = 0.;
	double s = 0.;

	//for (unsigned int i = 0; i < H.size_x() - 1; i++)
	//{
		//Wird zu 2*2 Matrix

		c = R(R.size_x()-2, R.size_y() - 2) / (sqrt((R(R.size_x() - 2, R.size_y() - 2) * R(R.size_x() - 2, R.size_y() - 2)) + (R(R.size_x() - 1, R.size_y() - 2) * R(R.size_x() - 1, R.size_y() - 2))));
		s = -1 * (R(R.size_x() - 1, R.size_y() - 2)) / (sqrt((R(R.size_x() - 2, R.size_y() - 2) * R(R.size_x() - 2, R.size_y() - 2)) + (R(R.size_x(), R.size_y() - 2) * R(R.size_x() - 1, R.size_y() - 2))));

		//std::cout << std::endl;
		//std::cout << c << "   " << s << std::endl;


		G_(0, 0) = c;		G_(0, 1) = (-1 * s);	  G_(1, 0) = s;		G_(1, 1) = c;
		//std::cout << "G_ matrix in above loop is: " << std::endl;
		//G_[i].print();

		std::cout << std::endl;


		for (unsigned int k = R.size_x() - 2; k <= R.size_x() - 1; k++)
		{
			for (unsigned int l = R.size_y() - 2; l <= R.size_y() - 1; l++)
			{
				H_rows(k , l) = R(k, l);
			}
		}

		//std::cout << "Matrix H_rows is" << std::endl;
		//H_rows[i].print();
		//Anwendung von 2*2 Matrix nur auf die beiden (2) entsprechenden Zeilen von H

		H_rows = mat_mat_prod(G_, H_rows);
		//std::cout << "Matrix H_rows after multiplication with G_ is" << std::endl;
		//H_rows[i].print();


		for (unsigned int k = R.size_x() - 2; k <= R.size_x() - 1; k++)
		{
			for (unsigned int l = R.size_y() - 2; l <= R.size_y() - 1; l++)
			{
				R(k, l) = H_rows(k , l);
			}
		}
		///std::cout << "Matrix R is" << std::endl;
		///R.print();
		
	std::cout << "R_Q matrix calculated inside qr_givens_rotation func is: " << std::endl;
	R.print();
	std::cout << std::endl;
	return (R);
}
*/