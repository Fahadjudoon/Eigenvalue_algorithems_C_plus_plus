#include <iostream>
#include <vector>
#include <assert.h>
#include "matrix.h"
#include <windows.h>
#include "vectors.h"
#include "linear_algebra.h"
#include <chrono>
#include<thread>
#include "matrix.cpp" 
using namespace std;

// O-Notation nachlesen (Landau-Symbolik, Algorithmik, Laufzeitanalyse)

int main()
{
	using namespace std::literals::chrono_literals;

	unsigned int e = 0;
	//matrix A = mat_five();
	//matrix A = mat_ten();
	//matrix A = mat_fifteen();
	//matrix A = mat_twenty();
	//matrix A = mat_thirty();
	//matrix A = mat_fifty();
	//matrix A = mat_seventy_five(); 
	matrix A = mat_hundred();
	//matrix A = mat_hundred_fifty();
	//matrix A = mat_two_hundred();   
	//std::cout << "Matrix A is " << std::endl;
	//A.print();
	//std::cout << std::endl;
	
	auto start = std::chrono::high_resolution_clock::now();
	matrix id = identity_mat(A.size_x(), A.size_y());
	matrix H = matrix(A.size_x(), A.size_y());
	matrix H_QR = matrix(A.size_x(), A.size_y());
	matrix R = matrix(A.size_x(), A.size_y());
	matrix Q = matrix(A.size_x(), A.size_y());


	std::vector <matrix> H_;
	H_.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		H_.push_back(matrix(A.size_x(), A.size_y()));
	}
	
	std::vector <matrix> I;
	I.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		I.push_back(matrix(A.size_x() - (i + 1), A.size_y() - (i + 1)));
		I[i] = identity_mat(A.size_x() - (i + 1), A.size_y() - (i + 1));
	}

	std::vector <matrix> P;
	P.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		P.push_back(matrix(A.size_x() - (i + 1), A.size_y() - (i + 1)));
	}

	std::vector <matrix> P_;
	P_.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		P_.push_back(matrix(A.size_x(), A.size_y()));
		P_[i] = id;
	}

	std::vector <vectors> x;
	x.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		x.push_back(vectors(A.size_x() - (i + 1)));
	}

	std::vector <vectors> u;
	u.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		u.push_back(vectors(A.size_x() - (i + 1)));
	}

	std::vector <vectors> uT;
	uT.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		uT.push_back(vectors(A.size_x() - (i + 1)));
	}

	std::vector <vectors> u_vec;
	u_vec.reserve(A.size_x() - 2);
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		u_vec.push_back(vectors(A.size_x() - (i + 1)));
	}
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		u_vec[i](0) = 1.;
	}

	/*
	for (unsigned int i = 0; i < A.size_x() - 2; i++)
	{
		if (i == 0)
		{
			//Picking elements from Matrix A for first x vector
			for (unsigned int l = 1; l < A.size_x(); l++)
			{
				x[0](l - 1) = A(l, 0);
			}
		}
		u[i] = x[i] + signum(x[i](0)) * vec_norm(x[i]) * u_vec[i];
		//u[i].print();

		P[i] = I[i] - ((2 / vec_trans_vec_prod(vec_trans(u[i]), u[i])) * (vec_vec_trans_prod(u[i], vec_trans(u[i]))));
		//P[i].print();

		P_[i] = id;

		for (unsigned int j = i + 1; j < P_[i].size_x(); j++)
		{
			for (unsigned int k = i + 1; k < P_[i].size_x(); k++)
			{
				//std::cout << j << "," << k << P_[i].size_x() << id.size_x() << std::endl;
				P_[i](j, k) = P[i](j - (1 + i), k - (1 + i));
				//P_[i].print();
			}
		}
		//P_[i].print();
		if (i == 0)
		{
			H_[i] = mat_mat_prod(mat_mat_prod(P_[i], A), P_[i]);
			//std::cout << i << "For H_0" << std::endl;
			//H_[i].print();
		}
		else
		{
			//std::cout << i << "For H_" << std::endl;
			H_[i] = mat_mat_prod(mat_mat_prod(P_[i], H_[i - 1]), P_[i]);
			//H_[i].print();
		}

		if (i < A.size_x() - 3)
		{
			unsigned int z = 1 + i;
			for (unsigned int j = 2 + i; j < A.size_x(); j++)
			{
				x[i + 1](j - (2 + i)) = H_[i](j, z);
			}
		}
		H = H_[i];

	}
	std::cout << "Hessenberg Matrix Form is" << std::endl;
	H.print();

	R = H;

	//auto start = std::chrono::high_resolution_clock::now();

	//Givens Rotation QR Algorithm

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

		c = R(i, i) / (sqrt((R(i, i) * R(i, i)) + (R(i + 1, i) * R(i + 1, i))));
		s = (-1 * R(i + 1, i)) / (sqrt((R(i, i) * R(i, i)) + (R(i + 1, i) * R(i + 1, i))));

		//std::cout << std::endl;
		//std::cout << c << "   " << s << std::endl;


		G_[i](0, 0) = c;		G_[i](0, 1) = (-1 * s);	  G_[i](1, 0) = s;		G_[i](1, 1) = c;
		//std::cout << "Given matrix is: " << std::endl;
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
		//std::cout << "Matrix at the end: (inside first loop) of Givens Rotation is" << std::endl;
		//R.print();

	}
	//std::cout << std::endl;
	//std::cout << "Matrix after first loop of Givens Rotation is" << std::endl;
	//R.print();
	matrix Hh = matrix(H.size_x(), H.size_y());
	Hh = R;

	std::vector <matrix> Hh_sub;
	Hh_sub.reserve(H.size_x() - 1);
	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		Hh_sub.push_back(matrix(i + 2, 2));
	}

	for (unsigned int i = 0; i < H.size_x() - 1 ; i++)
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
	std::cout << std::endl;
	std::cout << "Matrix after second loop of Givens Rotation is: " << std::endl;
	Hh.print();
/*
	matrix U = identity_mat(H.size_x(), H.size_y());
	//U = H;

	std::vector <matrix> U_sub;
	U_sub.reserve(H.size_x() - 1);
	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		U_sub.push_back(matrix(H.size_x(), 2));
	}


	for (unsigned int i = 0; i < H.size_x() - 1; i++)
	{
		c = H(i, i) / (sqrt((H(i, i) * H(i, i)) + (H(i + 1, i) * H(i + 1, i))));
		s = -1 * (H(i + 1, i)) / (sqrt((H(i, i) * H(i, i)) + (H(i + 1, i) * H(i + 1, i))));
		G_[i](0, 0) = c;	G_[i](0, 1) = (-1 * s);	G_[i](1, 0) = s;		G_[i](1, 1) = c;
		G_T[i] = transpose(G_[i]);


		for (unsigned int k = 0; k < H.size_x(); k++)
		{
			for (unsigned int l = i; l < i + 2; l++)
			{
				U_sub[i](k, l - i) = U(k, l);
			}
		}
		//std::cout << "U_sub matrix is: " << std::endl;
		//U_sub[i].print();

		U_sub[i] = mat_mat_prod(U_sub[i], G_T[i]);

		//std::cout << "U_sub matrix after multiplication is: " << std::endl;
		//U_sub[i].print();

		for (unsigned int k = 0; k < H.size_x(); k++)
		{
			for (unsigned int l = i; l < i + 2; l++)
			{
				U(k, l) = U_sub[i](k, l - i);
			}

		}
	}
	///std::cout << "U matrix is: " << std::endl;
	///U.print();
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> duration = end - start;
	cout << duration.count() << "s " << std::endl;                
	*/

	//auto start = std::chrono::high_resolution_clock::now();

	//Alg 4.4: The Hessenberg QR alg with Rayleigh quotient
	
	matrix Hr_save = matrix(H.size_x(), H.size_y());
	matrix R_Q = matrix(H.size_x(), H.size_y());

	//int k = 0;
	unsigned int f = 0;
	double sigma = 0.;
	double h_m_m = 0.;
	double h_m_m_1 = 0.;
	
	std::cout << std::endl;
	std::cout << "eigen values with Rayleigh quotient are: " << std::endl;

	for (unsigned int m = H.size_y() - 1; m >= 2; m--)
	{
		
		//std::cout << " m is" << m << std::endl;
		std::vector <matrix> Hr;
		Hr.reserve(50);
		for (unsigned int i = 0; i < 50; i++)
		{
			Hr.push_back(matrix(H.size_x() - f, H.size_y() - f));
		}
		matrix ident = identity_mat(H.size_x() - f, H.size_y() - f);
		matrix pp = matrix(H.size_x() - f, H.size_y() - f);
		f++;
		if (m == H.size_y() - 1)
		{
			Hr[0] = H;
		}
		else
		{
			for (unsigned int i = 0; i <= m ; i++)
			{
				for (unsigned int j = 0; j <= m ; j++)
				{
					Hr[0](i, j) = Hr_save(i, j);
				}
			}
		}
		
		unsigned int k = 0;
		do
		{
			k++;
			//std::cout << "Hr[k] at the beginning of loop is: " << std::endl;
			//std::cout << std::endl;
			//Hr[k - 1].print();

			h_m_m = Hr[k - 1](m, m);
			//std::cout << "h_m_m inside is: " << h_m_m << std::endl;
			//std::cout << std::endl;

			pp = Hr[k - 1] - h_m_m * ident;		//QR
			//std::cout << "matrix: (Hr[k - 1] - Lamb * I)  for qr_givens_rotation func is: " << std::endl;
			//pp.print();
			//std::cout << std::endl;

			R_Q = qr_givens_rotation(pp);	//RQ

			Hr[k] = R_Q + h_m_m * ident;

			//std::cout << "Hr[k] after (R*Q + Lamb * I) operation is: " << std::endl;
			//std::cout << std::endl;
			//Hr[k].print();

			h_m_m_1 = Hr[k](m, m - 1);
			//std::cout << std::endl;
			//std::cout << "h_m_m_1 inside at end of loop is: " << h_m_m_1 << std::endl;
			//std::cout << std::endl;

		} while (abs(h_m_m_1) > 1e-12);

		//std::cout << std::endl;
		std::cout << "h_m_m(eig) after loop is: " << Hr[k](m, m) << std::endl;
		//std::cout << std::endl;
		Hr[k](m, m) = 0.;
		
		Hr_save = Hr[k];
		//e = k;
	
		//std::cout << "Hr[0] at the end of loop after comparision is: " << std::endl;
		//std::cout << std::endl;
		//Hr[0].print();
	}
	
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> duration = end - start;
	cout << duration.count() << "s " << std::endl;
	/*
	//auto start = std::chrono::high_resolution_clock::now();
	//eig(A);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float> duration = end - start;
	cout << duration.count() << "s " << std::endl;
	//compute_all_eigenvectors(A);
	*/
	system("pause");
	return 0;
}


