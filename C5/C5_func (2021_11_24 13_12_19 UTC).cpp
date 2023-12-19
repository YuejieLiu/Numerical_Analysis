#include<iostream>
#include<vector>
#include"C4_func.h"
using namespace std;

double norm2(vector<double>& b) {
	//计算向量的二范数
	int n = b.size(), i; double norm=0;
	for (i = 0; i < n; i++) {
		norm += b[i] * b[i];
	}
	norm = sqrt(norm);
	return norm;
}

void Conjugate_Gradient(vector<vector<double>>& A, vector<double>&b,vector<double>&x) {
	//共轭梯度法求解方程组Ax=b,A对称正定
	int n = b.size(), k=0,i; double rho, rho_hat,alpha,beta;
	vector<double>x0(n, 0), r(n, 0), w(n, 0), p(n, 0);
	r = A_multiply_b(A, x0);
	for (i = 0; i < n; i++) {
		r[i] = b[i] - r[i];
	}
	rho = b_multiply_c(r, r);
	while (sqrt(rho) > 1e-7 * norm2(b) && k < 2 * n) {
		k++;
		if (k == 1)p = r;
		else {
			beta = rho / rho_hat;
			for (i = 0; i < n; i++) {
				p[i] = r[i] + beta * p[i];
			}
		}
		w = A_multiply_b(A, p);
		alpha = rho / b_multiply_c(p, w);
		for (i = 0; i < n; i++) {
			x0[i] = x0[i] + alpha *p[i];
		}
		for (i = 0; i < n; i++) {
			r[i] = r[i] - alpha * w[i];
		}
		rho_hat = rho;
		rho = b_multiply_c(r, r);
	}
	x = x0;
	cout << "共轭梯度法在限制条件下计算了" << k << "步" << endl;
}

double b_multiply_c_matrix(vector < vector<double>>& A, vector < vector<double>>& B) {
	//矩阵拉直后做内积
	int n = A.size(), i, j; double x = 0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			x += A[i][j] * B[i][j];
		}
	}
	return x;
}

