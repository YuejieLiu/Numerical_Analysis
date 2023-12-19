#include<iostream>
#include<vector>
#include <algorithm>
#include"C2_func.h"
#include"C1.h"
using namespace std;

int main() {
	int n, i, j, k; double kappa; vector<vector<double>>A,An,B; vector<double>x,b,r,delta;
	double niu, gamma, beta, miu,rho,rho1,a,a1;
	for (n = 5; n <= 20; n++) {//生成5-20阶的希尔伯特矩阵并计算其条件数
		A.resize(n);
		for (int p = 0; p < n; p++) {//生成希尔伯特矩阵
			A[p].resize(n);
			for (int q = 0; q < n; q++) {
				A[p][q] = 1.0 / (p + q + 1);
			}
		}
		a = norm1_matrix(A);
		a1= norm_infinite_inverseTmatrix(A);
		kappa = a * a1;
		//由于希尔伯特矩阵对称，无穷范数等于一范数
		cout << "n=" << n << " " <<"A的一范数为"<<a<< " 条件数为" << kappa << endl;
	}
	for (n = 5; n <=30; n++) {
		x = good_randVec(n);//生成随机向量
		An.resize(n); B.resize(n);//生成方阵An
		for (i = 0; i < n; i++) {
			An[i].resize(n); B[i].resize(n);
			An[i][i] =An[i][n-1] = 1;
			for (j = i + 1; j < n - 1; j++)An[i][j] = 0;
			for (j = 0; j < i; j++)An[i][j] = -1;
		}
		niu = norm_infinite_inverseTmatrix(An);//计算An的逆的无穷范数
		transpose(An,B);
		miu = norm1_matrix(B);//计算An的无穷范数
		b.resize(n); delta.resize(n);
		b = A_multiply_b(An, x);//计算b=An*x
		beta = norm_infinite_vector(b);
		vector<double>b1(b);//b的副本
		vector<vector<double>>A1(An);//An的副本
		Gauss_Col_Solve(An, b);//解方程An*x=b,解在b中（这时An已分解）
		r= A_multiply_b(A1, b);
		for (i = 0; i < n; i++) {
			r[i] = b1[i] - r[i]; 
		}//r是绝对误差向量
		gamma= norm_infinite_vector(r);
		rho = niu * miu * gamma / beta;//rho为解的精度估计
		for (i = 0; i < n; i++) {
			delta[i] = x[i] - b[i];
		}
		rho1 = norm_infinite_vector(delta) / norm_infinite_vector(x);//真实地相对误差
		cout << "n=" << n << " 随机向量x为"<<endl;
		for (i = 0; i < n; i++) {
			cout << x[i] << " ";
		}
		cout << endl << "计算得到的解为" << endl;
		for (i = 0; i < n; i++) {
			cout << b[i] << " ";
		}
		cout << endl<<"解的精度估计值为"<<rho<<" 真实地相对误差为" << rho1 << endl << endl;

	}
}
	