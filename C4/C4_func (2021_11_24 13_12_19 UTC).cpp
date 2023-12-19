#include<iostream>
#include<vector>
using namespace std;

double norm1_vector(vector<double>& b) {//计算向量b的1范数
	int n = b.size(); double a = 0;
	for (int i = 0; i < n; i++) {
		a += fabs(b[i]);
	}
	return a;
}
double b_multiply_c(vector<double>& b, vector<double>& c) {
	//计算向量内积
	int n = b.size(); double a = 0;
	for (int i = 0; i < n; i++)
		a += b[i] * c[i];
	return a;
}
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b) {
	//计算矩阵左乘向量
	int m = A.size(), n = A[0].size(); vector<double>c(m, 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c[i] += A[i][j] * b[j];
	return c;
}
void Jacobi(vector<vector<double>>&A, vector<double>&b, vector<double>& x,int&m) {
	//用Jacobi迭代法计算方程组Ax=b的解,返回迭代次数，输入x为初值
	int n = A.size(), i, j; double error;
	vector<double>inverseD(n,0),g(b),y(x),z(n,0); vector<vector<double>>B(A);
	for (i = 0; i < n; i++) {
		inverseD[i] = 1 / A[i][i];
		B[i][i] = 0.0;
		for (j = 0; j < n; j++) {
			B[i][j] = -B[i][j];
		}
	}
	for (i = 0; i < n; i++) {//计算迭代矩阵B和常数项g
		for (j = 0; j < n; j++) {
			B[i][j] *= inverseD[i];
		}
		g[i] *= inverseD[i];
	}
	m = 1;
	y = A_multiply_b(B, x);
	for (i = 0; i < n; i++)
		y[i] += g[i];
	for (i = 0; i < n; i++)z[i] = x[i] - y[i];
	error = norm1_vector(z);
	while (error >= 1e-7) {//迭代，终止条件为前后两项一范数只差小于1e-7
		x = y; m++;
		y = A_multiply_b(B, x);
		for (i = 0; i < n; i++)
			y[i] += g[i];
		for (i = 0; i < n; i++)z[i] = x[i] - y[i];
		error = norm1_vector(z);
	}
	x = y;
}

 void GS(vector<vector<double>>& A, vector<double>& b, vector<double>& x,int&m) {
	//用G-S迭代法计算方程组Ax=b的解，返回迭代次数，输入x为初值
	int n = A.size(), i, j; double error;
	vector<double>inverseD(n, 0), g(b),y(x),z(n,0); vector<vector<double>>B(A);
	for (i = 0; i < n; i++) {
		inverseD[i] = 1 / A[i][i];
		B[i][i] = 0.0;
		for (j = 0; j < n; j++) {
			B[i][j] = -B[i][j];
		}
	}
	for (i = 0; i < n; i++) {//计算迭代矩阵B和常数项g
		for (j = 0; j < n; j++) {
			B[i][j] *= inverseD[i];
		}
		g[i] *= inverseD[i];
	}
	for (i = 0; i < n; i++) {
		x[i] = b_multiply_c(B[i], x);
		x[i] = x[i] + g[i];
	}
	m = 1;
	for (i = 0; i < n; i++)z[i] = x[i] - y[i];
	error = norm1_vector(z);
	while (error >= 1e-7) {//迭代，终止条件为前后两项一范数只差小于1e-7
		y = x; m++;
		for (i = 0; i < n; i++) {
			x[i] = b_multiply_c(B[i], x);
			x[i] = x[i] + g[i];
		}
		for (i = 0; i < n; i++)z[i] = x[i] - y[i];
		error = norm1_vector(z);
	}
}

void SOR(vector<vector<double>>& A, vector<double>& b, vector<double>& x,double omega,int&m) {
	//用SOR方法求解方程组Ax=b,返回迭代次数，输入x为初值
	int n = A.size(), i, j; double error;
	vector<double>inverseD(n, 0), g(b), y(x),z(n,0); vector<vector<double>>B(A);
	for (i = 0; i < n; i++) {
		inverseD[i] = 1 / A[i][i];
		B[i][i] = 0.0;
		for (j = 0; j < n; j++) {
			B[i][j] = -B[i][j];
		}
	}
	for (i = 0; i < n; i++) {//计算迭代矩阵B和常数项g
		for (j = 0; j < n; j++) {
			B[i][j] *= inverseD[i];
		}
		g[i] *= inverseD[i];
	}
	m = 1;
	for (i = 0; i < n; i++) {
		x[i] = omega * (b_multiply_c(B[i], x) + g[i]) + (1 - omega) * x[i];
	}
	for (i = 0; i < n; i++)z[i] = x[i] - y[i];
	error = norm1_vector(z);
	while (error >= 1e-7) {//迭代，终止条件为前后两项一范数只差小于1e-7
		y = x; m++;
		for (i = 0; i < n; i++) {
			x[i] = omega * (b_multiply_c(B[i], x) + g[i]) + (1 - omega) * x[i];
		}
		for (i = 0; i < n; i++)z[i] = x[i] - y[i];
		error = norm1_vector(z);
	}
}
double Frobenius(vector < vector<double>>& A) {
	//计算矩阵的Frobenius范数
	int n = A.size(), i, j; double norm=0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			norm += A[i][j] * A[i][j];
		}
	}
	norm = sqrt(norm);
	return norm;
}

