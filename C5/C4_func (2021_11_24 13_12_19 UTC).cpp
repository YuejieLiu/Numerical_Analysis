#include<iostream>
#include<vector>
using namespace std;

double norm1_vector(vector<double>& b) {//��������b��1����
	int n = b.size(); double a = 0;
	for (int i = 0; i < n; i++) {
		a += fabs(b[i]);
	}
	return a;
}
double b_multiply_c(vector<double>& b, vector<double>& c) {
	//���������ڻ�
	int n = b.size(); double a = 0;
	for (int i = 0; i < n; i++)
		a += b[i] * c[i];
	return a;
}
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b) {
	//��������������
	int m = A.size(), n = A[0].size(); vector<double>c(m, 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c[i] += A[i][j] * b[j];
	return c;
}
void Jacobi(vector<vector<double>>&A, vector<double>&b, vector<double>& x,int&m) {
	//��Jacobi���������㷽����Ax=b�Ľ�,���ص�������������xΪ��ֵ
	int n = A.size(), i, j; double error;
	vector<double>inverseD(n,0),g(b),y(x),z(n,0); vector<vector<double>>B(A);
	for (i = 0; i < n; i++) {
		inverseD[i] = 1 / A[i][i];
		B[i][i] = 0.0;
		for (j = 0; j < n; j++) {
			B[i][j] = -B[i][j];
		}
	}
	for (i = 0; i < n; i++) {//�����������B�ͳ�����g
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
	while (error >= 1e-7) {//��������ֹ����Ϊǰ������һ����ֻ��С��1e-7
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
	//��G-S���������㷽����Ax=b�Ľ⣬���ص�������������xΪ��ֵ
	int n = A.size(), i, j; double error;
	vector<double>inverseD(n, 0), g(b),y(x),z(n,0); vector<vector<double>>B(A);
	for (i = 0; i < n; i++) {
		inverseD[i] = 1 / A[i][i];
		B[i][i] = 0.0;
		for (j = 0; j < n; j++) {
			B[i][j] = -B[i][j];
		}
	}
	for (i = 0; i < n; i++) {//�����������B�ͳ�����g
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
	while (error >= 1e-7) {//��������ֹ����Ϊǰ������һ����ֻ��С��1e-7
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
	//��SOR������ⷽ����Ax=b,���ص�������������xΪ��ֵ
	int n = A.size(), i, j; double error;
	vector<double>inverseD(n, 0), g(b), y(x),z(n,0); vector<vector<double>>B(A);
	for (i = 0; i < n; i++) {
		inverseD[i] = 1 / A[i][i];
		B[i][i] = 0.0;
		for (j = 0; j < n; j++) {
			B[i][j] = -B[i][j];
		}
	}
	for (i = 0; i < n; i++) {//�����������B�ͳ�����g
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
	while (error >= 1e-7) {//��������ֹ����Ϊǰ������һ����ֻ��С��1e-7
		y = x; m++;
		for (i = 0; i < n; i++) {
			x[i] = omega * (b_multiply_c(B[i], x) + g[i]) + (1 - omega) * x[i];
		}
		for (i = 0; i < n; i++)z[i] = x[i] - y[i];
		error = norm1_vector(z);
	}
}
double Frobenius(vector < vector<double>>& A) {
	//��������Frobenius����
	int n = A.size(), i, j; double norm=0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			norm += A[i][j] * A[i][j];
		}
	}
	norm = sqrt(norm);
	return norm;
}

