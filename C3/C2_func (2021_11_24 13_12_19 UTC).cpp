#include<iostream>
#include<vector>
#include<random>
using namespace std;
#include"C1.h"

vector<double>good_randVec(int n)
{	//生成随机向量(每调用一次函数使用生成一个固定的随机向量）
	// 由于我们希望引擎和分布对象保持状态，因此应该将他们定义为
	// static 的，从而每次调用都生成新的数
	static default_random_engine e;
	static uniform_int_distribution<unsigned> u(0, 10000);//0到100之间均匀分布的随机数
	vector<double> ret;
	for (int i = 0; i < n; i++)
		ret.push_back((double)u(e)/100);
	return ret;
}
double norm1_vector(vector<double>& b) {//计算向量b的1范数
	int n = b.size(); double a = 0;
	for (int i = 0; i < n; i++) {
		a += fabs(b[i]);
	}
	return a;
}

double norm_infinite_vector(vector<double>& b) {
	//计算向量b的无穷范数
	double norm = fabs(b[0]); int n = b.size();
	for (int i = 1; i < n; i++)
		if (fabs(b[i]) > norm)norm = fabs(b[i]);
	return norm;
}

void transpose(vector<vector<double>>& A,vector<vector<double>>&At){
	//计算转置矩阵
	int m = A.size(), n = A[0].size(); 
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			At[j][i] = A[i][j];
}

vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b) {
	//计算矩阵左乘向量
	int m = A.size(), n = A[0].size(); vector<double>c(m, 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c[i] += A[i][j] * b[j];
	return c;
}
double b_multiply_c(vector<double>& b, vector<double>& c) {
	//计算向量内积
	int n = b.size(); double a = 0;
	for (int i = 0; i < n; i++)
		a += b[i] * c[i];
	return a;
}
double norm1_matrix(vector<vector<double>>& A) {
	//计算方阵A的1范数
	int n = A.size(), i, j, k = 1,p; vector<double>x(n, 1 / (double)n), w(n, 0), v(n, 0), z(n, 0);
	vector<vector<double>>At(n, vector<double>(n)); double norm;
	while (k == 1) {
	again:for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				w[i] += A[i][j] * x[j];//w=Bx
			if (w[i] == 0) {//不可导，重新选取初值
				x.clear(); w.clear(); w.resize(n);
				x = good_randVec(n);
				norm = norm1_vector(x);
				for (p = 0; p < n; p++) {
					x[p]=x[p] / norm;
					w[p] = 0;
				}
				goto again;
			}
			else if (w[i] > 0)v[i] = 1;//v=sign(w)
			else v[i] = -1;
		}
		transpose(A,At);
		z = A_multiply_b(At, v);//计算z=At*v
		if (norm_infinite_vector(z) <= b_multiply_c(z, x))//找到了极大值
			k = 0;
		else {//未找到极大值，重新选取x
			for (j = 0; j < n; j++)
				if (fabs(z[j]) == norm_infinite_vector(z))break;
			for (i = 0; i < n; i++) {
				w[i] = 0;//注意要重置w
				if (i == j)x[i] = 1;
				else x[i] = 0;
			}
			k = 1;
		}
	}
	return norm1_vector(w);
}
double norm_infinite_inverseTmatrix(vector<vector<double>>& A) {
	//计算逆矩阵的无穷范数
	int n = A.size(), i, j,k=1,p,s=0; vector<double>x(n,0), w(n,0), v(n, 0),z(n,0);
	vector<vector<double>>At(n, vector<double>(n, 0)), B(n, vector<double>(n, 0));
	double norm;
	x[n-1] = 1;
	while (k == 1) {
		again:transpose(A, At);
		w = x;
		Gauss_Col_Solve(At, w);//计算At*w=x得到w
		for (i = 0; i < n; i++) {
			if (w[i] == 0) {//不可导，重新选取初值
				x.clear();
				x = good_randVec(n);
				for (j = 0; j < n; j++)x[j] = x[j] + s;
				s++;
				norm = norm1_vector(x);
				for (p = 0; p < n; p++) {
					x[p] = x[p] / norm;
					w[p] = 0;
				}
				goto again;
			}
		    else if (w[i] > 0)v[i] = 1;//v=sign(w)
			else v[i] = -1;
		}
		z = v;
		B = A;
		Gauss_Col_Solve(B, z);//解Az=v得到z
		if (norm_infinite_vector(z)-b_multiply_c(z, x)<=0) {//找到了极大值
			k = 0;
			return norm1_vector(w);
		}
		else {//未找到极大值，重新选取x
			for (i = 0; i < n; i++) {
				x[i] = 0;
			}
			for (j = 0; j < n; j++) {
				if (fabs(z[j])==norm_infinite_vector(z))break;
			}
			x[j] = 1;
			k = 1;	
		}
	}
}
