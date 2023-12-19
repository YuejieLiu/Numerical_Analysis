#pragma once
#pragma once
#include<iostream>
#include<vector>
using namespace std;

vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//计算矩阵左乘向量
double power_method_polynomial(vector<double>& a);
// 用幂法求首一多项式方程的模最大根, 输入系数矩阵a = { a_n - 1...a0 }
double norm_infinite_vector(vector<double>& b);//计算向量b的无穷范数
void Householder(vector<double>& x, vector<double>& v, double& beta);
//计算x的householder变换，v=x-norm2(x)e1,beta=vt*v
void transpose(vector<vector<double>>& A, vector<vector<double>>& At);//计算转置矩阵
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b);//计算两个向量相乘的矩阵
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta);
//计算一个Householder变换乘矩阵A;
void Hessenberg(vector<vector<double>>& A, vector<vector<double>>& Q);//计算A的上Hessenberg分解并储存在A中
void QR_iteration(vector<vector<double>>& H, vector<vector<double>>& P);
//双重步位移的QR迭代一次，变换矩阵记录在P中(这里要求n>=3)
int Implicit_QR(vector<vector<double>>& A);//隐式QR算法，最终结果存在A中，成功返回0，否则返回1
void Eigen_Value(vector<vector<double>>& A);//对一个拟上三角阵，输出其特征值
	