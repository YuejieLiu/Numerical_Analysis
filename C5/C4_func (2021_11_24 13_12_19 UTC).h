#pragma once
#include<iostream>
#include<vector>
using namespace std;

double norm1_vector(vector<double>& b);//计算向量b的1范数
double b_multiply_c(vector<double>& b, vector<double>& c);//计算向量内积
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//计算矩阵左乘向量
void Jacobi(vector<vector<double>>& A, vector<double>& b, vector<double>& x,int&m);
//用Jacobi迭代法计算方程组Ax=b的解,返回迭代次数，输入x为初值
void GS(vector<vector<double>>& A, vector<double>& b, vector<double>& x,int&m);
//用G-S迭代法计算方程组Ax=b的解，返回迭代次数，输入x为初值
void SOR(vector<vector<double>>& A, vector<double>& b, vector<double>& x, double omega,int&m);
//用SOR方法求解方程组Ax=b,返回迭代次数，输入x为初值
double Frobenius(vector < vector<double>>& A);//计算矩阵的Frobenius范数