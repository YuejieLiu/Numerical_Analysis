#pragma once
#include<iostream>
#include<vector>
using namespace std;

vector<double>good_randVec(int n);//生成随机向量
double norm1_vector(vector<double>& b);//计算向量b的1范数
double norm_infinite_vector(vector<double>& b);//计算向量b的无穷范数
void transpose(vector<vector<double>>& A, vector<vector<double>>& At);//计算转置矩阵
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//计算矩阵左乘向量
double b_multiply_c(vector<double>& b, vector<double>& c);//计算向量内积
double norm1_matrix(vector<vector<double>>& A);//计算方阵A的1范数
double norm_infinite_inverseTmatrix(vector<vector<double>>& A);//计算逆矩阵的无穷范数
void Gauss_Col_Solve(vector< vector<double> >& A, vector<double>& b);//用列主元高斯消去解方程组，解存于b中