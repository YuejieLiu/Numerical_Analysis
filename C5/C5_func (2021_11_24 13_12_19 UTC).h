#pragma once
#include<iostream>
#include<vector>
using namespace std;

double norm2(vector<double>& b);//计算向量的二范数
void Conjugate_Gradient(vector<vector<double>>& A, vector<double>& b, vector<double>& x);
//共轭梯度法求解方程组Ax=b,A对称正定
double b_multiply_c_matrix(vector < vector<double>>& A, vector < vector<double>>& B);//矩阵拉直后做内积