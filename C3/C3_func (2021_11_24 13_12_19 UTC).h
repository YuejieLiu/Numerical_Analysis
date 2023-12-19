#pragma once
#include<iostream>
#include<vector>
using namespace std;

void Householder(vector<double>& x, vector<double>& v, double& beta);//计算x的householder变换，v=x-norm2(x)e1,beta=vt*v
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b);//计算两个向量相乘的矩阵
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta);//计算一个Householder变换乘矩阵A;
void QR(vector<vector<double>>& A, vector<double>& d);//计算A的QR分解并储存在A中，d中存储的是beta_k
vector<double>LS(vector<vector<double>>& A, vector<double>& b);//计算LS问题的解并返回