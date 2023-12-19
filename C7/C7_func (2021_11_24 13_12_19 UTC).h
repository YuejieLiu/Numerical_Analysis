#pragma once
#include<iostream>
#include<vector>
using namespace std;

double EA(vector<vector<double>>& A);//计算对称方阵的非对角范数
void Jacobi_Threshold(vector<vector<double>>& A);
//用过关Jacobi方法求实对称矩阵的全部特征值，这里取sigma=n
double tri_norm(vector<double>& x, vector<double>& y);//计算三对角矩阵的无穷范数, x是对角元，y是次对角元
int signal(vector<double>& x, vector<double>& y, double miu);
//计算三对角矩阵对应多项式的变号数,x是对角元，y是次对角元，首项为0,miu是指定点
double bisect(vector<double>& x, vector<double>& y, int m);
	//用二分法计算实对称三对角矩阵指定特征值,这里指定精度为1e-12
	//x是对角元，y是次对角元，首项为0,m是第m个特征值