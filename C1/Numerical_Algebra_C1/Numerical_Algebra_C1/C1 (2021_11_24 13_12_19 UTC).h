#pragma once
#include<iostream>
#include<vector>
using namespace std;

void forward_sub(vector< vector<double> >& L, vector<double>& b);//前代法，解下三角形二维方阵Lx=b,解存储在b中
void backward_sub(vector< vector<double> >& U, vector<double>& b);//回代法，解下三角形二维方阵Ux=b,解存储在b中
void Gauss(vector< vector<double> >& A);//不选主元的Gauss消去,左下角存储L右上角存储U
void Gauss_Solve(vector<vector<double>>& A, vector<double>& b);//用不选主元的高斯消去法解方程组,解存储在b中
void Gauss_Col(vector< vector<double> >& A, vector<double>& u);
//列主元高斯消去.A是要分解的矩阵，L存在左下角，U存在左上角。用一个向量u存储置换矩阵P
void Gauss_Col_Solve(vector< vector<double> >& A, vector<double>& b);//用列主元高斯消去解方程组，解存于b中
void Square_root(vector< vector<double> >& A);//平方根法。输入方阵A，将L存储在下三角，L的转置存储在上三角
void Square_root_Solve(vector<vector<double>>& A, vector<double>& b);//用平方根法计算正定对称系数矩阵方程，结果存于b中
void Improved_Square_root(vector< vector<double> >& A);//改进的平方根法A=LDL^t,L储存在左下角，D储存在对角线，L^t储存在右上角
void Improved_Square_root_Solve(vector<vector<double>>& A, vector<double>& b);//用改进的平方根法计算正定对称系数矩阵方程，结果存于b中