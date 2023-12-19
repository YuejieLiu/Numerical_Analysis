#pragma once
#include<iostream>
#include<vector>
using namespace std;

double norm_infinite_vector(vector<double>& b);//计算向量b的无穷范数
void Householder(vector<double>& x, vector<double>& v, double& beta);
	//计算x的householder变换，v=x-norm2(x)e1,beta=vt*v
	//这里在beta的存储空间改变beta的值，因此写为&beta
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//计算矩阵左乘向量
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b);//计算两个向量相乘的矩阵
void bi_diagonal_house(vector<vector<double>>& A, vector<double>& b, vector<double>& c);
	//用householder变换法将矩阵阵A化为二对角矩阵,A的行数m大于列数n
	//计算后A中二对角元素下方是左乘Householder的v，右上角是右乘householder变换的v
	//b中是左乘的beta,至多n个元素，c中是右乘的beta，有n-2个元素
void transpose(vector<vector<double>>& A, vector<vector<double>>& At);//计算转置矩阵
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta);
//计算一个Householder变换乘矩阵A;
void accumulateUV(vector<vector<double>>& A, vector<double>& b, vector<double>& c, vector<vector<double>>& U, vector<vector<double>>& V);
	//对已经二对角化的矩阵计算U和V
	//b中是左乘的beta,有至多n个元素，c中是右乘的beta，有n-2个元素
void Givens(double a, double b, double& c, double& s);//计算givens变换
void Wilkinson_SVD(vector<vector<double>>& A, vector<double>& Qc, vector<double>& Qs, vector<double>& Pc, vector<double>& Ps);
	//对已经二对角化的带Wilkinson位移的SVD迭代（不可约）(n>=3)，
	//其中Qc Qs Ps Pc是n-1维向量，记录givens变换
void Wilkinson2(vector<vector<double>>& A, vector<double>& Qc, vector<double>& Qs, vector<double>& Pc, vector<double>& Ps);
	//对已经二对角化的带Wilkinson位移的SVD迭代（不可约）(n=2)，
	//其中Qc Qs Ps Pc是1维向量，记录givens变换
double matrix_infinite_norm(vector<vector<double>>& A);//计算矩阵的无穷范数
void print(vector<vector<double>>& A);//打印矩阵的所有元素
vector<vector<double>>multiply(vector<vector<double>>& A, vector<vector<double>>& B);//计算矩阵AB
vector<vector<double>>matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B);//计算矩阵相减A-B
double max_mod(vector<vector<double>>& A);//计算矩阵的最大模元素
void Givens2(double a, double b, double& c, double& s);//计算givens变换(第一个分量为0）
int SVD(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V, vector<vector<double>>& B);
//求m*n阶矩阵A的奇异值分解（m>=n)(B是0矩阵）