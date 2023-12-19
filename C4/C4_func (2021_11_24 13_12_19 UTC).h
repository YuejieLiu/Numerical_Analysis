#pragma once
#include<iostream>
#include<vector>
using namespace std;

double norm1_vector(vector<double>& b);//��������b��1����
double b_multiply_c(vector<double>& b, vector<double>& c);//���������ڻ�
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//��������������
void Jacobi(vector<vector<double>>& A, vector<double>& b, vector<double>& x,int&m);
//��Jacobi���������㷽����Ax=b�Ľ�,���ص�������������xΪ��ֵ
void GS(vector<vector<double>>& A, vector<double>& b, vector<double>& x,int&m);
//��G-S���������㷽����Ax=b�Ľ⣬���ص�������������xΪ��ֵ
void SOR(vector<vector<double>>& A, vector<double>& b, vector<double>& x, double omega,int&m);
//��SOR������ⷽ����Ax=b,���ص�������������xΪ��ֵ
double Frobenius(vector < vector<double>>& A);//��������Frobenius����