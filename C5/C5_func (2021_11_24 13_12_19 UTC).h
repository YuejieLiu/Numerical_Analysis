#pragma once
#include<iostream>
#include<vector>
using namespace std;

double norm2(vector<double>& b);//���������Ķ�����
void Conjugate_Gradient(vector<vector<double>>& A, vector<double>& b, vector<double>& x);
//�����ݶȷ���ⷽ����Ax=b,A�Գ�����
double b_multiply_c_matrix(vector < vector<double>>& A, vector < vector<double>>& B);//������ֱ�����ڻ�