#pragma once
#include<iostream>
#include<vector>
using namespace std;

vector<double>good_randVec(int n);//�����������
double norm1_vector(vector<double>& b);//��������b��1����
double norm_infinite_vector(vector<double>& b);//��������b�������
void transpose(vector<vector<double>>& A, vector<vector<double>>& At);//����ת�þ���
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//��������������
double b_multiply_c(vector<double>& b, vector<double>& c);//���������ڻ�
double norm1_matrix(vector<vector<double>>& A);//���㷽��A��1����
double norm_infinite_inverseTmatrix(vector<vector<double>>& A);//���������������
void Gauss_Col_Solve(vector< vector<double> >& A, vector<double>& b);//������Ԫ��˹��ȥ�ⷽ���飬�����b��