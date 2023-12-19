#pragma once
#include<iostream>
#include<vector>
using namespace std;

void forward_sub(vector< vector<double> >& L, vector<double>& b);//ǰ���������������ζ�ά����Lx=b,��洢��b��
void backward_sub(vector< vector<double> >& U, vector<double>& b);//�ش��������������ζ�ά����Ux=b,��洢��b��
void Gauss(vector< vector<double> >& A);//��ѡ��Ԫ��Gauss��ȥ,���½Ǵ洢L���ϽǴ洢U
void Gauss_Solve(vector<vector<double>>& A, vector<double>& b);//�ò�ѡ��Ԫ�ĸ�˹��ȥ���ⷽ����,��洢��b��
void Gauss_Col(vector< vector<double> >& A, vector<double>& u);
//����Ԫ��˹��ȥ.A��Ҫ�ֽ�ľ���L�������½ǣ�U�������Ͻǡ���һ������u�洢�û�����P
void Gauss_Col_Solve(vector< vector<double> >& A, vector<double>& b);//������Ԫ��˹��ȥ�ⷽ���飬�����b��
void Square_root(vector< vector<double> >& A);//ƽ�����������뷽��A����L�洢�������ǣ�L��ת�ô洢��������
void Square_root_Solve(vector<vector<double>>& A, vector<double>& b);//��ƽ���������������Գ�ϵ�����󷽳̣��������b��
void Improved_Square_root(vector< vector<double> >& A);//�Ľ���ƽ������A=LDL^t,L���������½ǣ�D�����ڶԽ��ߣ�L^t���������Ͻ�
void Improved_Square_root_Solve(vector<vector<double>>& A, vector<double>& b);//�øĽ���ƽ���������������Գ�ϵ�����󷽳̣��������b��