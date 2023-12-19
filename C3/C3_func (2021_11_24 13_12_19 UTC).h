#pragma once
#include<iostream>
#include<vector>
using namespace std;

void Householder(vector<double>& x, vector<double>& v, double& beta);//����x��householder�任��v=x-norm2(x)e1,beta=vt*v
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b);//��������������˵ľ���
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta);//����һ��Householder�任�˾���A;
void QR(vector<vector<double>>& A, vector<double>& d);//����A��QR�ֽⲢ������A�У�d�д洢����beta_k
vector<double>LS(vector<vector<double>>& A, vector<double>& b);//����LS����ĽⲢ����