#pragma once
#include<iostream>
#include<vector>
using namespace std;

double EA(vector<vector<double>>& A);//����ԳƷ���ķǶԽǷ���
void Jacobi_Threshold(vector<vector<double>>& A);
//�ù���Jacobi������ʵ�Գƾ����ȫ������ֵ������ȡsigma=n
double tri_norm(vector<double>& x, vector<double>& y);//�������ԽǾ���������, x�ǶԽ�Ԫ��y�ǴζԽ�Ԫ
int signal(vector<double>& x, vector<double>& y, double miu);
//�������ԽǾ����Ӧ����ʽ�ı����,x�ǶԽ�Ԫ��y�ǴζԽ�Ԫ������Ϊ0,miu��ָ����
double bisect(vector<double>& x, vector<double>& y, int m);
	//�ö��ַ�����ʵ�Գ����ԽǾ���ָ������ֵ,����ָ������Ϊ1e-12
	//x�ǶԽ�Ԫ��y�ǴζԽ�Ԫ������Ϊ0,m�ǵ�m������ֵ