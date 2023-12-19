#pragma once
#pragma once
#include<iostream>
#include<vector>
using namespace std;

vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//��������������
double power_method_polynomial(vector<double>& a);
// ���ݷ�����һ����ʽ���̵�ģ����, ����ϵ������a = { a_n - 1...a0 }
double norm_infinite_vector(vector<double>& b);//��������b�������
void Householder(vector<double>& x, vector<double>& v, double& beta);
//����x��householder�任��v=x-norm2(x)e1,beta=vt*v
void transpose(vector<vector<double>>& A, vector<vector<double>>& At);//����ת�þ���
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b);//��������������˵ľ���
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta);
//����һ��Householder�任�˾���A;
void Hessenberg(vector<vector<double>>& A, vector<vector<double>>& Q);//����A����Hessenberg�ֽⲢ������A��
void QR_iteration(vector<vector<double>>& H, vector<vector<double>>& P);
//˫�ز�λ�Ƶ�QR����һ�Σ��任�����¼��P��(����Ҫ��n>=3)
int Implicit_QR(vector<vector<double>>& A);//��ʽQR�㷨�����ս������A�У��ɹ�����0�����򷵻�1
void Eigen_Value(vector<vector<double>>& A);//��һ���������������������ֵ
	