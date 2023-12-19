#pragma once
#include<iostream>
#include<vector>
using namespace std;

double norm_infinite_vector(vector<double>& b);//��������b�������
void Householder(vector<double>& x, vector<double>& v, double& beta);
	//����x��householder�任��v=x-norm2(x)e1,beta=vt*v
	//������beta�Ĵ洢�ռ�ı�beta��ֵ�����дΪ&beta
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b);//��������������
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b);//��������������˵ľ���
void bi_diagonal_house(vector<vector<double>>& A, vector<double>& b, vector<double>& c);
	//��householder�任����������A��Ϊ���ԽǾ���,A������m��������n
	//�����A�ж��Խ�Ԫ���·������Householder��v�����Ͻ����ҳ�householder�任��v
	//b������˵�beta,����n��Ԫ�أ�c�����ҳ˵�beta����n-2��Ԫ��
void transpose(vector<vector<double>>& A, vector<vector<double>>& At);//����ת�þ���
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta);
//����һ��Householder�任�˾���A;
void accumulateUV(vector<vector<double>>& A, vector<double>& b, vector<double>& c, vector<vector<double>>& U, vector<vector<double>>& V);
	//���Ѿ����Խǻ��ľ������U��V
	//b������˵�beta,������n��Ԫ�أ�c�����ҳ˵�beta����n-2��Ԫ��
void Givens(double a, double b, double& c, double& s);//����givens�任
void Wilkinson_SVD(vector<vector<double>>& A, vector<double>& Qc, vector<double>& Qs, vector<double>& Pc, vector<double>& Ps);
	//���Ѿ����Խǻ��Ĵ�Wilkinsonλ�Ƶ�SVD����������Լ��(n>=3)��
	//����Qc Qs Ps Pc��n-1ά��������¼givens�任
void Wilkinson2(vector<vector<double>>& A, vector<double>& Qc, vector<double>& Qs, vector<double>& Pc, vector<double>& Ps);
	//���Ѿ����Խǻ��Ĵ�Wilkinsonλ�Ƶ�SVD����������Լ��(n=2)��
	//����Qc Qs Ps Pc��1ά��������¼givens�任
double matrix_infinite_norm(vector<vector<double>>& A);//�������������
void print(vector<vector<double>>& A);//��ӡ���������Ԫ��
vector<vector<double>>multiply(vector<vector<double>>& A, vector<vector<double>>& B);//�������AB
vector<vector<double>>matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B);//����������A-B
double max_mod(vector<vector<double>>& A);//�����������ģԪ��
void Givens2(double a, double b, double& c, double& s);//����givens�任(��һ������Ϊ0��
int SVD(vector<vector<double>>& A, vector<vector<double>>& U, vector<vector<double>>& V, vector<vector<double>>& B);
//��m*n�׾���A������ֵ�ֽ⣨m>=n)(B��0����