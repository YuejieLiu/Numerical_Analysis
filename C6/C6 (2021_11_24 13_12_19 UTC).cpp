#include<iostream>
#include<vector>
#include"C6_func.h" 
using namespace std;

int main() {
	vector<double>a = { 1,-5,3 }, b = { 0,-3,-1 }, c = { 101,208.01,10891.01,9802.08,79108.9,-99902,790,-1000 };
	double x;
	//��һ��
	x = power_method_polynomial(a);
	cout << "��һ������ʽ���̵�ģ����Ϊ" << x << endl;
	x = power_method_polynomial(b);
	cout << "�ڶ�������ʽ���̵�ģ����Ϊ" << x << endl;
	x = power_method_polynomial(c);
	cout << "����������ʽ���̵�ģ����Ϊ" << x << endl << endl;
	//�ڶ��������
	vector<vector<double>>A1 = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,0.9},{6.1,4.9,3.5,6.2} },
		A2 = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,1.0},{6.1,4.9,3.5,6.2} },
		A3 = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,1.1},{6.1,4.9,3.5,6.2} };
	int judge,i,j,n=41;
	judge = Implicit_QR(A1);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "x=0.9,��ʽQR�㷨����õ����������Ƿ���Ϊ" << endl;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			cout << A1[i][j] << " ";
		}
		cout << endl;
	}
	cout << "������ֵΪ" << endl;
	Eigen_Value(A1);
	cout << endl;

	judge = Implicit_QR(A2);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "x=1.0,��ʽQR�㷨����õ����������Ƿ���Ϊ" << endl;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			cout << A2[i][j] << " ";
		}
		cout << endl;
	}
	cout << "������ֵΪ" << endl;
	Eigen_Value(A2);
	cout << endl;

	judge = Implicit_QR(A3);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "x=1.1,��ʽQR�㷨����õ����������Ƿ���Ϊ" << endl;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			cout << A3[i][j] << " ";
		}
		cout << endl;
	}
	cout << "������ֵΪ" << endl;
	Eigen_Value(A3);
	cout << endl;
	//�ڶ���ڶ���
	vector<vector<double>>C(n, vector<double>(n, 0));
	C[0][n - 1] = -1;//��������
	C[0][n - 4] = -1;
	for (i = 1; i < n; i++) {
		C[i][i - 1] = 1;
	}
	judge = Implicit_QR(C);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "41�׷�������ֵΪ" << endl;
	Eigen_Value(C);
	cout << endl;
}