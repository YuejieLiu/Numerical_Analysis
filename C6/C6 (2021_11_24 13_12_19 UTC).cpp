#include<iostream>
#include<vector>
#include"C6_func.h" 
using namespace std;

int main() {
	vector<double>a = { 1,-5,3 }, b = { 0,-3,-1 }, c = { 101,208.01,10891.01,9802.08,79108.9,-99902,790,-1000 };
	double x;
	//第一题
	x = power_method_polynomial(a);
	cout << "第一个多项式方程的模最大根为" << x << endl;
	x = power_method_polynomial(b);
	cout << "第二个多项式方程的模最大根为" << x << endl;
	x = power_method_polynomial(c);
	cout << "第三个多项式方程的模最大根为" << x << endl << endl;
	//第二题第三问
	vector<vector<double>>A1 = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,0.9},{6.1,4.9,3.5,6.2} },
		A2 = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,1.0},{6.1,4.9,3.5,6.2} },
		A3 = { {9.1,3.0,2.6,4.0},{4.2,5.3,4.7,1.6},{3.2,1.7,9.4,1.1},{6.1,4.9,3.5,6.2} };
	int judge,i,j,n=41;
	judge = Implicit_QR(A1);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "x=0.9,隐式QR算法计算得到的拟上三角方阵为" << endl;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			cout << A1[i][j] << " ";
		}
		cout << endl;
	}
	cout << "其特征值为" << endl;
	Eigen_Value(A1);
	cout << endl;

	judge = Implicit_QR(A2);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "x=1.0,隐式QR算法计算得到的拟上三角方阵为" << endl;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			cout << A2[i][j] << " ";
		}
		cout << endl;
	}
	cout << "其特征值为" << endl;
	Eigen_Value(A2);
	cout << endl;

	judge = Implicit_QR(A3);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "x=1.1,隐式QR算法计算得到的拟上三角方阵为" << endl;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			cout << A3[i][j] << " ";
		}
		cout << endl;
	}
	cout << "其特征值为" << endl;
	Eigen_Value(A3);
	cout << endl;
	//第二题第二问
	vector<vector<double>>C(n, vector<double>(n, 0));
	C[0][n - 1] = -1;//生成友阵
	C[0][n - 4] = -1;
	for (i = 1; i < n; i++) {
		C[i][i - 1] = 1;
	}
	judge = Implicit_QR(C);
	if (judge) { cout << "error" << endl; return 0; }
	cout << "41阶方程特征值为" << endl;
	Eigen_Value(C);
	cout << endl;
}