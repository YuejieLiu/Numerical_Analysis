#include<iostream>
#include<vector>
#include <algorithm>
#include"C1.h"
using namespace std;

int main() {
	int a1 = 84 - 1, a2 = 100 - 1, a3 = 5 - 1, k; 
	vector< vector<double> >A1(a1+1, vector<double>(a1+1,0)), A2(a2+1,vector<double>(a2+1,0)), A3(a3+1,vector<double>(a3+1));
	vector<double>b1(a1+1,0), b2(a2+1,0), b3(a3+1,0);
	//生成第一题的84阶方程组A1
	for (k = 0; k <a1 ; k++) {
		A1[k][k] = 6;
		A1[k][k + 1] = 1;
		A1[k + 1][k] = 8;
	}
	A1[a1][a1] = 6;
	vector<vector<double>>B1(A1);//为第一题创建A1的副本
	//生成第一题的b1
	for (k = 1; k < a1; k++)
		b1[k] = 15;
	b1[0] = 7; b1[a1] = 14;
	vector<double>c1(b1);//为第一题创建b1的副本
	//生成第二题第一问的系数矩阵A2
	for (k = 0; k < a2; k++) {
		A2[k][k] = 10;
		A2[k][k + 1]=A2[k + 1][k] = 1;
	}
	A2[a2][a2] = 10;
	vector<vector<double>>B2(A2), C2(A2),D2(A2);//A2副本
	//生成第二问的随机向量b2
	for (k = 0; k <= a2; k++) {
		b2[k]=(double)(k + 1);
	}
	random_shuffle(b2.begin(), b2.end());
	vector<double>c2(b2), d2(b2), e2(b2);//b2的三个副本
	//生成第二题第二问的系数矩阵A3和b3
	for (int p = 0; p <= a3; p++) {//两个整数相除只保留整数，因此用double
		for (int q = 0; q<= a3; q++) {
			A3[p][q]=1/((double)p+ (double)q+1);
			b3[p] +=1/((double)p+ (double)q+1);
		}
	}
	vector<vector<double>>B3(A3), C3(A3), D3(A3);
	vector<double>c3(b3), d3(b3), e3(b3);//副本
	cout << "第一题,不选主元" << endl;
	Gauss_Solve(A1, b1);
	for (k = 0; k <= a1; k++) {
		cout << b1[k] << " ";
		if ((k+1 )%10==0 )cout << endl;
	}
	cout << endl;
	cout << "第一题，列主元" << endl;
	Gauss_Col_Solve(B1, c1);
	for (k = 0; k <= a1; k++) {
		cout << c1[k] << " ";
		if ((k+1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第二题，平方根法，第一个方程" << endl;
	cout << "系数矩阵b为" << endl;
	for (k = 0; k <= a2; k++) {
		cout << b2[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout<<"该矩阵的解为"<<endl;
	Square_root_Solve(A2, b2);
	for (k = 0; k <= a2; k++) {
		cout << b2[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第二题，改进的平方根法，第一个方程" << endl;
	cout << "系数矩阵b为" << endl;
	for (k = 0; k <= a2; k++) {
		cout << c2[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout<<"该矩阵的解为"<<endl;
	Improved_Square_root_Solve(B2, c2);
	for (k = 0; k <= a2; k++) {
		cout << c2[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第二题，平方根法，第二个方程" << endl;
	Square_root_Solve(A3, b3);
	cout << "由于矩阵病态，以上输出无效的Akk,以下是无效的解"<<endl;
	for (k = 0; k <= a3; k++) {
		cout << b3[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第二题，改进的平方根法，第二个方程" << endl;
	Improved_Square_root_Solve(B3, c3);
	for (k = 0; k <= a3; k++) {
		cout << c3[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第三题，不选主元Gauss消元，第一个方程" << endl;
	Gauss_Solve(C2,d2);
	for (k = 0; k <= a2; k++) {
		cout << d2[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第三题，列主元Gauss消元，第一个方程" << endl;
	Gauss_Col_Solve(D2, e2);
	for (k = 0; k <= a2; k++) {
		cout << e2[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第三题，不选主元Gauss消元，第二个方程" << endl;
	Gauss_Solve(C3, d3);
	for (k = 0; k <= a3; k++) {
		cout << d3[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	cout << "第三题，列主元Gauss消元，第二个方程" << endl;
	Gauss_Col_Solve(D3, e3);
	for (k = 0; k <= a3; k++) {
		cout << e3[k] << " ";
		if ((k + 1) % 10 == 0)cout << endl;
	}
}