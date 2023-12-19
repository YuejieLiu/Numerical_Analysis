#include<iostream>
#include<vector>
#include"C7_func.h" 
using namespace std;

int main() {
	int n = 100, m = 100, i;
	vector<vector<double>>A(n, vector<double>(n, 0));
	for (i = 0; i < n-1; i++) {
		A[i][i] = 4.0;
		A[i][i + 1] = 1.0;
		A[i + 1][i] = 1.0;
	}
	A[n - 1][n - 1] = 4.0;
	Jacobi_Threshold(A);
	cout << "过关Jacobi计算" << n << "阶矩阵的全部特征值为" << endl;
	for (i = 0; i < n; i++) {
		cout << A[i][i] << " ";
		if ((i + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	//第二题
	vector<double>x(m, 2), y(m, -1); double max, min;
	y[0] = 0;
	max = bisect(x, y, m);
	min = bisect(x, y, 1);
	cout << "二分法求得最大特征值为" << max << endl;
	cout << "二分法求得最小特征值为" << min << endl;
}