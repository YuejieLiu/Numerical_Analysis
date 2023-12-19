#include<vector>
#include"SVD_func.h" 
using namespace std;

int main() {
	vector<vector<double>> A = { {1.0000000000,4.9176000000,1.0000000000,3.4720000000,0.9980000000,1.0000000000,7.0000000000,4.0000000000,42.0000000000,3.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,5.0208000000,1.0000000000,3.5310000000,1.5000000000,2.0000000000,7.0000000000,4.0000000000,62.0000000000,1.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,4.5429000000,1.0000000000,2.2750000000,1.1750000000,1.0000000000,6.0000000000,3.0000000000,40.0000000000,2.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,4.5573000000,1.0000000000,4.0500000000,1.2320000000,1.0000000000,6.0000000000,3.0000000000,54.0000000000,4.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,5.0597000000,1.0000000000,4.4550000000,1.1210000000,1.0000000000,6.0000000000,3.0000000000,42.0000000000,3.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,3.8910000000,1.0000000000,4.4550000000,0.9880000000,1.0000000000,6.0000000000,3.0000000000,56.0000000000,2.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,5.8980000000,1.0000000000,5.8500000000,1.2400000000,1.0000000000,7.0000000000,3.0000000000,51.0000000000,2.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,5.6039000000,1.0000000000,9.5200000000,1.5010000000,0.0000000000,6.0000000000,3.0000000000,32.0000000000,1.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,15.4202000000,2.5000000000,9.8000000000,3.4200000000,2.0000000000,10.0000000000,5.0000000000,42.0000000000,2.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,14.4598000000,2.5000000000,12.8000000000,3.0000000000,2.0000000000,9.0000000000,5.0000000000,14.0000000000,4.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,5.8282000000,1.0000000000,6.4350000000,1.2250000000,2.0000000000,6.0000000000,3.0000000000,32.0000000000,1.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,5.3003000000,1.0000000000,4.9883000000,1.5520000000,1.0000000000,6.0000000000,3.0000000000,30.0000000000,1.0000000000,2.0000000000,0.0000000000},
		{1.0000000000,6.2712000000,1.0000000000,5.5200000000,0.9750000000,1.0000000000,5.0000000000,2.0000000000,30.0000000000,1.0000000000,2.0000000000,0.0000000000},
		{1.0000000000,5.9592000000,1.0000000000,6.6660000000,1.1210000000,2.0000000000,6.0000000000,3.0000000000,32.0000000000,2.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,5.0500000000,1.0000000000,5.0000000000,1.0200000000,0.0000000000,5.0000000000,2.0000000000,46.0000000000,4.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,5.6039000000,1.0000000000,9.5200000000,1.5010000000,0.0000000000,6.0000000000,3.0000000000,32.0000000000,1.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,8.2464000000,1.5000000000,5.1500000000,1.6640000000,2.0000000000,8.0000000000,4.0000000000,50.0000000000,4.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,6.6969000000,1.5000000000,6.0920000000,1.4880000000,1.5000000000,7.0000000000,3.0000000000,22.0000000000,1.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,7.7841000000,1.5000000000,7.1020000000,1.3760000000,1.0000000000,6.0000000000,3.0000000000,17.0000000000,2.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,9.0384000000,1.0000000000,7.8000000000,1.5000000000,1.5000000000,7.0000000000,3.0000000000,23.0000000000,3.0000000000,3.0000000000,0.0000000000},
		{1.0000000000,5.9894000000,1.0000000000,5.5200000000,1.2560000000,2.0000000000,6.0000000000,3.0000000000,40.0000000000,4.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,7.5422000000,1.5000000000,4.0000000000,1.6900000000,1.0000000000,6.0000000000,3.0000000000,22.0000000000,1.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,8.7951000000,1.5000000000,9.8900000000,1.8200000000,2.0000000000,8.0000000000,4.0000000000,50.0000000000,1.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,6.0931000000,1.5000000000,6.7265000000,1.6520000000,1.0000000000,6.0000000000,3.0000000000,44.0000000000,4.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,8.3607000000,1.5000000000,9.1500000000,1.7770000000,2.0000000000,8.0000000000,4.0000000000,48.0000000000,1.0000000000,1.0000000000,1.0000000000},
		{1.0000000000,8.1400000000,1.0000000000,8.0000000000,1.5040000000,2.0000000000,7.0000000000,3.0000000000,3.0000000000,1.0000000000,3.0000000000,0.0000000000},
		{1.0000000000,9.1416000000,1.5000000000,7.3262000000,1.8310000000,1.5000000000,8.0000000000,4.0000000000,31.0000000000,4.0000000000,1.0000000000,0.0000000000},
		{1.0000000000,12.0000000000,1.5000000000,5.0000000000,1.2000000000,2.0000000000,6.0000000000,3.0000000000,30.0000000000,3.0000000000,1.0000000000,1.0000000000} };
	int m = A.size(), n = A[0].size(), judge, i,j; double norm;
	vector<vector<double>>U(m, vector<double>(m, 0)),Ut(U),V(n, vector<double>(n, 0)),Vt(V),B(V),A1(A);
	judge=SVD(A,U,V,B);
	if (judge)cout << "分解错误" << endl;
	else {
		vector<vector<double>>Q(Ut), P(Vt),UI(Ut),VI(Vt);
		for (i = 0; i < m; i++) {
			UI[i][i] = 1.0;
		}
		for (i = 0; i < n; i++) {
			VI[i][i] = 1.0;
		}
		transpose(U, Ut);
		Q = multiply(U, Ut);
		Q=matrix_minus(Q, UI);
		norm = max_mod(Q);
		cout << "误差项U*Ut-I的模最大元素为" << norm << endl;
		transpose(V, Vt);
		P = multiply(V, Vt);
		P = matrix_minus(P, VI);
		norm = max_mod(P);
		cout << "误差项V*Vt-I的模最大元素为" << norm << endl;
		vector<vector<double>>sigma(m, vector<double>(n, 0));
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				sigma[i][j] = B[i][j];
			}
		}
		sigma = multiply(U, sigma);
		sigma = multiply(sigma, Vt);
		sigma = matrix_minus(sigma, A1);
		norm = max_mod(sigma);
		cout << "误差项U*sigma*Vt-A的模最大元素为" << norm << endl;
	}
}