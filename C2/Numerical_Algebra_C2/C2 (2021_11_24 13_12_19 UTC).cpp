#include<iostream>
#include<vector>
#include <algorithm>
#include"C2_func.h"
#include"C1.h"
using namespace std;

int main() {
	int n, i, j, k; double kappa; vector<vector<double>>A,An,B; vector<double>x,b,r,delta;
	double niu, gamma, beta, miu,rho,rho1,a,a1;
	for (n = 5; n <= 20; n++) {//����5-20�׵�ϣ�����ؾ��󲢼�����������
		A.resize(n);
		for (int p = 0; p < n; p++) {//����ϣ�����ؾ���
			A[p].resize(n);
			for (int q = 0; q < n; q++) {
				A[p][q] = 1.0 / (p + q + 1);
			}
		}
		a = norm1_matrix(A);
		a1= norm_infinite_inverseTmatrix(A);
		kappa = a * a1;
		//����ϣ�����ؾ���Գƣ����������һ����
		cout << "n=" << n << " " <<"A��һ����Ϊ"<<a<< " ������Ϊ" << kappa << endl;
	}
	for (n = 5; n <=30; n++) {
		x = good_randVec(n);//�����������
		An.resize(n); B.resize(n);//���ɷ���An
		for (i = 0; i < n; i++) {
			An[i].resize(n); B[i].resize(n);
			An[i][i] =An[i][n-1] = 1;
			for (j = i + 1; j < n - 1; j++)An[i][j] = 0;
			for (j = 0; j < i; j++)An[i][j] = -1;
		}
		niu = norm_infinite_inverseTmatrix(An);//����An����������
		transpose(An,B);
		miu = norm1_matrix(B);//����An�������
		b.resize(n); delta.resize(n);
		b = A_multiply_b(An, x);//����b=An*x
		beta = norm_infinite_vector(b);
		vector<double>b1(b);//b�ĸ���
		vector<vector<double>>A1(An);//An�ĸ���
		Gauss_Col_Solve(An, b);//�ⷽ��An*x=b,����b�У���ʱAn�ѷֽ⣩
		r= A_multiply_b(A1, b);
		for (i = 0; i < n; i++) {
			r[i] = b1[i] - r[i]; 
		}//r�Ǿ����������
		gamma= norm_infinite_vector(r);
		rho = niu * miu * gamma / beta;//rhoΪ��ľ��ȹ���
		for (i = 0; i < n; i++) {
			delta[i] = x[i] - b[i];
		}
		rho1 = norm_infinite_vector(delta) / norm_infinite_vector(x);//��ʵ��������
		cout << "n=" << n << " �������xΪ"<<endl;
		for (i = 0; i < n; i++) {
			cout << x[i] << " ";
		}
		cout << endl << "����õ��Ľ�Ϊ" << endl;
		for (i = 0; i < n; i++) {
			cout << b[i] << " ";
		}
		cout << endl<<"��ľ��ȹ���ֵΪ"<<rho<<" ��ʵ��������Ϊ" << rho1 << endl << endl;

	}
}
	