#include<vector>
#include<random>
using namespace std;
#include"C1.h"
#include"C2_func.h"

void Householder(vector<double>& x, vector<double>& v, double&beta) {
	//����x��householder�任��v=x-norm2(x)e1,beta=vt*v
	//������beta�Ĵ洢�ռ�ı�beta��ֵ�����дΪ&beta
	int n = x.size(), i; double iota = norm_infinite_vector(x), sigma = 0, alpha;
	for (i = 0; i < n; i++) {
		x[i] /= iota;
	}
	for (i = 1; i < n; i++) {
		sigma += x[i] * x[i];
	}
	v = x;//������Ȼֻ����2-n������������v1�����¸�ֵ������ֱ�Ӹ�������
	if (sigma == 0)beta = 0;
	else {
		alpha = sqrt(x[0] * x[0] + sigma);
		if (x[0] <= 0)v[0] = x[0] - alpha;
		else v[0] = -sigma / (x[0] + alpha);
	}
	beta = 2 * v[0] * v[0] / (sigma + v[0] * v[0]);
	for (i = 1; i < n; i++) {
		v[i] = v[i] / v[0];
	}
	v[0] = 1;
	
}

vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b){
	//��������������˵ľ���
	int i,j,m = a.size(),n = b.size(); vector <vector<double>>A(m,vector<double>(n));
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			A[i][j] = a[i] * b[j];
	return A;
}
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta) {
	//����һ��Householder�任�˾���A;
	int m = A.size(), n = A[0].size();
	vector<vector<double>>At(n, vector<double>(m)),B(m, vector<double>(n));
	transpose(A, At);
	vector<double>w = A_multiply_b(At, v);
	for (int i = 0; i < n; i++)w[i] *= beta;
	B=a_mutiply_bt(v, w);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			B[i][j] = A[i][j] - B[i][j];
	return B;
}
void QR(vector<vector<double>>& A, vector<double>& d) {
	//����A��QR�ֽⲢ������A�У�d�д洢����beta_k
	int m = A.size(), n = A[0].size(), i, j, k;
	vector < vector<double>>A0; vector<double>v,a0; double beta;
	for (j = 0; j < n; j++) {
		if (j < m) {
			a0.resize(m - j); v.resize(m - j);
			for (i = j; i < m; i++)a0[i - j] = A[i][j];
			Householder(a0, v, beta);//�����j+1�е�householder�任
			//����Ϊ����household�任���A�����½Ǿ���
			A0.resize(m - j);
			for (i = 0; i < m - j; i++)A0[i].resize(n - j);
			for (i = 0; i < m - j; i++)//�����½Ǿ����Ƶ�A0��
				for (k = 0; k < n - j; k++)
					A0[i][k] = A[i+j][k+j];
			A0=HouseMultiply(A0, v, beta);
			for (i = 0; i < m - j; i++)//��household���ù���A0���ƻ�ȥ
				for (k = 0; k < n - j; k++)
					A[i + j][k + j] = A0[i][k];
			d[j] = beta;
			for (i = 1; i < m - j; i++)A[i + j][j] = v[i];//��v�浽A��
			A0.clear(); v.clear(); a0.clear();
		}
	}
}

vector<double>LS(vector<vector<double>>& A, vector<double>& b) {
	//����LS����ĽⲢ����
	int m = A.size(), n = A[0].size(),i,j,k;
	vector<double>d(n),vn(m,0),v(m,0),c1; 
	vector<vector<double>>H(m,vector<double>(m,0)),B(m, vector<double>(m, 0)),Q1(m,vector<double>(n)),Q1t(n, vector<double>(m));
	QR(A, d);//����A��QR�ֽ�
	//����Ϊ����Q
	for (i = n; i < m; i++)vn[i] = A[i][n-1];
	vn[n-1] = 1;
	for (i = 0; i < m; i++)H[i][i] = 1;
	B = a_mutiply_bt(vn, vn);
	for (i = n - 1; i < m; i++)
		for (j = n - 1; j < m; j++)
			H[i][j] -= d[n - 1] * B[i][j];//����Hn
	for (k = n - 2; k >= 0; k--) {
		for (i = k+1; i < m; i++)v[i] = A[i][k];
		v[k] = 1;
		H = HouseMultiply(H, v, d[k]);
	}//ѭ������Hk�˻����յõ�Q
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			Q1[i][j] = H[i][j];//ȡ�ֿ�
	transpose(Q1, Q1t);
	c1 = A_multiply_b(Q1t, b);
	vector<vector<double>>R(n, vector<double>(n, 0));
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)
			R[i][j] = A[i][j];//����R
	backward_sub(R, c1);
	return c1;
}

