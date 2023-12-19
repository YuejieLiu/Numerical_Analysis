#include<iostream>
#include<vector>
using namespace std;

double EA(vector<vector<double>>& A) {
	//����ԳƷ���ķǶԽǷ���
	int n = A.size(), i, j; double norm=0;
	for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
			norm += A[i][j] * A[i][j];
	norm =2.0*norm;
	norm = sqrt(norm);
	return norm;
}

void Jacobi_Threshold(vector<vector<double>>& A) {
	//�ù���Jacobi������ʵ�Գƾ����ȫ������ֵ������ȡsigma=n
	int n = A.size(), j,p,q; double delta,tau,t,s,c;
	vector<vector<double>>B(A);
	delta = EA(A);
	while(delta>1e-20){//ÿѭ��һ��Ϊһ��ɨ��
		for (p = 0; p < n - 1; p++) {
			for (q = p+1; q < n; q++) {
				if (fabs(A[p][q]) <= delta)continue;//ȷ��p,q
				//ȷ��c,s
				tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
				if (tau > 0)t = 1.0 / (tau + sqrt(1.0 + tau * tau));
				else t = 1.0 / (tau - sqrt(1.0 + tau * tau));
				c = 1.0 / sqrt(1.0 + t * t);
				s = t * c;
				//����һ�ε���
				for (j = 0; j < n; j++) {
					B[p][j] = c * A[p][j] - s * A[q][j];
					B[q][j] = s * A[p][j] + c * A[q][j];
				}
				A = B;
				for (j = 0; j < n; j++) {
					B[j][p] = c * A[j][p] - s * A[j][q];
					B[j][q] = s * A[j][p] + c * A[j][q];
				}
				A = B;
			}
		}
		delta = delta / n;//�ı���ֵ
	}
}

double tri_norm(vector<double>& x, vector<double>& y) {
	//�������ԽǾ���������, x�ǶԽ�Ԫ��y�ǴζԽ�Ԫ
	int n = x.size(), i, j; double norm, max;
	max = fabs(x[0]) + fabs(y[1]);
	for (i = 1; i < n - 1; i++) {
		norm = fabs(x[i]) + fabs(y[i + 1]) + fabs(y[i]);
		if (norm > max)max = norm;
	}
	norm = fabs(x[n - 1]) + fabs(y[n - 1]);
	if (norm > max)return norm;
	else return max;
}

int signal(vector<double>& x, vector<double>& y,double miu) {
	//�������ԽǾ����Ӧ����ʽ�ı����,x�ǶԽ�Ԫ��y�ǴζԽ�Ԫ������Ϊ0,miu��ָ����
	int s = 0,n=x.size(),k; double q = x[0] - miu;
	for (k = 0; k < n; k++) {
		if (q < 0)s++;
		if (k < n - 1) {
			if (q == 0) q = fabs(y[k + 1]) * 1e-10;
			q = x[k + 1] - miu - y[k + 1] * y[k + 1] / q;
		}
	}
	return s;
}

double bisect(vector<double>& x, vector<double>& y,int m) {
	//�ö��ַ�����ʵ�Գ����ԽǾ���ָ������ֵ,����ָ������Ϊ1e-12
	//x�ǶԽ�Ԫ��y�ǴζԽ�Ԫ������Ϊ0,m�ǵ�m������ֵ
	int n = x.size(), i,s; double l, u, r;
	u = tri_norm(x, y);
	l = 0.0 - u;
	while ((u - l) > 1e-12) {
		r = (l + u) / 2;
		s = signal(x, y, r);
		if (s >= m)u = r;
		else l = r;
	}
	return (l + u) / 2;
}
	