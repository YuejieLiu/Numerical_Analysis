#include<iostream>
#include<vector>
#include<math.h>
#include <windows.h>
#include"C4_func.h"   
using namespace std;

int main() {
	int n = 100, i, j, m, k, N =20; double epsilon = 0.0001, a = 0.5, h = 1.0 / n, omega =0.9,omega1=1.7,
		h1 = 1.0 / N;
	vector<vector<double>>A(n - 1, vector<double>(n - 1, 0)), h2f((N - 1), vector<double>(N - 1, 0)),
		inverseD(h2f), u(h2f), v(u), w(u),u1(u);
	vector<double>b(n - 1, 0), x(n - 1, 0), x1(x), y(n - 1, 0);
	double error, begin, end;
	for (i = 0; i < n - 2; i++) {//���ɵ�һ��ϵ������A
			A[i][i] = -(2 * epsilon + h);
			A[i][i + 1] = epsilon + h;
			A[i + 1][i] = epsilon;
		}
		A[n - 2][n - 2] = -(2 * epsilon + h);
		for (i = 0; i < n - 2; i++) {//���ɵ�һ�������b
			b[i] = a * h * h;
		}
		b[n - 2] = a * h * h - epsilon - h;
		for (i = 0; i < n - 1; i++) {//���㾫ȷ��y1��yn-1
			y[i] = (1 - a) / (1 - exp(-1.0 / epsilon)) * (1 - exp(-(i + 1) * h / epsilon)) + a * (i + 1) * h;
		}
		Jacobi(A, b, x, m);
		cout << "epsilonΪ" << epsilon << "ʱ��Jacobi��������þ���Ľ�Ϊ" << endl;
		for (i = 0; i < n - 1; i++) {
			cout << x[i] << " ";
			if ((i + 1) % 10 == 0)cout << endl;
		}
		cout << endl;
		cout << "ÿһ���뾫ȷ������Ϊ" << endl;
		for (i = 0; i < n - 1; i++) {
			cout << x[i] - y[i] << " ";
			if ((i + 1) % 10 == 0)cout << endl;
		}
		cout << endl;
		cout << "��������Ϊ" << m << endl;
		x = x1;
		GS(A, b, x, m);
		cout << "epsilonΪ" << epsilon << "ʱ��G-S��������þ���Ľ�Ϊ" << endl;
		for (i = 0; i < n - 1; i++) {
			cout << x[i] << " ";
			if ((i + 1) % 10 == 0)cout << endl;
		}
		cout << endl;
		cout << "ÿһ���뾫ȷ������Ϊ" << endl;
		for (i = 0; i < n - 1; i++) {
			cout << x[i] - y[i] << " ";
			if ((i + 1) % 10 == 0)cout << endl;
		}
		cout << endl;
		cout << "��������Ϊ" << m << endl;
		x = x1;
		SOR(A, b, x, omega, m);
		cout << "epsilonΪ" << epsilon << "ʱ��SOR��������þ���Ľ�Ϊ" << endl;
		for (i = 0; i < n - 1; i++) {
			cout << x[i] << " ";
			if ((i + 1) % 10 == 0)cout << endl;
		}
		cout << endl;
		cout << "ÿһ���뾫ȷ������Ϊ" << endl;
		for (i = 0; i < n - 1; i++) {
			cout << x[i] - y[i] << " ";
			if ((i + 1) % 10 == 0)cout << endl;
		}
		cout << endl;
		cout << "��������Ϊ" << m << endl;

	//�ڶ���
	for (i = 0; i < N - 1; i++) {//����ϵ������ĶԽ�ԪD����
		for (j = 0; j < N - 1; j++) {
			inverseD[i][j] = 1.0 / (4 + h1 * h1 * exp((i + 1) * (j + 1) * h1 * h1));
		}
	}
	for (i = 0; i < N - 1; i++) {//���ɳ�����
		for (j = 0; j < N - 1; j++) {
			h2f[i][j] = h1 * h1 * ((i + 1) * h1 + (j + 1) * h1)*inverseD[i][j];
		}
	}
	begin = GetTickCount64();//��ʼ��ʱ
	m = 1;//Jacobi����
	for (j = 0; j < N - 1; j++) {
		for (i = 0; i < N - 1; i++) {
			if (i == 0) {
				if (j == 0)
					v[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				else if (j == N - 2)v[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
				else v[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
			}
			else if (i == N - 2) {
				if (j == 0)
					v[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
				else if (j == N - 2)v[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j - 1]) + h2f[i][j];
				else v[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
			}
			else {
				if (j == 0)
					v[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				else if (j == N - 2)v[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
				else v[i][j] = inverseD[i][j] * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
			}
		}
	}
	for (i = 0; i < N - 1; i++) {
		for (j = 0; j < N - 1; j++) {
			w[i][j] = v[i][j] - u[i][j];
		}
	}
	error = Frobenius(w);
	while (error >= 1e-7) {//��������ֹ����Ϊǰ�����������֮��С��1e-7
		u = v; m++;
		for (j = 0; j < N - 1; j++) {
			for (i = 0; i < N - 1; i++) {
				if (i == 0) {
					if (j == 0)
						v[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
					else if (j == N - 2)v[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
					else v[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				}
				else if (i == N - 2) {
					if (j == 0)
						v[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
					else if (j == N - 2)v[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j - 1]) + h2f[i][j];
					else v[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
				}
				else {
					if (j == 0)
						v[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
					else if (j == N - 2)v[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
					else v[i][j] = inverseD[i][j] * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				}
			}
		}
		for (i = 0; i < N - 1; i++) {
			for (j = 0; j < N - 1; j++) {
				w[i][j] = v[i][j] - u[i][j];
			}
		}
		error = Frobenius(w);
	}
	end = GetTickCount64();
	cout << "NΪ" << N << "ʱ��Jacobi��������Ϊ" << m << " CPU��ʱΪ" << end - begin <<"ms"<<endl;
	u = u1; v = u1;
	begin = GetTickCount64();//��ʼ��ʱ
	m = 1;//GS����
	for (j = 0; j < N - 1; j++) {
		for (i = 0; i < N - 1; i++) {
			if (i == 0) {
				if (j == 0)
					u[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				else if (j == N - 2)u[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
				else u[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
			}
			else if (i == N - 2) {
				if (j == 0)
					u[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
				else if (j == N - 2)u[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j - 1]) + h2f[i][j];
				else u[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
			}
			else {
				if (j == 0)
					u[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				else if (j == N - 2)u[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
				else u[i][j] = inverseD[i][j] * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
			}
		}
	}
	for (i = 0; i < N - 1; i++) {
		for (j = 0; j < N - 1; j++) {
			w[i][j] = v[i][j] - u[i][j];
		}
	}
	error = Frobenius(w);
	while (error >= 1e-7) {//��������ֹ����Ϊǰ�����������֮��С��1e-7
		v = u; m++;
		for (j = 0; j < N - 1; j++) {
			for (i = 0; i < N - 1; i++) {
				if (i == 0) {
					if (j == 0)
						u[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
					else if (j == N - 2)u[i][j] = inverseD[i][j] * (2 + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
					else u[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				}
				else if (i == N - 2) {
					if (j == 0)
						u[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
					else if (j == N - 2)u[i][j] = inverseD[i][j] * (2 + u[i - 1][j] + u[i][j - 1]) + h2f[i][j];
					else u[i][j] = inverseD[i][j] * (1 + u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + h2f[i][j];
				}
				else {
					if (j == 0)
						u[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
					else if (j == N - 2)u[i][j] = inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + h2f[i][j];
					else u[i][j] = inverseD[i][j] * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j];
				}
			}
		}
		for (i = 0; i < N - 1; i++) {
			for (j = 0; j < N - 1; j++) {
				w[i][j] = v[i][j] - u[i][j];
			}
		}
		error = Frobenius(w);
	}
	end = GetTickCount64();
	cout << "NΪ" << N << "ʱ��G-S��������Ϊ" << m << " CPU��ʱΪ" << end - begin << "ms" << endl;
	u = u1; v = u1;
	begin = GetTickCount64();//��ʼ��ʱ
	m = 1;//SOR����
	for (j = 0; j < N - 1; j++) {
		for (i = 0; i < N - 1; i++) {
			if (i == 0) {
				if (j == 0)
					u[i][j] = omega1*(inverseD[i][j] * (2 + u[i + 1][j] + u[i][j + 1]) + h2f[i][j])+(1-omega1)*u[i][j];
				else if (j == N - 2)u[i][j] = omega1*(inverseD[i][j] * (2 + u[i + 1][j] + u[i][j - 1]) + h2f[i][j])+ (1 - omega1) * u[i][j];
				else u[i][j] = omega1*(inverseD[i][j] * (1 + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j])+ (1 - omega1) * u[i][j];
			}
			else if (i == N - 2) {
				if (j == 0)
					u[i][j] =omega1*(inverseD[i][j] * (2 + u[i - 1][j] + u[i][j + 1]) + h2f[i][j])+ (1 - omega1) * u[i][j];
				else if (j == N - 2)u[i][j] = omega1*(inverseD[i][j] * (2 + u[i - 1][j] + u[i][j - 1]) + h2f[i][j])+(1 - omega1) * u[i][j];
				else u[i][j] = omega1*(inverseD[i][j] * (1 + u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + h2f[i][j])+ (1 - omega1) * u[i][j];
			}
			else {
				if (j == 0)
					u[i][j] = omega1 * (inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
				else if (j == N - 2)u[i][j] = omega1 * (inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
				else u[i][j] = omega1 * (inverseD[i][j] * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
			}
		}
	}
	for (i = 0; i < N - 1; i++) {
		for (j = 0; j < N - 1; j++) {
			w[i][j] = v[i][j] - u[i][j];
		}
	}
	error = Frobenius(w);
	while (error >= 1e-7) {//��������ֹ����Ϊǰ�����������֮��С��1e-7
		v = u; m++;
		for (j = 0; j < N - 1; j++) {
			for (i = 0; i < N - 1; i++) {
				if (i == 0) {
					if (j == 0)
						u[i][j] = omega1 * (inverseD[i][j] * (2 + u[i + 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
					else if (j == N - 2)u[i][j] = omega1 * (inverseD[i][j] * (2 + u[i + 1][j] + u[i][j - 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
					else u[i][j] = omega1 * (inverseD[i][j] * (1 + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
				}
				else if (i == N - 2) {
					if (j == 0)
						u[i][j] = omega1 * (inverseD[i][j] * (2 + u[i - 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
					else if (j == N - 2)u[i][j] = omega1 * (inverseD[i][j] * (2 + u[i - 1][j] + u[i][j - 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
					else u[i][j] = omega1 * (inverseD[i][j] * (1 + u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
				}
				else {
					if (j == 0)
						u[i][j] = omega1 * (inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
					else if (j == N - 2)u[i][j] = omega1 * (inverseD[i][j] * (1 + u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
					else u[i][j] = omega1 * (inverseD[i][j] * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + h2f[i][j]) + (1 - omega1) * u[i][j];
				}
			}
		}
		for (i = 0; i < N - 1; i++) {
			for (j = 0; j < N - 1; j++) {
				w[i][j] = v[i][j] - u[i][j];
			}
		}
		error = Frobenius(w);
	}
	end = GetTickCount64();
	cout << "NΪ" << N << "ʱ��SOR��������Ϊ" << m << " CPU��ʱΪ" << end - begin << "ms" << endl;
}