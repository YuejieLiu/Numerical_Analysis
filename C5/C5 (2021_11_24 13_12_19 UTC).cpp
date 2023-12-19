#include<iostream>
#include<vector>
#include"C5_func.h" 
#include"C4_func.h"
using namespace std;

int main() {
	int n = 20, i, j, k = 0,m=1; double rho, rho_hat, alpha, beta,h=1.0/n,omega=1.73,error;
	vector<vector<double>>b(n - 1, vector<double>(n - 1, 0)), u(b), r(b), w(b), p(b),inverseD(b),u1(b),v(b);
	//��һ��
	//��������b
	for (i = 1; i < n - 2; i++) {//ȥ����ΧһȦ��Ԫ��
		for (j = 1; j < n - 2; j++) {
			b[i][j] = h * h / 4 * sin((i+1) *h* (j+1)*h);
		}
	}
	i = 0;//�±ߵ�Ԫ��
	for (j = 1; j < n - 2; j++) {
		b[i][j]= h * h / 4 * sin((j+1)*h*h)+0.25*(j+1)*(j+1)*h*h;
	}
	i = n - 2;//�ϱߵ�Ԫ��
	for (j = 1; j < n - 2; j++) {
		b[i][j] = h * h / 4 * sin((i+1)*(j + 1)*h*h) + 0.25 *((i+2)*h*(i+2)*h+(j + 1) * (j + 1)*h*h);
	}
	j = 0;//��ߵ�Ԫ��
	for (i = 1; i < n - 2; i++) {
		b[i][j] = h * h / 4 * sin((i+1)*h*h) + 0.25 * (i + 1) * (i + 1)*h*h;
	}
	j = n - 2;//�ұߵ�Ԫ��
	for (i = 1; i < n - 2; i++) {
		b[i][j] = h * h / 4 * sin((i + 1) * (j + 1)*h*h) + 0.25 * ((j + 2) * (j + 2)*h*h + (i + 1) * (i + 1)*h*h);
	}
	//�ĸ��ǵ�Ԫ��
	b[0][0]= h * h / 4 * sin(h*h) + 0.5*h*h;
	b[0][n - 2] = h * h / 4 * sin((n - 1)*h*h) + 0.25 * ((n - 1) * (n - 1)*h*h + (h*h + n * n*h*h));
	b[n - 2][0] = b[0][n - 2];
	b[n-2][n-2]= h * h / 4 * sin((n - 1)* (n - 1)*h*h) + 0.5 * ((n - 1) * n*h*h);
	//�����ݶȷ�
	r = b;//ѡȡ�ĳ�ֵΪ0
	rho = b_multiply_c_matrix(r, r);
	while (sqrt(rho) > 1e-7 * Frobenius(b) && k < n) {
		k++;
		if (k == 1)p = r;
		else {
			beta = rho / rho_hat;
			for (i = 0; i < n-1; i++) {
				for (j = 0; j < n - 1; j++) {
					p[i][j] = r[i][j] + beta * p[i][j];
				}			
			}
		}
		//����w=A*p_k
		for (j = 0; j < n - 1; j++) {
			for (i = 0; i < n - 1; i++) {
				if (i == 0) {
					if (j == 0)w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i + 1][j]+ p[i][j + 1]);
					else if (j == n - 2)w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i + 1][j]+ p[i][j - 1]);
					else w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i + 1][j]+ p[i][j + 1] + p[i][j - 1]);
				}
				else if (i == n - 2) {
					if (j == 0)w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i - 1][j] + p[i][j + 1]);
					else if (j == n - 2)w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25* ( p[i - 1][j]+ p[i][j - 1]);
					else w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i - 1][j] + p[i][j + 1] + p[i][j - 1]);
				}
				else {
					if (j == 0)w[i][j] = (1 + h * h / 4) * p[i][j] -0.25* (p[i + 1][j] + p[i - 1][j] + p[i][j + 1]);
					else if (j == n - 2)w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i + 1][j] + p[i - 1][j]+ p[i][j - 1]);
					else w[i][j] = (1 + h * h / 4) * p[i][j] - 0.25 * (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1]);
				}
			}
		}
		alpha = rho / b_multiply_c_matrix(p, w);
		for (i = 0; i < n-1; i++) {
			for (j = 0; j < n - 1; j++) {
				u[i][j] = u[i][j] + alpha * p[i][j];
			}
		}
		for (i = 0; i < n-1; i++) {
			for (j = 0; j < n - 1; j++) {
				r[i][j] = r[i][j] - alpha * w[i][j];
			}
		}
		rho_hat = rho;
		rho = b_multiply_c_matrix(r, r);
	}
	cout << "�����ݶȷ����ַ��̵Ľ�Ϊ��j��С����,i��С����" << endl;
	for (j = 0; j < n - 1; j++) {
		for (i = 0; i < n - 1; i++) {
			cout << u[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	u = u1;
	//˼����
	for (i = 0; i < n - 1; i++) {//����ϵ������ĶԽ�ԪD����
		for (j = 0; j < n - 1; j++) {
			inverseD[i][j] = 4 / (4 + h * h);
		}
	}
	//SOR����
	for (j = 0; j < n - 1; j++) {
		for (i = 0; i < n - 1; i++) {
			if (i == 0) {
				if (j == 0)u[i][j] = omega * inverseD[i][j] * ( 0.25*(u[i + 1][j] + u[i][j + 1]) +b[i][j]) + (1 - omega) * u[i][j];
				else if (j == n - 2)u[i][j] = omega * inverseD[i][j] * (0.25*(u[i + 1][j] + u[i][j - 1])+b[i][j]) + (1 - omega) * u[i][j];
				else u[i][j] = omega * inverseD[i][j] * (0.25*(u[i][j - 1] + u[i + 1][j] + u[i][j + 1])+b[i][j]) + (1 - omega) * u[i][j];
			}
			else if (i == n - 2) {
				if (j == 0)u[i][j] = omega*inverseD[i][j] * (0.25*(u[i - 1][j] + u[i][j + 1]) +b[i][j]) + (1 - omega) * u[i][j];
				else if (j == n - 2)u[i][j] = omega* inverseD[i][j] * (0.25*(u[i - 1][j] + u[i][j - 1]) + b[i][j]) + (1 - omega) * u[i][j];
				else u[i][j] = omega*inverseD[i][j] * (0.25*(u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
			}
			else {
				if (j == 0)u[i][j] = omega*inverseD[i][j] * (0.25*(u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
				else if (j == n - 2)u[i][j] = omega* inverseD[i][j] *(0.25*( u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + b[i][j]) + (1 - omega) * u[i][j];
				else u[i][j] = omega*inverseD[i][j] * (0.25*(u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
			}
		}
	}
	for (i = 0; i < n - 1; i++) {
		for (j = 0; j < n - 1; j++) {
			w[i][j] = v[i][j] - u[i][j];
		}
	}
	error = Frobenius(w);
	while (error >= 1e-7) {//��������ֹ����Ϊǰ�����������֮��С��1e-7
		v = u; m++;
		for (j = 0; j < n - 1; j++) {
			for (i = 0; i < n - 1; i++) {
				if (i == 0) {
					if (j == 0)u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i + 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
					else if (j == n - 2)u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i + 1][j] + u[i][j - 1]) + b[i][j]) + (1 - omega) * u[i][j];
					else u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
				}
				else if (i == n - 2) {
					if (j == 0)u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i - 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
					else if (j == n - 2)u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i - 1][j] + u[i][j - 1]) + b[i][j]) + (1 - omega) * u[i][j];
					else u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i][j - 1] + u[i - 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
				}
				else {
					if (j == 0)u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
					else if (j == n - 2)u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1]) + b[i][j]) + (1 - omega) * u[i][j];
					else u[i][j] = omega * inverseD[i][j] * (0.25 * (u[i - 1][j] + u[i][j - 1] + u[i + 1][j] + u[i][j + 1]) + b[i][j]) + (1 - omega) * u[i][j];
				}
			}
		}
		for (i = 0; i < n - 1; i++) {
			for (j = 0; j < n - 1; j++) {
				w[i][j] = v[i][j] - u[i][j];
			}
		}
		error = Frobenius(w);
	}
	cout << "SOR���������ַ��̵Ľ�Ϊ��j��С����,i��С����" << endl;
	for (j = 0; j < n - 1; j++) {
		for (i = 0; i < n - 1; i++) {
			cout << u[i][j] << " ";
		}
		cout << endl;
	}
	cout << "SOR��������Ϊ" << m << endl;
	cout << endl;
	//�ڶ���
	n = 20;
	vector<vector<double>>H(n, vector<double>(n, 0));
	vector<double>a(n, 0),x(a);
	for (i = 0; i < n; i++) {//����ϣ�����ؾ���
		for (j = 0; j < n ; j++) {
			H[i][j]=1.0/(i+j+1);
			a[i] += H[i][j] / 3;
		}
	}
	Conjugate_Gradient(H, a, x);
	cout << "�ڶ��ⷽ�̵Ľ�Ϊ" << endl;
	for (i = 0; i < n; i++) {
		cout<< x[i] << " ";
		if ((i + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	//������
	vector<vector<double>>B = { {10,1,2,3,4},{1,9,-1,2,-3},{2,-1,7,3,-5},{3,2,3,12,-1},{4,-3,-5,-1,15} };
	vector<double>c = { 12,-27,14,-17,12 }, x1(5, 0);
	Jacobi(B, c, x1, m);
	cout <<"Jacobi��������þ���Ľ�Ϊ" << endl;
	for (i = 0; i < 5; i++) {
		cout << x1[i] << " ";
	}
	cout << endl;
	cout << "��������Ϊ" << m << endl;
	cout << endl;
	GS(B, c, x1, m);
	cout << "GS��������þ���Ľ�Ϊ" << endl;
	for (i = 0; i < 5; i++) {
		cout << x1[i] << " ";
	}
	cout << endl;
	cout << "��������Ϊ" << m << endl;
	cout << endl;
	Conjugate_Gradient(B, c, x1);
	cout << "�����ݶȷ���þ���Ľ�Ϊ" << endl;
	for (i = 0; i < 5; i++) {
		cout << x1[i] << " ";
	}
}