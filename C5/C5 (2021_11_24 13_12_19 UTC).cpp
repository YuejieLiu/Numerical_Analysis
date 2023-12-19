#include<iostream>
#include<vector>
#include"C5_func.h" 
#include"C4_func.h"
using namespace std;

int main() {
	int n = 20, i, j, k = 0,m=1; double rho, rho_hat, alpha, beta,h=1.0/n,omega=1.73,error;
	vector<vector<double>>b(n - 1, vector<double>(n - 1, 0)), u(b), r(b), w(b), p(b),inverseD(b),u1(b),v(b);
	//第一题
	//生成向量b
	for (i = 1; i < n - 2; i++) {//去掉外围一圈的元素
		for (j = 1; j < n - 2; j++) {
			b[i][j] = h * h / 4 * sin((i+1) *h* (j+1)*h);
		}
	}
	i = 0;//下边的元素
	for (j = 1; j < n - 2; j++) {
		b[i][j]= h * h / 4 * sin((j+1)*h*h)+0.25*(j+1)*(j+1)*h*h;
	}
	i = n - 2;//上边的元素
	for (j = 1; j < n - 2; j++) {
		b[i][j] = h * h / 4 * sin((i+1)*(j + 1)*h*h) + 0.25 *((i+2)*h*(i+2)*h+(j + 1) * (j + 1)*h*h);
	}
	j = 0;//左边的元素
	for (i = 1; i < n - 2; i++) {
		b[i][j] = h * h / 4 * sin((i+1)*h*h) + 0.25 * (i + 1) * (i + 1)*h*h;
	}
	j = n - 2;//右边的元素
	for (i = 1; i < n - 2; i++) {
		b[i][j] = h * h / 4 * sin((i + 1) * (j + 1)*h*h) + 0.25 * ((j + 2) * (j + 2)*h*h + (i + 1) * (i + 1)*h*h);
	}
	//四个角的元素
	b[0][0]= h * h / 4 * sin(h*h) + 0.5*h*h;
	b[0][n - 2] = h * h / 4 * sin((n - 1)*h*h) + 0.25 * ((n - 1) * (n - 1)*h*h + (h*h + n * n*h*h));
	b[n - 2][0] = b[0][n - 2];
	b[n-2][n-2]= h * h / 4 * sin((n - 1)* (n - 1)*h*h) + 0.5 * ((n - 1) * n*h*h);
	//共轭梯度法
	r = b;//选取的初值为0
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
		//计算w=A*p_k
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
	cout << "共轭梯度法解差分方程的解为（j从小到大,i从小到大）" << endl;
	for (j = 0; j < n - 1; j++) {
		for (i = 0; i < n - 1; i++) {
			cout << u[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	u = u1;
	//思考题
	for (i = 0; i < n - 1; i++) {//生成系数矩阵的对角元D的逆
		for (j = 0; j < n - 1; j++) {
			inverseD[i][j] = 4 / (4 + h * h);
		}
	}
	//SOR迭代
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
	while (error >= 1e-7) {//迭代，终止条件为前后两项二范数之差小于1e-7
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
	cout << "SOR迭代法解差分方程的解为（j从小到大,i从小到大）" << endl;
	for (j = 0; j < n - 1; j++) {
		for (i = 0; i < n - 1; i++) {
			cout << u[i][j] << " ";
		}
		cout << endl;
	}
	cout << "SOR迭代次数为" << m << endl;
	cout << endl;
	//第二题
	n = 20;
	vector<vector<double>>H(n, vector<double>(n, 0));
	vector<double>a(n, 0),x(a);
	for (i = 0; i < n; i++) {//生成希尔伯特矩阵
		for (j = 0; j < n ; j++) {
			H[i][j]=1.0/(i+j+1);
			a[i] += H[i][j] / 3;
		}
	}
	Conjugate_Gradient(H, a, x);
	cout << "第二题方程的解为" << endl;
	for (i = 0; i < n; i++) {
		cout<< x[i] << " ";
		if ((i + 1) % 10 == 0)cout << endl;
	}
	cout << endl;
	//第三题
	vector<vector<double>>B = { {10,1,2,3,4},{1,9,-1,2,-3},{2,-1,7,3,-5},{3,2,3,12,-1},{4,-3,-5,-1,15} };
	vector<double>c = { 12,-27,14,-17,12 }, x1(5, 0);
	Jacobi(B, c, x1, m);
	cout <<"Jacobi迭代法求得矩阵的解为" << endl;
	for (i = 0; i < 5; i++) {
		cout << x1[i] << " ";
	}
	cout << endl;
	cout << "迭代次数为" << m << endl;
	cout << endl;
	GS(B, c, x1, m);
	cout << "GS迭代法求得矩阵的解为" << endl;
	for (i = 0; i < 5; i++) {
		cout << x1[i] << " ";
	}
	cout << endl;
	cout << "迭代次数为" << m << endl;
	cout << endl;
	Conjugate_Gradient(B, c, x1);
	cout << "共轭梯度法求得矩阵的解为" << endl;
	for (i = 0; i < 5; i++) {
		cout << x1[i] << " ";
	}
}