#include<iostream>
#include<vector>
using namespace std;

vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b) {
	//计算矩阵左乘向量
	int m = A.size(), n = A[0].size(); vector<double>c(m, 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c[i] += A[i][j] * b[j];
	return c;
}

double power_method_polynomial(vector<double>& a) {
	//用幂法求首一多项式方程的模最大根,输入系数矩阵a={a_n-1...a0}
	int n = a.size(), i, k=0,j; double norm;
	vector<vector<double>>C(n, vector<double>(n, 0)); vector<double>u(n, 0), y(u);
	u[0] = 1;//初值取为e1
	for (i = 0; i < n; i++) {//生成友阵
		C[0][i] = -a[i];
	}
	for (i = 1; i < n; i++) {
		C[i][i - 1] = 1;
	}
	while (k < 100) {//迭代次数设为100
		y = A_multiply_b(C, u);
		norm = fabs(y[0]); j = 0;
		for (i = 1; i < n; i++) {
			if (fabs(y[i])-norm>1e-7) {
				norm = fabs(y[i]);
				j = i;
			}
		}
		for (i = 0; i < n; i++) {
			u[i] = y[i] / y[j];
		}
		k++;
	}
	return y[j];
}

double norm_infinite_vector(vector<double>& b) {
	//计算向量b的无穷范数
	double norm = fabs(b[0]); int n = b.size();
	for (int i = 1; i < n; i++)
		if (fabs(b[i]) > norm)norm = fabs(b[i]);
	return norm;
}
void Householder(vector<double>& x, vector<double>& v, double& beta) {
	//计算x的householder变换，v=x-norm2(x)e1,beta=vt*v
	//这里在beta的存储空间改变beta的值，因此写为&beta
	int n = x.size(), i; double iota = norm_infinite_vector(x), sigma = 0, alpha;
	for (i = 0; i < n; i++) {
		x[i] /= iota;
	}
	for (i = 1; i < n; i++) {
		sigma += x[i] * x[i];
	}
	v = x;//这里虽然只复制2-n的向量，但是v1后重新赋值，可以直接复制向量
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
void transpose(vector<vector<double>>& A, vector<vector<double>>& At) {
	//计算转置矩阵
	int m = A.size(), n = A[0].size();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			At[j][i] = A[i][j];
}
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b) {
	//计算两个向量相乘的矩阵
	int i, j, m = a.size(), n = b.size(); vector <vector<double>>A(m, vector<double>(n));
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			A[i][j] = a[i] * b[j];
	return A;
}
vector<vector<double>>HouseMultiply(vector<vector<double>>& A, vector<double>& v, double beta) {
	//计算一个Householder变换乘矩阵A;
	int m = A.size(), n = A[0].size();
	vector<vector<double>>At(n, vector<double>(m)), B(m, vector<double>(n));
	transpose(A, At);
	vector<double>w = A_multiply_b(At, v);
	for (int i = 0; i < n; i++)w[i] *= beta;
	B = a_mutiply_bt(v, w);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			B[i][j] = A[i][j] - B[i][j];
	return B;
}

void Hessenberg(vector<vector<double>>& A, vector<vector<double>>&Q) {
	//计算A的上Hessenberg分解并储存在A中(householder变换法）,并记录变换的正交矩阵Q
	int n = A.size(), k, i, j; vector<double>v, x, v1(n, 0),v2(v1);
	double beta; vector<vector<double>>A0, A1, H(n, vector<double>(n, 0));
	for (i = 0; i < n; i++)H[i][i] = 1;//H置为单位阵
	for (k = 0; k < n - 2; k++) {
		x.resize(n - k-1); v.resize(n - k-1);
		for (i = k + 1; i < n; i++) {
			x[i-k-1] = A[i][k];//取A要打洞的一列
		}
		Householder(x, v, beta);//计算该列的householder变换
		A0.resize(n - k - 1);
		for (i = 0; i < n-k-1; i++)A0[i].resize(n - k);
		for (i = k + 1; i < n; i++)
			for (j = k; j < n; j++)
				A0[i - k - 1][j - k]=A[i][j];
		A0=HouseMultiply(A0, v, beta);//householder变换左乘矩阵A(A0)
		for (i = k + 1; i < n; i++)
			for (j = k ; j < n; j++)
				A[i][j] = A0[i-k-1][j-k];
		A1.resize(n-k-1);
		for (i = 0; i < n-k-1; i++)A1[i].resize(n);
		for (i = k + 1; i < n; i++)
			for (j = 0; j < n; j++)
				A1[i - k - 1][j] = A[j][i];
		A1 = HouseMultiply(A1, v, beta);//转置乘Householder变换，将右乘变为左乘
		for (i = k + 1; i < n; i++)
			for (j = 0; j < n; j++)
				A[j][i]=A1[i - k - 1][j];
		//以下为记录矩阵Q
		v1 = v2;
		for (i = k + 1; i < n; i++)v1[i] = v[i-k-1];
		H = HouseMultiply(H, v1, beta);//H=Hk……H1
	}
	transpose(H, Q);
	for (i = 2; i < n; i++)//把下三角强制置0
		for (j = 0; j < i - 1; j++)
			A[i][j] = 0;
}

void QR_iteration(vector<vector<double>>& H, vector<vector<double>>& P) {
	//双重步位移的QR迭代一次，变换矩阵记录在P中(这里要求n>=3)
	int n = H.size(), m = n - 1,q,k,i,j,r; double s, t, x, y, z,beta;
	vector<double>v(3, 0), xi(v),v1(n,0),v2(v1); 
	vector<vector<double>>H1(3),P1(n,vector<double>(n,0));
	for (i = 0; i < n; i++)P1[i][i] = 1;
	s = H[m-1][m-1] + H[n-1][n-1];
	t = H[m - 1][m - 1] * H[n - 1][n - 1] - H[m - 1][n - 1]*H[n-1][m-1];
	x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
	y = H[1][0] * (H[0][0] + H[1][1] - s);
	z = H[1][0] * H[2][1];
	for (k = 0; k < n - 2; k++) {
		xi[0] = x; xi[1] = y; xi[2] = z;
		Householder(xi, v, beta);
		if (k > 1)q = k;
		else q = 1;
		for (i = 0; i < 3; i++)H1[i].resize(n - q + 1);
		for (i = 0; i < 3; i++)
			for (j = 0; j < n - q + 1; j++)
				H1[i][j] = H[k + i][j + q - 1];
		H1 = HouseMultiply(H1, v, beta);//左乘一次householder变换
		for (i = 0; i < 3; i++)
			for (j = 0; j < n - q + 1; j++)
				H[k + i][j + q - 1] = H1[i][j];
		if (k + 4 < n)r = k + 4;
		else r = n;
		for (i = 0; i < 3; i++)H1[i].resize(r);
		for (i = 0; i < 3; i++)
			for (j = 0; j < r; j++)
				H1[i][j] = H[j][k + i];
		H1 = HouseMultiply(H1, v, beta);//右乘一次householder变换
		for (i = 0; i < 3; i++)
			for (j = 0; j < r; j++)
				H[j][k + i] = H1[i][j];
		x = H[k + 1][k];
		y = H[k + 2][k];
		if (k < n - 3)z = H[k + 3][k];
		//以下记录矩阵P
		v1 = v2;
		for (i = k; i < k + 3; i++)v1[i] = v[i - k];
		P1 = HouseMultiply(P1, v1, beta);
	}
	xi.resize(2);
	xi[0] = x; xi[1] = y;
	v.resize(2);
	H1.resize(2);
	for (i = 0; i < 2; i++)H1[i].resize(3);
	Householder(xi, v, beta);
	for (i = n - 2; i < n; i++)
		for (j = n - 3; j < n; j++)
			H1[i - n + 2][j - n + 3] = H[i][j];
	H1 = HouseMultiply(H1, v, beta);
	for (i = n - 2; i < n; i++)
		for (j = n - 3; j < n; j++)
			H[i][j]=H1[i - n + 2][j - n + 3];
	for (i = 0; i < 2; i++)H1[i].resize(n);
	for (i = n - 2; i < n; i++)
		for (j = 0; j < n; j++)
			H1[i - n + 2][j] = H[j][i];
	H1 = HouseMultiply(H1, v, beta);
	for (i = n - 2; i < n; i++)
		for (j = 0; j < n; j++)
			H[j][i]=H1[i - n + 2][j];
	v1 = v2;
	for (i = n-2; i < n; i++)v1[i] = v[i - n+2];
	P1 = HouseMultiply(P1, v1, beta);
	transpose(P1, P);
	for (i = 2; i < n; i++)//把下三角强制置0
		for (j = 0; j < i - 1; j++)
			H[i][j] = 0;
}

int Implicit_QR(vector<vector<double>>& A) {
	//隐式QR算法，最终结果存在A中，成功返回0，否则返回1
	int n = A.size(), i, j, m, l, k;
	vector<vector<double>>Q(n, vector<double>(n, 0)), H22, P, Q1, H12, H23;
	Hessenberg(A, Q);
step3:
	for (i = 1; i < n; i++) { //次对角元忽略判定
		if (fabs(A[i][i - 1]) <= (fabs(A[i][i]) + fabs(A[i - 1][i - 1])) * 1e-10)
			A[i][i - 1] = 0;
	}
	//以下找不可约上Hessenberg阵
	for (i = n - 1; i > 0; i--) {
		if (A[i][i - 1] == 0) m = i - 1;//第i行单个实特征值
		else if (i == 1) m = 0;//左上角二阶矩阵
		else {
			if (A[i - 1][i - 2] == 0) {
				i--; m = i - 1;
			}//一个二阶方阵,移动两格
			else {
				m = i; break;
			}
		}
	}//m为H22右下角角标
	if (m == 0) return 0;//已经找到所有特征值
	else {//寻找l为H22左上角标并迭代H22一次
		for (i = m - 1; i > 0; i--) {
			if (A[i][i - 1] == 0) { l = i; break; }
			else l = 0;
		}
		H22.resize(m - l + 1); P.resize(m - l + 1);
		for (i = 0; i < m - l + 1; i++) {
			H22[i].resize(m - l + 1); P[i].resize(m - l + 1);
		}
		for (i = 0; i < m - l + 1; i++)
			for (j = 0; j < m - l + 1; j++)
				H22[i][j] = A[l + i][l + j];
		QR_iteration(H22, P);
	}
	//以下计算新的矩阵H
	Q1.resize(n);
	for (i = 0; i < n; i++)Q1[i].resize(m - l + 1);//计算新的Q
	for (i = 0; i < n; i++)
		for (j = 0; j < m - l + 1; j++) {
			Q1[i][j] = 0;
			for (k = l; k <= m; k++)
				Q1[i][j] += Q[i][k] * P[k - l][j];
		}
	for (i = 0; i < n; i++)
		for (j = l; j <= m; j++)
			Q[i][j] = Q1[i][j-l];
	if (l != 0) {
		H12.resize(l);
		for (i = 0; i < l; i++)H12[i].resize(m - l + 1);//计算新的H12
		for (i = 0; i < l; i++)
			for (j = 0; j < m - l + 1; j++) {
				H12[i][j] = 0;
				for (k = l; k <= m; k++)
					H12[i][j] += A[i][k] * P[k - l][j];
			}
		for (i = 0; i < l; i++)
			for (j = l; j <= m; j++)
				A[i][j] = H12[i][j - l];
	}
	if (m != n - 1) {
		H23.resize(m - l + 1);
		for (i = 0; i < m - l + 1; i++)H23[i].resize(n - m - 1);//计算新的H23
		for (i = 0; i < m - l + 1; i++)
			for (j = 0; j < n - m - 1; j++) {
				H23[i][j] = 0;
				for (k = l; k <= m; k++)
					H23[i][j] += P[k - l][i] * A[k][j];
			}
		for (i = l; i <= m; i++)
			for (j = m + 1; j < n; j++)
				A[i][j] = H23[i - l][j - m - 1];
	}
	for (i = l; i <= m; i++)//复制H22
		for (j = l; j <= m; j++)
			A[i][j] = H22[i - l][j - l];
	Q1.clear(); P.clear(); H22.clear(); H12.clear(); H23.clear();
	goto step3;
	return 1;
}                                                                                                                                                                                                                      

void Eigen_Value(vector<vector<double>>& A) {
	//对一个拟上三角阵，输出其特征值
	int n = A.size(), i; double a, b;
	for (i = 0; i < n; i++) {
		if (i == n - 1)cout << A[i][i] << endl;
		else if(A[i + 1][i] == 0) {
			cout << A[i][i] << endl;
		}
		else {
			a = A[i][i] + A[i + 1][i + 1];
			b = sqrt(fabs(a * a - 4 * (A[i][i] * A[i + 1][i + 1] - A[i + 1][i] * A[i][i + 1])));
			a = a / 2;
			b = b / 2;
			cout << a << "+" << b << "i" << " " << a << "-" << b << "i" << endl;
			i++;
		}

	}
}