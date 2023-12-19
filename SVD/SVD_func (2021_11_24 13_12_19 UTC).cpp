#include<iostream>
#include<vector>
using namespace std;

double norm_infinite_vector(vector<double>& b) {
	//计算向量b的无穷范数
	double norm = fabs(b[0]); int n = b.size();
	for (int i = 1; i < n; i++)
		if (fabs(b[i]) > norm)norm = fabs(b[i]);
	return norm;
}
void print(vector<vector<double>>& A) {
	//打印矩阵的所有元素
	int m = A.size(), n = A[0].size(), i, j;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
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
vector<double>A_multiply_b(vector<vector<double>>& A, vector<double>& b) {
	//计算矩阵左乘向量
	int m = A.size(), n = A[0].size(); vector<double>c(m, 0);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			c[i] += A[i][j] * b[j];
	return c;
}
vector<vector<double>>a_mutiply_bt(vector<double>& a, vector<double>& b) {
	//计算两个向量相乘的矩阵
	int i, j, m = a.size(), n = b.size(); vector <vector<double>>A(m, vector<double>(n));
	for (i = 0; i < m; i++)
		for (j = 0; j < n; j++)
			A[i][j] = a[i] * b[j];
	return A;
}
void bi_diagonal_house(vector<vector<double>>& A,vector<double>&b,vector<double>&c) {
	//用householder变换法将矩阵阵A化为二对角矩阵,A的行数m大于等于列数n
	//计算后A中二对角元素下方是左乘Householder的v，右上角是右乘householder变换的v
	//b中是左乘的beta,有至多n个元素，c中是右乘的beta，有n-2个元素
	int m = A.size(),n=A[0].size(), i, j,k; vector<double>x, v,u; double beta;
	vector<vector<double>>B,C;
	for (j = 0; j < n; j++) {
		if (j < m - 1) {//计算左乘householder
			x.resize(m - j); v.resize(m - j); u.resize(n - j);
			for (i = j; i < m; i++) {
				x[i - j] = A[i][j];
			}
			Householder(x, v, beta);
			B.resize(n - j); C.resize(n - j);
			for (i = j; i < n; i++) {
				B[i - j].resize(m - j);
				C[i - j].resize(m - j);
			}
			for (i = j; i < n; i++)
				for (k = j; k < m; k++) {
					B[i - j][k - j] = A[k][i];
				}
			u = A_multiply_b(B, v);
			for (i = j; i < n; i++)u[i - j] = beta * u[i - j];
			C = a_mutiply_bt(u, v);
			for (i = j; i < n; i++)
				for (k = j; k < m; k++) {
					B[i - j][k - j] -= C[i - j][k - j];
				}
			for (i = j; i < n; i++)//householder对矩阵的作用计算完毕
				for (k = j; k < m; k++) {
					A[k][i] = B[i - j][k - j];
				}
			for (i = j + 1; i < m; i++) {//记录变换的v和beta
				A[i][j] = v[i - j];
			}
			b[j] = beta;
			x.clear(); v.clear(); u.clear(); B.clear(); C.clear();
		}
		if (j < n - 2) {//计算右乘Householder变换
			x.resize(n - j - 1); v.resize(n - j - 1); u.resize(m - j);
			for (i = j+1; i < n; i++) {
				x[i - j-1] = A[j][i];
			}
			Householder(x, v, beta);
			B.resize(m - j); C.resize(m - j);
			for (i = j; i < m; i++) {
				B[i - j].resize(n - j-1);
				C[i - j].resize(n - j-1);
			}
			for (i = j; i < m; i++)
				for (k = j+1; k < n; k++) {
					B[i - j][k - j-1] = A[i][k];
				}
			u = A_multiply_b(B, v);
			for (i = j; i < m; i++)u[i - j] = beta * u[i - j];
			C = a_mutiply_bt(u, v);
			for (i = j; i < m; i++)
				for (k = j+1; k < n; k++) {
					B[i - j][k - j-1] -= C[i - j][k - j-1];
				}
			for (i = j; i < m; i++)//householder对矩阵的作用计算完毕
				for (k = j+1; k < n; k++) {
					A[i][k] = B[i - j][k - j-1];
				}
			for (i = j + 2; i < n; i++)A[j][i] = v[i - j - 1];
			c[j] = beta;
			x.clear(); v.clear(); u.clear(); B.clear(); C.clear();
		}
	}
	//注意，以上我们作用的householder变换是>=2阶的，因此无法保证最后一个元素也非负数，而SVD分解要求全正数
}
void transpose(vector<vector<double>>& A, vector<vector<double>>& At) {
	//计算转置矩阵
	int m = A.size(), n = A[0].size();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			At[j][i] = A[i][j];
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
void accumulateUV(vector<vector<double>>& A,vector<double>& b, vector<double>& c, vector<vector<double>>&U, vector<vector<double>>&V) {
	//对已经二对角化的矩阵计算U和V
	//b中是左乘的beta,有至多n个元素，c中是右乘的beta，有n-2个元素
	int m = A.size(), n = A[0].size(), i, j;
	vector<double>v(m, 0),v1(v), u(n, 0),u1(u);
	vector<vector<double>>U1(m, vector<double>(m, 0)), V1(n, vector<double>(n, 0));
	for (i = 0; i < m; i++)U1[i][i] = 1.0;
	for (i = 0; i < n; i++)V1[i][i] = 1.0;
	for (j = 0; j < n; j++) {//计算左乘的矩阵U,注意U对称
		if (j < m - 1) {
			v[j] = 1;
			for (i = j + 1; i < m; i++)v[i] = A[i][j];
			U1=HouseMultiply(U1, v, b[j]);
			v = v1;
		}
	}
	for (j = 0; j < n; j++) {//计算右乘的矩阵V
		if (j < n - 2) {
			u[j+1] = 1;
			for (i = j + 2; i < n; i++)u[i] = A[j][i];
			V1=HouseMultiply(V1, u, c[j]);
			u = u1;
		}
	}
	transpose(V1, V);
	transpose(U1, U);
}

void Givens(double a, double b, double& c, double& s) {
	//计算givens变换(第二个分量为0）
	double tau;
	if (b == 0) {
		c = 1; s = 0;
	}
	else {
		if (fabs(b) > fabs(a)) {
			tau = a / b;
			s = 1 / sqrt(1 + tau * tau);
			c = s * tau;
		}
		else {
			tau = b / a;
			c=1/ sqrt(1 + tau * tau);
			s = c * tau;
		}
	}
}
void Wilkinson_SVD(vector<vector<double>>& A,vector<double>&Qc,vector<double>&Qs,vector<double>&Pc,vector<double>&Ps) {
	//对已经二对角化的带Wilkinson位移的SVD迭代（不可约）(n>=3)，
	//其中Qc Qs Ps Pc是n-1维向量，记录givens变换
	int n = A.size(), i, j; double alpha, delta, beta, miu, y, z, c, s,a,b;
	vector<double>d(n, 0), gamma(n - 1, 0);
	for (i = 0; i < n - 1; i++) {//取出二对角的元素
		d[i] = A[i][i];
		gamma[i] = A[i][i + 1];
	}
	d[n - 1] = A[n - 1][n - 1];//计算位移
	alpha = d[n - 1] * d[n - 1] + gamma[n - 2] * gamma[n - 2];
	delta = (d[n - 2] * d[n - 2] + gamma[n - 3] * gamma[n - 3] - alpha)/2;
	beta = d[n - 2] * gamma[n - 2];
	if (delta > 0)miu = alpha - beta * beta / (delta + sqrt(delta * delta + beta * beta));
	else miu = alpha - beta * beta / (delta - sqrt(delta * delta + beta * beta));
	y = d[0] * d[0] - miu;
	z = d[0] * gamma[0];
	for (i = 0; i < n - 1; i++) {//givens变换
		Givens(y, z, c, s);//计算右乘的givens变换,这里取得是转置
		Qc[i] = c; Qs[i] = -s;
		if (i != 0)gamma[i - 1] = c * y + s * z;
		y = c * d[i] + s * gamma[i];
		gamma[i] = -s * d[i] + c * gamma[i];
		z = s * d[i + 1];
		d[i + 1] = c * d[i + 1];
		Givens(y, z, c, s);//计算左乘的givens变换,这里取得是转置
		Pc[i] = c; Ps[i] = -s;
		d[i] = c * y + s * z;
		if (i < n - 2) {
			y = c * gamma[i] + s * d[i + 1];
			z = s * gamma[i + 1];
			d[i + 1] = -s * gamma[i] + c * d[i + 1];
			gamma[i + 1] = c * gamma[i + 1];
		}
		else {
			a = gamma[i]; b = d[i + 1];
			gamma[i] = c * a + s * b;
			d[i + 1] = -s * a + c * b;
		}
	}
	for (i = 0; i < n - 1; i++) {//返回二对角的元素
		A[i][i]=d[i];
		A[i][i + 1]= gamma[i];
	}
	A[n - 1][n - 1]=d[n - 1];
}
void Wilkinson2(vector<vector<double>>& A, vector<double>& Qc, vector<double>& Qs, vector<double>& Pc, vector<double>& Ps) {
	//对已经二对角化的带Wilkinson位移的SVD迭代（不可约）(n=2)，
	//其中Qc Qs Ps Pc是1维向量，记录givens变换
	int i, j; double alpha, delta, beta, miu, y, z, c, s,gamma,a,b;
	vector<double>d(2, 0);
	//取出二对角的元素
	d[0] = A[0][0];
	gamma = A[0][1];
	d[1] = A[1][1];
	//计算位移
	alpha = d[1] * d[1] + gamma * gamma;
	delta = (d[0] * d[0] - alpha) / 2;
	beta = d[0] * gamma;
	if (delta > 0)miu = alpha - beta * beta / (delta + sqrt(delta * delta + beta * beta));
	else miu = alpha - beta * beta / (delta - sqrt(delta * delta + beta * beta));
	y = d[0] * d[0] - miu;
	z = d[0] * gamma;
	//givens变换
	Givens(y, z, c, s);//计算右乘的givens变换
	s = -s;//这里取得是转置
	Qc[0] = c; Qs[0] = s;
	y = c * d[0] - s * gamma;
	gamma = s * d[0] + c * gamma;
	z = -s * d[1];
	d[1] = c * d[1];
	Givens(y, z, c, s);//计算左乘的givens变换
	s = -s;//这里取得是转置
	Pc[0] = c; Ps[0] = s;
	d[0] = c * y - s * z;
	a = gamma; b = d[1];
	gamma = c * a- s * b;
	d[1] = s * a + c * b;
	//返回二对角的元素
	A[0][0]=d[0];
	A[0][1]=gamma;
	A[1][1]=d[1];
}
double matrix_infinite_norm(vector<vector<double>>& A) {
	//计算矩阵的无穷范数
	int m = A.size(), n = A[0].size(), i, j; double max=0, norm=0;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++)max += fabs(A[i][j]);
		if (max > norm)norm = max;
		max = 0;
	}
	return norm;
}
vector<vector<double>>multiply(vector<vector<double>>& A,vector<vector<double>>& B) {
	//计算矩阵AB
	int m = A.size(), n = A[0].size(), l = B[0].size(), i, j, k;
	vector<vector<double>>C(m, vector<double>(l, 0));
	for (i = 0; i < m; i++) {
		for (j = 0; j < l; j++) {
			for (k = 0; k < n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return C;
}
vector<vector<double>>matrix_minus(vector<vector<double>>& A, vector<vector<double>>& B) {
	//计算矩阵相减A-B
	int m = A.size(), n = A[0].size(), i, j;
	vector<vector<double>>C(m, vector<double>(n, 0));
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			C[i][j] = A[i][j] - B[i][j];
		}
	}
	return C;
}
double max_mod(vector<vector<double>>& A) {
	//计算矩阵的最大模元素
	int m = A.size(), n = A[0].size(), i, j,p,q; double mod=0;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			if (fabs(A[i][j]) > mod) {
				mod = fabs(A[i][j]);
				p = i; q = j;
			}
		}
	}
	return A[p][q];
}
void Givens2(double a, double b, double& c, double& s) {
	//计算givens变换(第一个分量为0）
	double tau;
	if (a == 0) {
		c = 1; s = 0;
	}
	else {
		if (fabs(a) > fabs(b)) {
			tau = b / a; 
			s = 1 / sqrt(1 + tau * tau);
			c = -s * tau;
		}
		else {
			tau = a / b;
			c = -1 / sqrt(1 + tau * tau);
			s = -c * tau;
		}
	}
}
int SVD(vector<vector<double>>& A,vector<vector<double>>&U,vector<vector<double>>&V,vector<vector<double>>&B) {
	//求m*n阶矩阵A的奇异值分解（m>=n)(B是0矩阵）
	int m = A.size(), n = A[0].size(), i,p,q,k,j; double norm,c2,s2,a2,b2;
	vector<double>b(n, 0), c(n - 2, 0),Qc,Qs,Pc,Ps;
	vector<vector<double>>Ut(U),B22;
	bi_diagonal_house(A, b, c);
	accumulateUV(A, b, c, U, V);
	for (i = 0; i < n-1; i++) {//B是A中存储的二对角矩阵
		B[i][i] = A[i][i];
		B[i][i + 1] = A[i][i + 1];
	}
	B[n - 1][n - 1] = A[n - 1][n - 1];
	//收敛性判定
step3:	
	norm = matrix_infinite_norm(B);
	for (i = 0; i < n - 1; i++) {
		if (fabs(B[i][i + 1]) <= (fabs(B[i][i]) + fabs(B[i + 1][i + 1])) * 1e-12)B[i][i + 1] = 0;
	}
	for (i = 0; i < n - 1; i++) {
		if (fabs(B[i][i]) <= norm * 1e-12)B[i][i] = 0;
	}
	for (i = n - 1; i > 0; i--) {
		if (B[i - 1][i] == 0) {
			q = 0;
		}
		else {
			q = i; break;
		}//找到q为B22的最后一行坐标
	}
	if (q == 0) {//已得到分解Ut*A*V=(B0),现在通过行列变换和正负号变换得到奇异值分解中非负递减排列的要求
		for (i = 0; i < n; i++) {//把对角元变成正数，变号矩阵累积的V上
			if (B[i][i] < 0) {
				B[i][i] = fabs(B[i][i]);
				for (j = 0; j < n; j++)V[j][i] *= -1;
			}
		}
		for (i = 0; i < n; i++) {//排序
			double max = B[i][i]; int label = i;
			for (j = i + 1; j < n; j++) {//找到最大元
				if (B[j][j] > max) {
					max = B[j][j]; label = j;
				}
			}
			double mid = B[i][i]; B[i][i] = B[label][label]; B[label][label] = mid;//互换位置
			for (j = 0; j < m; j++) {//互换U的两列（Ut的两行）
				mid = U[j][i]; U[j][i] = U[j][label]; U[j][label] = mid;
			}
			for (j = 0; j < n; j++) {//互换V的两列
				mid = V[j][i]; V[j][i] = V[j][label]; V[j][label] = mid;
			}
		}
		cout << "正交矩阵U为" << endl;
		print(U);
		cout << "正交矩阵V为" << endl;
		print(V);
		cout << "奇异值矩阵B为" << endl;
		print(B);
		return 0;
	}
	for (i = q-1; i >=0; i--) {//找到p为B22的第一行坐标
		if (i == 0) {
			p = 0; break;
		}
		if (B[i - 1][i] == 0) {
			p = i; break;
		}
	}
	//SVD迭代
	k = q;
	for (i = p; i < q; i++) {
		if (B[i][i] == 0) { k = i; break; }
	}//判断B22在B中对角元为0的位置k
	if (k != q) {//有对角元为0,打洞，累加，gotostep3
		transpose(U, Ut);
		for (i = k; i < q; i++) {
			a2 = B[k][i + 1]; b2 = B[i + 1][i + 1];
			Givens2(a2,b2, c2, s2);
			B[k][i + 1] = c2 * a2 + s2 * b2;
			B[i + 1][i + 1] = -s2 * a2 + c2 * b2;
			if (i != q - 1) {
				a2 = B[k][i + 2]; b2 = B[i + 1][i + 2];
				B[k][i+2]= c2 * a2 + s2 * b2;
				B[i+1][i+2]= -s2 * a2 + c2 * b2;
			}
			//把givens变换累积到Ut上
			for (j = 0; j < m; j++) {
				a2 = Ut[k][j]; b2 = Ut[i+1][j];
				Ut[k][j] = c2 * a2 + s2 * b2;
				Ut[i + 1][j] = -s2 * a2 + c2 * b2;
			}
		}
		transpose(Ut, U);
		goto step3;
	}
	else{
		B22.resize(q - p + 1);
		for (i = 0; i < q - p + 1; i++) {
			B22[i].resize(q - p + 1);
		}
		for (i = 0; i < q - p; i++) {//单独复制从B中复制B22
			B22[i][i] = B[i + p][i + p];
			B22[i][i + 1] = B[i + p][i + p + 1];
		}
		B22[q - p][q - p] = B[q][q];
		Qc.resize(q - p); Qs.resize(q - p); Pc.resize(q - p); Ps.resize(q - p);
		if (q - p + 1 >= 3) {//对>=3的进行一次SVD迭代
			Wilkinson_SVD(B22, Qc, Qs, Pc, Ps);
		}
		else {//对2阶迭代单独写函数讨论
			Wilkinson2(B22, Qc, Qs, Pc, Ps);
		}
		for (i = 0; i < q - p; i++) {//把B22并清除
			B[i + p][i + p]=B22[i][i];
			B[i + p][i + p + 1]=B22[i][i + 1];
		}
		B[q][q]=B22[q - p][q - p];
		for (i = 0; i < q - p; i++) {//累计正交矩阵U和V
			for (j = 0; j < m; j++) {//累计U
				a2 = U[j][p + i]; b2 = U[j][p + i + 1];
				U[j][p + i] = Pc[i] * a2 - Ps[i] * b2;
				U[j][p + i + 1] = Ps[i] * a2 + Pc[i] * b2;
			}
			for (j = 0; j < n; j++) {//累计V
				a2 = V[j][p + i]; b2 = V[j][p + i + 1];
				V[j][p + i] = Qc[i] * a2 - Qs[i] * b2;
				V[j][p + i + 1] = Qs[i] * a2 + Qc[i] * b2;
			}
		}
		Qc.clear(); Qs.clear(); Pc.clear(); Ps.clear();
		goto step3;
	}
	return 1;
}