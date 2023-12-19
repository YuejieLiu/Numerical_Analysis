#include<iostream>
#include<vector>
using namespace std;

void forward_sub(vector< vector<double> >& L,vector<double>&b) {
	//前代法，解下三角形二维方阵Lx=b,解存储在b中
	int n = L.size()-1;
	for (int j = 0; j <= n - 1;j++) {
		b[j] = b[j] / L[j][j];
		for (int i = j + 1; i <= n; i++) {
			b[i] = b[i] - b[j] * L[i][j];
		}
	}
	b[n] = b[n] / L[n][n];
}

void backward_sub(vector< vector<double> >& U, vector<double>& b) {
	//回代法，解下三角形二维方阵Ux=b,解存储在b中
	int n = U.size() - 1;
	for (int j = n; j >= 1; j--) {
		b[j] = b[j] / U[j][j];
		for (int i = j -1; i>=0; i--) {
			b[i] = b[i] - b[j] * U[i][j];
		}
	}
	b[0] = b[0] / U[0][0];
}

void Gauss(vector< vector<double> >& A) {
	//不选主元的Gauss消去,左下角存储L右上角存储U
	int n = A.size() - 1, k, j, i;
	for (k = 0; k <= n - 1; k++) 
		for (j = k + 1; j <= n; j++) {
			A[j][k] = A[j][k] / A[k][k];
			for (i = k + 1; i <= n; i++)
				A[j][i] = A[j][i] - A[j][k] * A[k][i];
		}
}

void Gauss_Solve(vector<vector<double>>& A,vector<double>&b) {
	//用不选主元的高斯消去法解方程组,解存储在b中
	int k, i, j, n = A.size() - 1;
	Gauss(A);
	vector<vector<double>>c(A);
	for (k = 0; k <= n; k++)
		c[k][k] = 1;
	forward_sub(c, b);
	backward_sub(A, b);
}

void Gauss_Col(vector< vector<double> >& A, vector<int>& u) {
	//列主元高斯消去.A是要分解的矩阵，L存在左下角，U存在左上角。用一个向量u存储置换矩阵P
	int k, p1, p, i, j, n = A.size() - 1; double max, c;
	for (k = 0; k <= n - 1; k++) {
		p1 = k;
		for (p = k + 1; p <= n; p++) {//寻找列主元
			if (fabs(A[p][k])-fabs(A[p1][k])>=1e-7) {
				p1 = p;	
			}
		}
		if (p1 != k) {
			for (i = 0; i <= n; i++) {//交换第k和第p1行
				c = A[k][i];
				A[k][i] = A[p1][i];
				A[p1][i] = c;
			}
		}
		u[k] = p1;//记录置换矩阵Pk
		if (A[k][k]) {
			for (j = k + 1; j <= n; j++) {
				A[j][k] = A[j][k] / A[k][k];
				for (i = k + 1; i <= n; i++)
					A[j][i] = A[j][i] - A[j][k] * A[k][i];
			}
		}	
		else {
			cout << "矩阵奇异";
			break;
		}
	}
}

void Gauss_Col_Solve(vector< vector<double> >& A, vector<double>& b) {
	//用列主元高斯消去解方程组，解存于b中
	int k, n = A.size() - 1; vector<int>u(n, 0); double c; //注意这里置换矩阵只有n-1个
	Gauss_Col(A, u);
	for (k = 0; k <= n-1; k++) {//计算Pb
		c = b[k]; 
		b[k] = b[u[k]];
		b[u[k]] = c;
	}
	vector<vector<double>>B(A);
	for (k = 0; k <= n; k++)
		B[k][k] = 1;
	forward_sub(B, b);//前代法计算Ly=Pb
	backward_sub(A, b);
}

void Square_root(vector< vector<double> >& A) {
	//平方根法。输入方阵A，将L存储在下三角，L的转置存储在上三角
	int k, i, j,n = A.size() - 1;
	for (k = 0; k <= n; k++) {
		if (A[k][k] <= 0) {
			cout << A[k][k] << " ";
			if ((k + 1) % 10 == 0)cout << endl;
		}
		else {
			A[k][k] = sqrt(A[k][k]);
			for (i = k + 1; i <= n; i++)
				A[i][k] = A[i][k] / A[k][k];
			for (j = k + 1; j <= n; j++)
				for (i = j; i <= n; i++)
					A[i][j] = A[i][j] - A[i][k] * A[j][k];
		}
	}
	for (i = 0; i <= n; i++)//将矩阵作对称
		for (j = i + 1; j<=n; j++)
			A[i][j] = A[j][i];
}

void Square_root_Solve(vector<vector<double>>& A, vector<double>& b) {
	//用平方根法计算正定对称系数矩阵方程，结果存于b中
	Square_root(A);
	forward_sub(A, b);
	backward_sub(A, b);
}
void Improved_Square_root(vector< vector<double> >& A) {
	//改进的平方根法A=LDL^t,L储存在左下角，D储存在对角线，L^t储存在右上角
	int j, i,k, n = A.size() - 1; vector<double> v(n);
	for (j = 0; j <= n; j++) {
		for (k = 0; k <= j - 1; k++) {
			v[k] = A[j][k] * A[k][k];
			A[j][j] = A[j][j] - A[j][k] * v[k];
		}
		for (i = j + 1; i <= n; i++) {
			for(k=0;k<=j-1;k++)
				A[i][j] = A[i][j] - A[i][k] * v[k];
			A[i][j] = A[i][j] / A[j][j];
		}
	}
	for (i = 0; i <= n; i++)//将矩阵作对称
		for (j = i + 1; j <= n; j++)
			A[i][j] = A[j][i];
}

void Improved_Square_root_Solve(vector<vector<double>>& A, vector<double>& b) {
	//用改进的平方根法计算正定对称系数矩阵方程，结果存于b中
	int k, n = A.size() - 1;
	Improved_Square_root(A);
	vector<vector<double>>B(A);
	for (k = 0; k <= n; k++)
		B[k][k] = 1;
	forward_sub(B, b);
	for (k = 0; k <= n; k++)
		b[k]=b[k] / A[k][k];
	backward_sub(B, b);
}