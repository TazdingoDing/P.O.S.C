#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <functions.h>

#define NR_END 1
#define RADIX 2.0
#define FREE_ARG char*
#define SWAP(a,b) do { double t = (a); (a) = (b); (b) = t; } while (0)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) ((a)*(a))

const int length = 3;

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

static void balanc(double **a, int n)
{
	int             i, j, last = 0;
	double          s, r, g, f, c, sqrdx;

	sqrdx = RADIX * RADIX;
	while (last == 0) {
		last = 1;
		for (i = 0; i < n; i++) {
			r = c = 0.0;
			for (j = 0; j < n; j++)
				if (j != i) {
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
				}
			if (c != 0.0 && r != 0.0) {
				g = r / RADIX;
				f = 1.0;
				s = c + r;
				while (c < g) {
					f *= RADIX;
					c *= sqrdx;
				}
				g = r * RADIX;
				while (c > g) {
					f /= RADIX;
					c /= sqrdx;
				}
				if ((c + r) / f < 0.95 * s) {
					last = 0;
					g = 1.0 / f;
					for (j = 0; j < n; j++)
						a[i][j] *= g;
					for (j = 0; j < n; j++)
						a[j][i] *= f;
				}
			}
		}
	}
}

static void elmhes(double **a, int n)
{
	int             i, j, m;
	double          y, x;

	for (m = 1; m < n - 1; m++) {
		x = 0.0;
		i = m;
		for (j = m; j < n; j++) {
			if (fabs(a[j][m - 1]) > fabs(x)) {
				x = a[j][m - 1];
				i = j;
			}
		}
		if (i != m) {
			for (j = m - 1; j < n; j++)
				SWAP(a[i][j], a[m][j]);
			for (j = 0; j < n; j++)
				SWAP(a[j][i], a[j][m]);
		}
		if (x != 0.0) {
			for (i = m + 1; i < n; i++) {
				if ((y = a[i][m - 1]) != 0.0) {
					y /= x;
					a[i][m - 1] = y;
					for (j = m; j < n; j++)
						a[i][j] -= y * a[m][j];
					for (j = 0; j < n; j++)
						a[j][m] += y * a[j][i];
				}
			}
		}
	}
}

static void hqr(double **a, int n, double *wr, double *wi)
{
	int             nn, m, l, k, j, its, i, mmin;
	double          z, y, x, w, v, u, t, s, r, q, p, anorm;

	p = q = r = 0.0;
	anorm = 0.0;
	for (i = 0; i < n; i++)
		for (j = i - 1 > 0 ? i - 1 : 0; j < n; j++)
			anorm += fabs(a[i][j]);
	nn = n - 1;
	t = 0.0;
	while (nn >= 0) {
		its = 0;
		do {
			for (l = nn; l > 0; l--) {
				s = fabs(a[l - 1][l - 1]) + fabs(a[l][l]);
				if (s == 0.0)
					s = anorm;
				if (fabs(a[l][l - 1]) + s == s) {
					a[l][l - 1] = 0.0;
					break;
				}
			}
			x = a[nn][nn];
			if (l == nn) {
				wr[nn] = x + t;
				wi[nn--] = 0.0;
			} else {
				y = a[nn - 1][nn - 1];
				w = a[nn][nn - 1] * a[nn - 1][nn];
				if (l == nn - 1) {
					p = 0.5 * (y - x);
					q = p * p + w;
					z = sqrt(fabs(q));
					x += t;
					if (q >= 0.0) {
						z = p + SIGN(z, p);
						wr[nn - 1] = wr[nn] = x + z;
						if (z != 0.0)
							wr[nn] = x - w / z;
						wi[nn - 1] = wi[nn] = 0.0;
					} else {
						wr[nn - 1] = wr[nn] = x + p;
						wi[nn - 1] = -(wi[nn] = z);
					}
					nn -= 2;
				} else {
					if (its == 30) {
						fprintf(stderr, "[hqr] too many iterations.\n");
						break;
					}
					if (its == 10 || its == 20) {
						t += x;
						for (i = 0; i < nn + 1; i++)
							a[i][i] -= x;
						s = fabs(a[nn][nn - 1]) + fabs(a[nn - 1][nn - 2]);
						y = x = 0.75 * s;
						w = -0.4375 * s * s;
					}
					++its;
					for (m = nn - 2; m >= l; m--) {
						z = a[m][m];
						r = x - z;
						s = y - z;
						p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
						q = a[m + 1][m + 1] - z - r - s;
						r = a[m + 2][m + 1];
						s = fabs(p) + fabs(q) + fabs(r);
						p /= s;
						q /= s;
						r /= s;
						if (m == l)
							break;
						u = fabs(a[m][m - 1]) * (fabs(q) + fabs(r));
						v = fabs(p) * (fabs(a[m - 1][m - 1]) + fabs(z) + fabs(a[m + 1][m + 1]));
						if (u + v == v)
							break;
					}
					for (i = m; i < nn - 1; i++) {
						a[i + 2][i] = 0.0;
						if (i != m)
							a[i + 2][i - 1] = 0.0;
					}
					for (k = m; k < nn; k++) {
						if (k != m) {
							p = a[k][k - 1];
							q = a[k + 1][k - 1];
							r = 0.0;
							if (k + 1 != nn)
								r = a[k + 2][k - 1];
							if ((x = fabs(p) + fabs(q) + fabs(r)) != 0.0) {
								p /= x;
								q /= x;
								r /= x;
							}
						}
						if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0) {
							if (k == m) {
								if (l != m)
									a[k][k - 1] = -a[k][k - 1];
							} else
								a[k][k - 1] = -s * x;
							p += s;
							x = p / s;
							y = q / s;
							z = r / s;
							q /= p;
							r /= p;
							for (j = k; j < nn + 1; j++) {
								p = a[k][j] + q * a[k + 1][j];
								if (k + 1 != nn) {
									p += r * a[k + 2][j];
									a[k + 2][j] -= p * z;
								}
								a[k + 1][j] -= p * y;
								a[k][j] -= p * x;
							}
							mmin = nn < k + 3 ? nn : k + 3;
							for (i = l; i < mmin + 1; i++) {
								p = x * a[i][k] + y * a[i][k + 1];
								if (k != (nn)) {
									p += z * a[i][k + 2];
									a[i][k + 2] -= p * r;
								}
								a[i][k + 1] -= p * q;
								a[i][k] -= p;
							}
						}
					}
				}
			}
		} while (l + 1 < nn);
	}
}

void n_eigen(double *_a, int n, double *wr, double *wi)
{
	int             i;
	double        **a = (double **) calloc(n, sizeof(void *));
	for (i = 0; i < n; ++i)
		a[i] = _a + i * n;
	balanc(a, n);

	elmhes(a, n);
	/*
	int o, p; 

	for (o = 0; o < n; o++) {
		for (p = 0; p < n; p++)
			printf("%13.7e ", a[o][p]);
		printf("\n");
	} 
	*/

	hqr(a, n, wr, wi);

	free(a);
}

static double pythag(double a, double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + SQR(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) printf("allocation failure in vector()");
	return v-nl+NR_END;
}


void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

double** newDoubleArray(int m, int n)
{
	int i, j;
	double **NewMetrix;
	NewMetrix = new double*[m];
	for (i=0;i<m;i++)
	{
		NewMetrix[i] = new double[n];
		for (j=0;j<n;j++)
		{
			NewMetrix[i][j] = 0;
		}
	}
	return NewMetrix;
}

double** dot(double** A, double** B, int Am, int An, int Bm, int Bn)
{
	//if (An != Bm)
	//	return ;
	int i, j, k;
	double tmp;
	double** M = newDoubleArray(Am, Bn);
	for (i=0;i<Am;i++)
	{
		for (j=0;j<Bn;j++)
		{
			tmp=0;
			for(k=0;k<Bm;k++)
			{
				tmp += A[i][k] * B[k][j];
			}
			M[i][j] = tmp;
		}
	}
	
	return M;
}

void svdcmp(double **a, int m,int n,double w[],double **v)
{

	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=vector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_vector(rv1,1,n);
}

double** inv(double **M)
{
	int i,j,k;
	double tmp, *w, **u, **v, **uT, **vT, **dia, **res;
	
	w=vector(1,length);
	u = newDoubleArray(length+1,length+1);
	v = newDoubleArray(length+1,length+1);
	uT = newDoubleArray(length,length);
	vT = newDoubleArray(length,length);
	dia = newDoubleArray(length,length);
	res = newDoubleArray(length,length);
	
	//svd
	for (i=0;i<length;i++)
	{
		for (j=0;j<length;j++)
		{
			u[i+1][j+1] = M[i][j];
		}
	}
	svdcmp(u, length, length, w, v);
	//svd

	for (i=0;i<length;i++)
	{
		for (j=0;j<length;j++)
		{
			uT[i][j] = u[j+1][i+1];
			vT[i][j] = v[i+1][j+1];
			//v return from svd is not vT.
			//no need to transpose here.
		}
		if (w[i+1] > 0)
		{
			dia[i][i] = 1.0/w[i+1];
		}
	}
	
	res = dot(vT, dia, length,length,length,length);
	M = dot(res, uT, length,length,length,length);
	return M;	
}

double norm(double** M, int m, int n)
{
	int i, j;
	double total = 0;
	for (i=0;i<m;i++)
	{
		for (j=0;j<n;j++)
		{
			total += pow(M[i][j], 2);
		}
	}
	return sqrt(total);
}



void MatMinusLambda(double** mat, double lambda)
{
	int i;
	for(i=0;i<length;i++)
	{
		mat[i][i] = mat[i][i] - lambda;
	}
}

void eigTest(double **cov)
{
    int length = 3;
    double u[3], v[3];
    double matrix[3][3];
    double **mat = newDoubleMatrix(length, length);
    for(int i =0;i<3;i++)
    {
        for(int j =0;j<3;j++)
        {
            matrix[i][j] = cov[i][j];
            mat[i][j] = cov[i][j];
        }
    }
    n_eigen(matrix[0], 3, u, v);
    printf("\n");
    for (int i = 0; i < 3; i++)
    {
        printf("eigen value %d: %13.7e + %13.7e J\n",i+1, u[i], v[i]);
    }
    printf("\n");






	for (int i = 0; i < length; i++)
	{
		double **A = newDoubleMatrix(length, length);
		copyMatrix(mat, A, length, length);
		MatMinusLambda(A, u[i]);
		A = inv(A);
		
		double ** b = newDoubleArray(3,1);
		b[0][0] = 1;
		b[1][0] = 2;
		b[2][0] = 3;
		for(int j=0;j<30;j++)
		{
			double ** WB = dot(A, b,3,3,3,1);
			double nor = norm(WB,3,1);
			for(int k=0;k<3;k++)
			{
				b[k][0] = WB[k][0] / nor;
			}
		}
        printf("eigen vector %d:\n", i+1);
		for(int j=0;j<3;j++)
		{
            printf("  %13.7e \n", b[j][0]);
		}
	}
}

void test2()
{
	int             i, j, k;

	static double   u[length], v[length];
	static double   mat_[length][length] = {{3,2,1},{6,5,4}, {9,8,7}};
	static double   ori[length][length] = {{1.0, 6.0, 8.0},{8.0, -15, 17}, {18.0, -3, -24}};

	
	
	double **mat = newDoubleArray(length, length);
	mat[0][0] = mat_[0][0];
	mat[0][1] = mat_[0][1];
	mat[0][2] = mat_[0][2];
	mat[1][0] = mat_[1][0];
	mat[1][1] = mat_[1][1];
	mat[1][2] = mat_[1][2];
	mat[2][0] = mat_[2][0];
	mat[2][1] = mat_[2][1];
	mat[2][2] = mat_[2][2];
	

	n_eigen(mat_[0], length, u, v);
	for (i = 0; i < length; i++)
		printf("%13.7e +J %13.7e\n", u[i], v[i]);
	printf("\n");
	
	
	

	for (i = 0; i < length; i++)
	{
		double **A = newDoubleArray(length, length);
		copyMatrix(mat, A, length, length);
		MatMinusLambda(A, u[i]);
		A = inv(A);
		
		double ** b = newDoubleArray(3,1);
		b[0][0] = 1;
		b[1][0] = 2;
		b[2][0] = 3;
		for(j=0;j<30;j++)
		{
			double ** WB = dot(A, b,3,3,3,1);
			double nor = norm(WB,3,1);
			for(k=0;k<3;k++)
			{
				b[k][0] = WB[k][0] / nor;
			}
		}
		for(j=0;j<3;j++)
		{
			printf("%13.7e \n", b[j][0]);
		}
	}
	
	
	
	/*
	a[0][0] = -6.68465844;
	a[0][1] = 4;
	a[0][2] = 1;
	a[1][0] = 8;
	a[1][1] = -8.68465844;
	a[1][2] = 2;
	a[2][0] = 9;
	a[2][1] = 6;
	a[2][2] = -10.68465844;
	

	
	
	printf("MAT H IS:\n");
	for (i = 0; i < length; i++) {
		for (j = 0; j < length; j++)
			printf("%13.7e ", a[i][j]);
		printf("\n");
	}
	printf("\n");
	
	
	a = inv(a);
	
	double ** b = newDoubleArray(3,1);
	b[0][0] = 1;
	b[1][0] = 2;
	b[2][0] = 3;
	

	
	for(i=0;i<30;i++)
	{
		double ** WB = dot(a, b,3,3,3,1);
		double nor = norm(WB,3,1);
		for(j=0;j<3;j++)
		{
			b[j][0] = WB[j][0] / nor;
		}
	}
	for(j=0;j<3;j++)
	{
		printf("%13.7e \n", b[j][0]);
	}
	*/
	
	
	/*
	printf("res IS:\n");
	for (i = 0; i < length; i++) {
		for (j = 0; j < length; j++)
			printf("%13.7e ", a[i][j]);
		printf("\n");
	}
	printf("\n");
	*/
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/*
	n_eigen(a[0], length, u, v);
	for (i = 0; i < length; i++)
		printf("%13.7e +J %13.7e\n", u[i], v[i]);
	printf("\n");

	//init eigenvector array
	double **vec = new double*[length];
	for (i = 0;i<length;i++)
	{
		vec[i] = new double[length];
		for (j = 0;j<length;j++)
		{
			vec[i][j] = 1.0;
		}
	}
	

	
	
	
	for (i = 0;i<length;i++)
	{
		for (j = 0;j<length;j++)
		{
			printf("%13.7e ", vec[i][j]);
		}
		printf("\n");
	}
	
*/
	//return 0;
}


