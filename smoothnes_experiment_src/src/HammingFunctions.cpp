/*
 * hamming_functions.cpp
 *
 *  Created on: Jun 6, 2018
 */
#include "HammingFunctions.h"
#include <iostream>
#include <cmath>
#include "Lap.h"
#include "QAP.h"
#include <cstdlib>
#include <float.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include "Tools.h"
#include <assert.h>

#define MAXIT 200

//offet between first perm and first index
#define FIRST_ITEM 1

#define ALLOW_0_DISTANCE_SAMPLING false

#define USE_FAST_FITNESS_IF_DISTANCE_LOWER_OR_EQUAL_TO 2

// random number generation

//void HammingFunctions::seed(void) {
//	int fd, buf;
//	if ((fd = open("/dev/urandom", O_RDONLY)) < 0) {
//
//		perror("/dev/urandom");
//	}
//	read(fd, &buf, sizeof(buf));
//	close(fd);
//	srand(buf);
//}

void HammingFunctions::generate_random_permutation(int *sigma)
{
	//implements knuth shuffle (Fisher–Yates shuffle)
	//other option is to code Feller coupling
	//seed();
	for (int i = 0; i < n_; i++)
	{
		sigma[i] = i + 1;
	}
	for (int i = 0; i < n_ - 1; i++)
	{
		int pos = rand() % (n_ - i) + i;
		//int pos = (int) (unif_rand() * (len-i) + i);
		int aux = sigma[i];
		sigma[i] = sigma[pos];
		sigma[pos] = aux;
	}
}

void HammingFunctions::generate_random_sample(int m, int **sample)
{
	for (int i = 0; i < m; i++)
	{
		generate_random_permutation(sample[i]);
	}
}

long double HammingFunctions::count_permus_with_at_least_k_unfixed_points(int n, int k)
{
	//if ( facts_ == NULL ) {
	//   init_factorials(n, facts_);}
	//else if (facts_n_ < n ) {cout<<"Check n in Generic::count_permus_no_fixed_points. ";exit(0);}
	long double sum = 0, aux = 0;
	int multi = -1;
	for (int i = 1; i <= k; i++)
	{
		//num = num + (-1)^j * factorial(n-l) * factorial(l) / (factorial(j) * factorial(l-j) );
		aux = (long double)multi * (long double)facts_[k] * facts_[n - i] / (facts_[i] * facts_[k - i]);
		sum += aux;
		multi *= -1;
	}
	return facts_[n] + sum;
}

void HammingFunctions::random_derangement(int n, int *sigma)
{
	

	if (n == 2)
	{
		sigma[0] = 2;
		sigma[1] = 1;
	}
	else if ((n - 1) * deran_num_[n - 1] / deran_num_[n] > (double)rand() / RAND_MAX || n == 3)
	{
		//}else if ( (n-1)*deran_num_[n-1] / deran_num_[n] > unif_rand() ){
		random_derangement(n - 1, sigma);
		int ran = rand() % (n - 1);
		//int ran = (int) (unif_rand() * (n - 1 ));
		sigma[n - 1] = sigma[ran];
		sigma[ran] = n;
	}
	else
	{
		int *deran = new int[n - 2], *conv = new int[n - 1];
		random_derangement(n - 2, deran);
		int ran = rand() % (n - 1);
		//int ran =(int)( unif_rand() * (n - 1 ));
		int j = 0;
		for (int i = 0; i < n - 1; i++)
			if (i != ran)
				conv[j++] = i + 1;
		j = 0;
		for (int i = 0; i < n - 1; i++)
			if (i != ran)
				sigma[i] = conv[deran[j++] - 1];
		sigma[ran] = n;
		sigma[n - 1] = ran + 1;
		delete[] deran;
		delete[] conv;
	}
}

void HammingFunctions::generate_permu_from_list(int *ran, int dist, int *sigma)
{
	//the first d items in ran will be deranged. for the rest, sigma[i]=i
	if (dist == 0)
	{
		for (int i = 0; i < n_; i++)
			sigma[i] = i + 1;
		return;
	}
	int *deran = new int[n_];
	if (dist > 1)
		random_derangement(dist, deran);
	for (int i = 0; i < dist; i++)
		sigma[ran[i] - 1] = ran[deran[i] - 1];
	for (int i = dist; i < n_; i++)
		sigma[ran[i] - 1] = ran[i];
	delete[] deran;
}

void HammingFunctions::dist_decomp_vector2perm(int *vec, int *sigma)
{
	int last = n_ - 1;
	int first = 0;
	for (int i = 0; i < n_; i++)
	{
		if (vec[i] == 0)
			ran_[last--] = i + 1;
		else
			ran_[first++] = i + 1;
	}
	//if ( first == 1 )bool traza = true;
	generate_permu_from_list(ran_, first, sigma);
	//Generic gen;cout<< "h ";gen.print_int_vector(h, n_); cout<< "p ";gen.print_int_vector(permu, n_);
}

long double HammingFunctions::compute_marginal_iterative(int *h, double *theta,
														 int marginal_order)
{
	// must be initialized :
	//esp_ini: elementary_symmetric_polynomial( theta, n_ ,..., esp_ini_ );
	//t_sampling_[i]=exp(theta[i]-1
	int a = 0; //a : |A|; b = |B| is global
	long double result = 0;
	int current_var = marginal_order - 1;
	int num_vars = n_ - marginal_order;
	long double res = 0;

	if (marginal_order == 1)
	{ //the first time it is called

		theta_acum_not_in_A = 0;
		b_ = 0;
		for (int i = 0; i < n_; i++)
		{

			theta_acum_not_in_A += (long double)theta[i];

			esp_red_[i] = esp_ini_[i];
		}
		esp_red_[n_] = esp_ini_[n_];
	}
	if (current_var > 0)
	{
		if (h[current_var - 1] == 0)
			theta_acum_not_in_A -= (long double)theta[current_var - 1]; //the set of fixed points is A. If the last position was sampled as Unfix , update set
		else
			b_++; //otherwise, (current_var-1) \in B
	}

	a = marginal_order - b_; // ojo: h[current var] = 0 => current_var \in A
	//split the ESP by variable current_var:
	//esp -> esp_no_a + esp_yes_a
	//since esp in the next iteration we be esp_no_a of the current iter (esp=esp_no_a) then
	//esp_no_a is directly stored in esp
	esp_red_yes_a_[1] = t_sampling_[current_var];
	for (int k = 1; k < num_vars; k++)
	{
		esp_red_[k] = esp_red_[k] - esp_red_yes_a_[k];
		esp_red_yes_a_[k + 1] = esp_red_[k] * t_sampling_[current_var];
		res += g_n_[n_ - a - k][b_] * esp_red_[k];
	}
	esp_red_[num_vars] = esp_red_[num_vars] - esp_red_yes_a_[num_vars];
	res += g_n_[n_ - a][b_]; //* esp_red_[0]; // iter k=0
	if (num_vars != 0)
		res += g_n_[n_ - a - num_vars][b_] * esp_red_[num_vars]; // iter k= num_vars

	result = (long double)exp(-theta_acum_not_in_A + theta[current_var]) * res;
	if (result < 0)
	{
		//cout<<"ERROR Negative marginal probability, maybe theta is too large?"<<endl;
		//exit(0);
	}
	return result;
}

void HammingFunctions::elementary_symmetric_polynomial(double *theta, int n,
													   long double *theta_exp_aux, long double **esp_aux, long double *esp)
{
	//esp[j][n]: j-th elementarySymmetricPolynomials of n items
	//theta_exp_aux , esp_aux: are defined outside because the learning process (NewtonRaphson) calls this func lots of return;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= n; j++)
			esp_aux[i][j] = 0;

		theta_exp_aux[i + 1] = (long double)exp(theta[i]) - 1;
	}

	for (int j = 0; j <= n; j++)
		esp_aux[n][j] = 0;
	for (int j = 1; j <= n; j++)
		for (int k = 1; k <= j; k++)
			esp_aux[1][j] += theta_exp_aux[k]; //la suma de los primeros
	for (int i = 2; i <= n; i++)
		for (int j = i; j <= n; j++)
			esp_aux[i][j] = esp_aux[i][j - 1] + esp_aux[i - 1][j - 1] * theta_exp_aux[j]; //theta va de 0..n-1 y esp va de 1..n
	esp[0] = 1;

	for (int i = 1; i < n + 1; i++)
		esp[i] = esp_aux[i][n];
}

void HammingFunctions::split_elementary_symmetric_polynomial(long double *esp,
															 double *theta, int n, long double **esp_no_a, long double **esp_yes_a)
{
	for (int k = 0; k <= n; k++)
	{
		for (int i = 0; i < n; i++)
		{
			esp_no_a[k][i] = 0;
			esp_yes_a[k][i] = 0;
		}
	}
	for (int i = 0; i < n; i++)
	{
		esp_no_a[0][i] = 1;
		esp_yes_a[0][i] = 1; //default
		esp_yes_a[1][i] = (exp(theta[i]) - 1);
	}
	for (int k = 1; k < n; k++)
		for (int i = 0; i < n; i++)
		{
			esp_no_a[k][i] = esp[k] - esp_yes_a[k][i];
			esp_yes_a[k + 1][i] = esp_no_a[k][i] * (exp(theta[i]) - 1);
		}
	for (int i = 0; i < n; i++)
		esp_no_a[n][i] = esp[n] - esp_yes_a[n][i];
}

double HammingFunctions::psi_whm(double *theta)
{
	long double res = 0, sum_theta = 0;
	//long double *esp = new long double[n_ + 1];
	for (int k = 0; k < n_; k++)
		sum_theta += theta[k];

	//#for (int i = 0 ; i< n_ + 1; i++) {
	//#for (int j = 0 ; j< n_ + 1; j++) std::cout <<  aux_esp_[i][j] << "\n";
	//#}

	elementary_symmetric_polynomial(theta, n_, t_, aux_esp_, esp_);

	for (int k = 0; k <= n_; k++)
		res += facts_[n_ - k] * esp_[k];

	//delete[] esp_;
	return (res * exp(-sum_theta));
}

bool HammingFunctions::is_theta_array_uniform(double *theta_array)
{

	for (int i = 1; i < n_; i++)
	{
		if (theta_array[0] + 0.00001 < theta_array[i] || theta_array[0] - 0.00001 > theta_array[i])
		{
			return false;
		}
	}
	return true;
}

void HammingFunctions::multistage_sampling(double theta, int *sample)
{
	random_shuffle_sampling_MM(theta, sample);
}

void HammingFunctions::compose(int *s1, int *s2, int *res)
{
	for (int i = 0; i < n_; i++)
	{
		res[i] = s1[s2[i] - 1];
	}
}

void HammingFunctions::multistage_sampling_consensus(double theta, int *sigma_0, int *sample)
{

	multistage_sampling(theta, e_centered_sample_);

	compose(e_centered_sample_, sigma_0, sample);

	// if (fastfitness_can_be_used)
	// {

	// 	// cout << "-" << endl;
	// 	// PrintPythonArray(h_, n_);
	// 	// PrintPythonArray(e_centered_sample_, n_);
	// 	// cout << "sigma_0";
	// 	// PrintPythonArray(sigma_0, n_);
	// 	// cout << "Pos1 " << swap_pos1 << ", pos2" << swap_pos2 << endl;
	// 	// cout << "composed sample: ";
	// 	// PrintPythonArray(sample, n_);
	// 	// cout << "-" << endl;

	// 	last_fitness_value_from_fastfitness_delta = qap->update_fitness_on_swap(sigma_0, sample);
	// 	//cout << "--> " << last_fitness_value_from_fastfitness_delta << "<-- DELTA" << endl;
	// }
}

// find sigma mode permutation mm
void HammingFunctions::estimate_consensus_exact_mm(int m, int **samples, int *sigma_0)
{
	Lap lap;

	for (int i = 0; i < n_; i++)
	{
		for (int j = 0; j < n_; j++)
			freq_[i][j] = 0;
	}

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n_; j++)
			freq_[j][samples[i][j] - 1]--;

	//int cost = -1 * lap.lap(n_, freq, rows, cols, u, v);
	lap.lap(n_, freq_, rows_, cols_, u_, v_);
	for (int i = 0; i < n_; i++)
		sigma_0[i] = rows_[i] + 1;
}

#define NR_END 1
#define FREE_ARG char *

double *vector_hf(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v = (double *)malloc((int)((nh - nl + 1 + NR_END) * sizeof(double)));
	if (!v)
	{
		//ERRORchar str_msg[100];strcpy(str_msg, "allocation failure in vector()");nrerror(str_msg);
	}
	return v - nl + NR_END;
}

void HammingFunctions::free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}

#define TOL 2.0e-4

double (*nrfunc_qap)(double[]);
void (*nrdfun_qap)(double[], double[]);

int ncom_qap = 0;
double *pcom_qap = 0, *xicom_qap = 0;

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a, b, c, d) \
	(a) = (b);           \
	(b) = (c);           \
	(c) = (d);

void HammingFunctions::mnbrak(double *ax, double *bx, double *cx, double *fa,
							  double *fb, double *fc)
{

	//double (HammingFunctions::*)(double)
	//double (*)(double)

	double ulim, u, r, q, fu, dum;

	*fa = f1dim(*ax);
	*fb = f1dim(*bx);
	/*if (*fb > *fa && *ax == 0 ) {
	 *ax = *bx = *cx = 0.1 ;
	 return ;
	 }*/
	if (*fb > *fa)
	{
		SHFT(dum, *ax, *bx, dum)
		SHFT(dum, *fb, *fa, dum)
	}
	*cx = (*bx) + GOLD * (*bx - *ax);
	*fc = f1dim(*cx);
	while (*fb > *fc)
	{
		r = (*bx - *ax) * (*fb - *fc);
		q = (*bx - *cx) * (*fb - *fa);
		u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) / (2.0 * SIGN(MAX(fabs(q - r), TINY), q - r));
		ulim = (*bx) + GLIMIT * (*cx - *bx);
		if ((*bx - u) * (u - *cx) > 0.0)
		{
			fu = f1dim(u);
			if (fu < *fc)
			{
				*ax = (*bx);
				*bx = u;
				*fa = (*fb);
				*fb = fu;
				return;
			}
			else if (fu > *fb)
			{
				*cx = u;
				*fc = fu;
				return;
			}
			u = (*cx) + GOLD * (*cx - *bx);
			fu = f1dim(u);
		}
		else if ((*cx - u) * (u - ulim) > 0.0)
		{
			fu = f1dim(u);
			if (fu < *fc)
			{
				SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
				SHFT(*fb, *fc, fu, f1dim(u))
			}
		}
		else if ((u - ulim) * (ulim - *cx) >= 0.0)
		{
			u = ulim;
			fu = f1dim(u);
		}
		else
		{
			u = (*cx) + GOLD * (*cx - *bx);
			fu = f1dim(u);
		}
		SHFT(*ax, *bx, *cx, u)
		SHFT(*fa, *fb, *fc, fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT

#define ITMAX 100
#define ZEPS 1.0e-10
#define SIGN(a, b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MOV3(a, b, c, d, e, f) \
	(a) = (d);                 \
	(b) = (e);                 \
	(c) = (f);

double HammingFunctions::dbrent(double ax, double bx, double cx, double tol,
								double *xmin)
{
	int iter, ok1, ok2;
	double a, b, d, d1, d2, du, dv, dw, dx, e = 0.0;
	double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

	//&f1dim,&df1dim

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	fw = fv = fx = f1dim(x);  //f1
	dw = dv = dx = df1dim(x); //df1

	for (iter = 1; iter <= ITMAX; iter++)
	{
		xm = 0.5 * (a + b);
		tol1 = tol * fabs(x) + ZEPS;
		tol2 = 2.0 * tol1;
		if (fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
		{
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1)
		{
			d1 = 2.0 * (b - a);
			d2 = d1;
			if (dw != dx)
				d1 = (w - x) * dx / (dx - dw);
			if (dv != dx)
				d2 = (v - x) * dx / (dx - dv);
			u1 = x + d1;
			u2 = x + d2;
			ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
			ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
			olde = e;
			e = d;
			if (ok1 || ok2)
			{
				if (ok1 && ok2)
					d = (fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d = d1;
				else
				{
					d = d2;
				}
				if (fabs(d) <= fabs(0.5 * olde))
				{
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = SIGN(tol1, xm - x);
				}
				else
				{
					d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
				}
			}
			else
			{
				d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
			}
		}
		else
		{
			d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
		}
		if (fabs(d) >= tol1)
		{
			u = x + d;
			fu = f1dim(u);
		}
		else
		{
			u = x + SIGN(tol1, d);
			fu = f1dim(u);
			if (fu > fx)
			{
				*xmin = x;
				return fx;
			}
		}
		du = df1dim(u);
		if (fu <= fx)
		{
			if (u >= x)
				a = x;
			else
				b = x;
			MOV3(v, fv, dv, w, fw, dw)
			MOV3(w, fw, dw, x, fx, dx)
			MOV3(x, fx, dx, u, fu, du)
		}
		else
		{
			if (u < x)
				a = u;
			else
				b = u;
			if (fu <= fw || w == x)
			{
				MOV3(v, fv, dv, w, fw, dw)
				MOV3(w, fw, dw, u, fu, du)
			}
			else if (fu < fv || v == x || v == w)
			{
				MOV3(v, fv, dv, u, fu, du)
			}
		}
	}
	//ERROR char str_msg[100];strcpy(str_msg, "Too many iterations in routine DBRENT");nrerror(str_msg);
	return 0;
}

#undef ITMAX
#undef ZEPS
#undef SIGN
#undef MOV3

double HammingFunctions::expectation(long double theta, int n)
{
	long double x_n = 0, x_n_1 = 0, aux = 0;
	for (int k = 0; k <= n; k++)
	{
		aux = powl(expl(theta) - 1, k) / facts_[k];
		x_n += aux;
		if (k < n)
			x_n_1 += aux; //pow (exp(theta )-1, k) / facts_[ k ];
	}
	return (double)(n * x_n - x_n_1 * expl(theta)) / x_n;
}

double HammingFunctions::f1dim(double x)
{
	int j;
	double f, *xt; //,*vector();

	xt = vector_hf(1, ncom_qap);
	for (j = 1; j <= ncom_qap; j++)
		xt[j] = pcom_qap[j] + x * xicom_qap[j];
	f = likeli_wmh(xt);
	free_vector(xt, 1, ncom_qap);
	return f;
}

double HammingFunctions::df1dim(double x)
{
	int j;
	double df1 = 0.0;
	double *xt, *df; //,*vector();
	//void free_vector();

	xt = vector_hf(1, ncom_qap);
	df = vector_hf(1, ncom_qap);
	for (j = 1; j <= ncom_qap; j++)
		xt[j] = pcom_qap[j] + x * xicom_qap[j];
	dlikeli_wmh(xt, df);
	for (j = 1; j <= ncom_qap; j++)
		df1 += df[j] * xicom_qap[j];
	free_vector(df, 1, ncom_qap);
	free_vector(xt, 1, ncom_qap);

	return df1;
}

void HammingFunctions::dlinmin(double p[], double xi[], int n, double *fret)
{
	int j;
	double xx, xmin, fx, fb, fa, bx, ax;
	ncom_qap = n;
	pcom_qap = vector_hf(1, n);
	xicom_qap = vector_hf(1, n);
	//nrfunc_qap=likeli_wmh;
	//nrdfun_qap=dlikeli_wmh;
	for (j = 1; j <= n; j++)
	{
		pcom_qap[j] = p[j];
		xicom_qap[j] = xi[j];
	}
	ax = 0.0;
	xx = 1.0;
	bx = 2.0;

	mnbrak(&ax, &xx, &bx, &fa, &fx, &fb);
	*fret = dbrent(ax, xx, bx, TOL, &xmin);
	for (j = 1; j <= n; j++)
	{
		xi[j] *= xmin;
		p[j] += xi[j];
		//if (p[ j ] < 0 ) p[j] = 0 ;//test1
		//if (p[ j ] > 10 ) p[j] = 10 ;
	}
	free_vector(xicom_qap, 1, n);
	free_vector(pcom_qap, 1, n);
}

#undef TOL

#include <math.h>
#define ITMAX 100
#define EPS1 1.0e-10
#define FREEALL            \
	free_vector(xi, 1, n); \
	free_vector(h, 1, n);  \
	free_vector(g, 1, n);

void HammingFunctions::frprmn(double p[], int n, double ftol, int *iter,
							  double *fret)
{

	//Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a function func, using its gradient as calculated by a routine dfunc. The convergence tolerance on the function value is input as ftol. Returned quantities are p (the location of the minimum), iter (the number of iterations that were performed), and fret (the minimum value of the function). The routine linmin is called to perform line minimizations.
	int j, its;
	long double gg, gam, fp, dgg;
	double *g, *h, *xi; //,*vector_hf();
	g = vector_hf(1, n);
	h = vector_hf(1, n);
	xi = vector_hf(1, n);
	fp = likeli_wmh(p);
	dlikeli_wmh(p, xi);

	for (j = 1; j <= n; j++)
	{
		g[j] = -xi[j];
		xi[j] = h[j] = g[j];
	}

	for (its = 1; its <= ITMAX; its++)
	{
		//cout<<"p[i]: ";for (j=1;j<=n;j++) cout <<p[j]<<" ";cout<<" point from frprmn"<<endl;
		//cout<<"h_avg: ";for (j=1;j<=n;j++) cout <<h_avg_[j-1]<<" ";cout<<" h from frprmn"<<endl;
		*iter = its;
		dlinmin(p, xi, n, fret);
		//for (int i = 0 ; i < n ; i++) if (p[ i +1 ]<0) p[i+1] = 0 ;
		//linmin(p, xi, n, fret, func);
		if (2.0 * fabs(*fret - fp) <= ftol * (fabs(*fret) + fabs(fp) + EPS1))
		{
			FREEALL
			return;
		}
		fp = likeli_wmh(p);
		dlikeli_wmh(p, xi);
		dgg = gg = 0.0;
		for (j = 1; j <= n; j++)
		{
			gg += g[j] * g[j];
			//dgg += xi[j]*xi[j];	//or this or the next(see numerical receipes)
			dgg += (xi[j] + g[j]) * xi[j];
		}
		if (gg == 0.0)
		{
			FREEALL
			return;
		}
		gam = dgg / gg;
		for (j = 1; j <= n; j++)
		{
			g[j] = -xi[j];
			xi[j] = h[j] = g[j] + gam * h[j];
		}
	}
	//char str_msg[100];strcpy(str_msg, "Too many iterations in FRPRMN");nrerror(str_msg);
	//ERROR
}

#undef ITMAX
#undef EPS1
#undef FREEALL

double HammingFunctions::likeli_wmh(double x[])
{
	//x is a vector_hf from 1..n
	//both  likeli_wmh and dlikeli_wmh return the result *(-1) because frprmn is for minimization and we need maximization
	long double like = 0, aux1 = 0, aux2 = 0, sum_theta = 0;
	bool penalty = false;
	for (int i = 0; i < n_; i++)
	{
		do_x[i] = (long double)x[i + 1];
		sum_theta += x[i + 1];
		if (do_x[i] < 0 || do_x[i] > 10)
			penalty = true;
	}
	elementary_symmetric_polynomial(do_x, n_, t_, aux_esp_, esp_);
	for (int i = 0; i < n_; i++)
		aux1 += do_x[i] * h_avg_[i];

	for (int k = 0; k <= n_; k++)
		aux2 += facts_[n_ - k] * esp_[k];
	aux2 = aux2 * exp(-sum_theta); //psi
	like = -m_ * (aux1 + log(aux2));

	if (like != like || penalty)
	{ //|| penalty != 0
		return DBL_MAX;
	}
	return -like; //* (penalty+1);
}

void HammingFunctions::dlikeli_wmh(double x[], double deriv[])
{
	//derivada de la verosimilitud
	//x y deriv [1..n]
	long double psi = 0;
	double sum_theta = 0, aux = 0;
	for (int i = 0; i < n_; i++)
	{
		do_x[i] = (double)x[i + 1];
		sum_theta += x[i + 1];
	}
	elementary_symmetric_polynomial(do_x, n_, t_, aux_esp_, esp_);
	split_elementary_symmetric_polynomial(esp_, do_x, n_, esp_no_a_,
										  esp_yes_a_);
	psi = 0;
	for (int k = 0; k <= n_; k++)
		psi += facts_[n_ - k] * esp_[k];
	psi = psi * exp(-sum_theta); //psi
	for (int i = 0; i < n_; i++)
	{
		psi_der[i] = 0;
		aux = 0;
		for (int k = 1; k <= n_; k++)
			aux += facts_[n_ - k] * esp_no_a_[k - 1][i];
		aux = aux * exp(-sum_theta + do_x[i]); //psi
		psi_der[i] = -psi + aux;
		//both  likeli_wmh and dlikeli_wmh return the result *(-1) because frprmn is for minimization and we need maximization
		deriv[i + 1] = (-1) * (psi_der[i] / psi + (long double)h_avg_[i]);
	}
}

void HammingFunctions::mle_theta_weighted_mallows_hamming(int m, double *h_avg,
														  double *theta)
{
	//OJO : Given a starting point p[1..n]!!!
	//minimize a function on multiple variables

	m_ = m;
	h_avg_ = h_avg;
	//double a1 = 0;

	int iter = 0;
	double fmin = 0;

	for (int i = 0; i < n_; i++)
	{
		point[i + 1] = 0.2;
	}
	for (int j = 0; j < n_; j++)
	{

		// if (h_avg_[ j ] == 0) h_avg_[ j ] = 0.001;//los ceros dan problemas, txapu
		// if (h_avg_[ j ] >= 1) h_avg_[ j ] =  0.99;
	}

	frprmn(point, n_, 0.0001, &iter, &fmin);

	for (int j = 0; j < n_; j++)
	{
		theta[j] = (double)point[j + 1];
	}
}

// find sigma_0 mode permutation
void HammingFunctions::estimate_consensus_approx_gmm(int m, int **samples, int *sigma_0)
{

	Lap lap;
	double a1 = 0;

	for (int i = 0; i < n_; i++)
	{
		for (int j = 0; j < n_; j++)
			freq_[i][j] = 0;
	}

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n_; j++)
			freq_[j][samples[i][j] - FIRST_ITEM]--; //for the lap, minimize

	//freq[i,j]=sample[i,j]

	//int cost = -1 * lap.lap(n_, freq, rows, cols, u, v);
	lap.lap(n_, freq_, rows_, cols_, u_, v_);
	for (int i = 0; i < n_; i++)
	{
		sigma_0[i] = rows_[i] + FIRST_ITEM;
		sigma_0_inv_[(rows_[i] + FIRST_ITEM) - FIRST_ITEM] = i + FIRST_ITEM;
	}

	for (int i = 0; i < n_; i++)
		for (int j = 0; j < n_; j++)
			freq_[i][j] = freq_[i][j] * -1;
	for (int i = 0; i < n_; i++)
	{
		h_avg_[i] = (double)1 - (double)freq_[sigma_0_inv_[i] - FIRST_ITEM][i] / m;
	}
	mle_theta_weighted_mallows_hamming(m, h_avg_, theta_);
	for (int i = 0; i < n_; i++)
		a1 += h_avg_[i] * theta_[i];
}

// find theta on mm
int HammingFunctions::distance_to_sample(int **samples, int m, int *sigma)
{
	int dist = 0;
	for (int s = 0; s < m; s++)
	{
		for (int i = 0; i < n_; i++)
			if (samples[s][i] != sigma[i])
				dist++;
	}
	return dist;
}

// * Theta parameter estimation function.
double HammingFunctions::f(double theta)
{

	long double x_j = 0, sum_to_n = 0, sum_to_n_1 = 0, psi = 0, psi_der = 0;

	for (int j = 0; j <= n_; j++)
	{

		x_j = (long double)pow(exp(theta) - 1, j) / facts_[j];
		if (j < n_)
			sum_to_n_1 += x_j;
		sum_to_n += x_j;
		if (sum_to_n > DBL_MAX || sum_to_n != sum_to_n)
			return DBL_MAX;
	}
	psi = facts_[n_] * exp(-theta * n_) * sum_to_n;
	psi_der = -n_ * psi + facts_[n_] * exp(theta * (1 - n_)) * sum_to_n_1;

	//cout << "psi : " << psi << " psi_der: " << psi_der<<endl;
	double f_fligner = (double)(psi_der / psi + dist_avg_);
	//cout<<"theta "<<theta<<" dist "<<dist_avg_<<" fu "<<f_fligner<<endl;
	if (f_fligner != f_fligner)
		f_fligner = 0; //trace
	return f_fligner;
}

// Theta parameter estimation function derivation.
double HammingFunctions::fdev(double theta)
{

	if (model_ == 0)
	{
		long double x_j = 0, sum_to_n = 0, sum_to_n_1 = 0, sum_to_n_2 = 0, psi = 0, psi_der = 0, psi_der_2 = 0;

		for (int j = 0; j <= n_; j++)
		{
			x_j = (long double)pow(exp(theta) - 1, j) / facts_[j];
			if (j < n_ - 1)
				sum_to_n_2 += x_j;
			if (j < n_)
				sum_to_n_1 += x_j;
			sum_to_n += x_j;
			if (sum_to_n_1 > DBL_MAX || sum_to_n_1 != sum_to_n_1)
				return DBL_MAX - 1;
		}
		psi = facts_[n_] * exp(-theta * n_) * sum_to_n;
		psi_der = -n_ * psi + facts_[n_] * exp(theta * (1 - n_)) * sum_to_n_1;
		psi_der_2 = -n_ * psi_der + facts_[n_] * exp(theta * (1 - n_)) * ((1 - n_) * sum_to_n_1 + sum_to_n_2);
		double res = (double)(-psi_der_2 * psi_der - psi_der * psi_der) / (psi * psi);
		//cout<<"theta "<<theta<<" dist "<<dist_avg_<<"ps... "<<psi<<" "<<psi_der<<" "<<psi_der_2<<" fd "<<res<<endl;
		if (res > DBL_MAX)
			return DBL_MAX - 1;
		return res;
	}
	else if (model_ == 1)
	{
		cout
			<< "Solve Weigthed Hamming Mallows with Newton_raphson::mle_theta_weighted_mallows_hamming "
			<< endl;
		exit(1);
	}
	return 0;
}

void HammingFunctions::funcd(double theta, double *ff, double *ffdev)
{
	//*ff = f(theta, n, Vjs);
	//*ffdev = fdev(theta, n);

	*ff = f(theta);

	*ffdev = fdev(theta);
}

// Newton - Rapshon execution algorithm.
double HammingFunctions::rtsafe(double x1, double x2, double xacc)
{
	int j;
	double dx, dxold;
	double temp, xh, xl, rts;
	double f, df, fl, fh;

	funcd(x1, &fl, &df); //params,: theta, f, fdev
	funcd(x2, &fh, &df);
	//if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
	//	cout<<"Root must be bracketed in rtsafe"<<endl;
	if (fl == 0.0)
		return x1;
	if (fh == 0.0)
		return x2;
	if (fl < 0.0)
	{
		xl = x1;
		xh = x2;
	}
	else
	{
		xh = x1;
		xl = x2;
	}
	rts = 0.5 * (x1 + x2);
	rts = x1; //<-fijamos un valor inicial para el theta.
	dxold = fabs(x2 - x1);
	dx = dxold;
	funcd(rts, &f, &df);
	//cout << "\n rts: " << rts << " f: " << f << " df: " << df << endl;
	for (j = 1; j <= MAXIT; j++)
	{
		//cout<<"rts: "<<rts<<". Function val: "<<f<<" f_dev "<<df<<endl;
		//Initialize the guess for root, the “stepsize before last,” and the last step.
		//Loop over allowed iterations.
		if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) || (fabs(2.0 * f) > fabs(dxold * df)))
		{
			//Bisect if Newton out of range, //or not decreasing fast enough.
			dxold = dx;
			dx = 0.5 * (xh - xl);
			rts = xl + dx;
			if (xl == rts)
				return rts;
		}
		else
		{
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts)
				return rts;
			//Change in root is negligible. Newton step acceptable. Take it.
		}
		//cout<<"DX: "<<dx<<endl;
		if (fabs(dx) < xacc)
		{
			return rts;
		}
		funcd(rts, &f, &df); //The one new function evaluation per iteration.
		//Orient the search so that f(xl) < 0.
		//Convergence criterion.

		if (f < 0.0) //Maintain the bracket on the root.
			xl = rts;
		else
			xh = rts;

		//cout<<" rts: "<<rts<<" Function val: "<<f<<" f_dev "<<df<<endl;
	}
	//  cout << "Maximum number of iterations exceeded in rtsafe" << endl;

	return 0.0; //Never get here.
}

double HammingFunctions::Newton_raphson_method_qap(double dAvg_val,
											   double initialGuess, int model, int j_index, long double *count)
{
	//count_ = count;//count of num permus for ulam

	//    n_ = n_val;
	//j_index_ = j_index; // for kendall GMM
	model_ = model;
	dist_avg_ = dAvg_val;
	//distance_id_=distanceModel_val;
	double theta;
	double xacc = 0.000001;

	//cout << "initial guess: " <<initialGuess << "  upper_theta: " << UPPER_THETA << " acc: " << xacc << endl;
	//TEST //for (double i = -5.1 ; i< 5.1 ; i++) cout<<"Theta "<<i<< " f: "<<f(i)<<" dev: "<<fdev(i)<<endl;
	//cout << "dist_average: " << dist_avg_ << endl;
	theta = rtsafe(initialGuess, UPPER_THETA, xacc);
	//cout<<"newton "<<j_index_<<" "<<dist_avg_<<" "<<theta<<endl;
	return theta;
}

double HammingFunctions::estimate_theta_mm(int m, int *sigma_0, int **samples)
{

	//Newton_raphson nr(n_);
	double dist = 0;

	dist = distance_to_sample(samples, m, sigma_0);

	int model = 0; // mm --> model = 0
	model_ = 0;

	return this->Newton_raphson_method_qap((double)dist / m, 0.0, model, -1, NULL);
}

// find theta on gmm

void HammingFunctions::sample_to_h_vector(int **samples, int m, int *sigma,
										  double *h_avg)
{
	//if sigma != NULL => h( samples sigma^{-1} )
	for (int i = 0; i < n_; i++)
		h_avg[i] = 0;
	for (int s = 0; s < m; s++)
		for (int i = 0; i < n_; i++)
			if (sigma == NULL)
			{
				if (samples[s][i] != i + 1)
					h_avg[i]++;
			}
			else
			{ //right compose with the inverse of the sample
				//h_j = 0  <=> sigma^{-1}(j) = sigma_0^{-1}(j) <=> sigma(i) = sigma_0(i) = j
				if (sigma[i] != samples[s][i])
					h_avg[samples[s][i] - 1]++;
			}
	for (int i = 0; i < n_; i++)
		h_avg[i] = (double)h_avg[i] / m;
}

void HammingFunctions::estimate_theta_gmm(int m, int *sigma_0, int **samples,
										  double *theta)
{
	//Newton_raphson nr(n_);
	//Generic gen;

	//double  *h_avg         = new double [ n_ ];

	sample_to_h_vector(samples, m, sigma_0, h_avg_);

	mle_theta_weighted_mallows_hamming(m, h_avg_, theta);

	//for ( int j = 0 ; j < n_ ; j ++) cout<< theta[j]<<" "; cout<<" theta (estimate_theta_gmm) "<<endl;
	//delete [] h_avg;
}

int HammingFunctions::which_distance_to_sample_MM(double theta)
{
	if (last_theta - 0.000001 > theta || last_theta + 0.000001 < theta) // recalculate probabilities
	{
		last_theta = theta;
		compute_normalization_constant((long double)theta);
		compute_sampling_probabilities((long double)theta);
	}
	return choose_index_given_probabilities(sampling_probabilities, n_);
}
void HammingFunctions::random_shuffle_sampling_MM(double theta, int *sample)
{
	fastfitness_can_be_used = false;

	int n_of_deranged_positions = which_distance_to_sample_MM(theta);

	//cout << n_of_deranged_positions << " ";

	for (int i = 0; i < n_of_deranged_positions; i++)
	{
		h_[i] = 1;
	}

	for (int i = n_of_deranged_positions; i < n_; i++)
	{
		h_[i] = 0;
	}

	//  Yates shuffle
	for (int i = n_ - 1; i > 0; --i)
	{
		int r_index = random_range_integer_uniform(i + 1);
		int temp = h_[i];
		h_[i] = h_[r_index];
		h_[r_index] = temp;
	}

	dist_decomp_vector2perm(h_, sample);

	if (n_of_deranged_positions == 2)
	{
		fastfitness_can_be_used = true;

		// get deranged positions
		// for (int i = 0; i < n_; i++)
		// {
		// 	if (h_[i] == 1)
		// 	{
		// 		swap_pos1 = i;
		// 		for (int j = i + 1; j < n_; j++)
		// 		{
		// 			if (h_[j] == 1)
		// 			{
		// 				swap_pos2 = j;
		// 				break;
		// 			}
		// 		}
		// 		break;
		// 	}
		// }
		//cout << endl;
		//PrintPythonArray(h_, n_);
		//cout << swap_pos1 << " - " << swap_pos2 << endl;
	}
	//cout << "-" << endl;
	//PrintPythonArray(h_, n_);
	//PrintPythonArray(sample, n_);
	//cout << "Pos1 " << swap_pos1 << ", pos2" << swap_pos2 << endl;
	//cout << "-" << endl;
}

void HammingFunctions::compute_normalization_constant(long double theta)
{
	//compute first term
	long double first_term = facts_[n_] * expl(-theta * (long double)n_);

	//compute second term
	long double exp_theta_minus_one = expl(theta) - (long double)1;
	long double second_term = 1; // second_term[0] = (exp - 1)^0 / 0! 0 1
	long double last_second_term_sumand = 1;

	for (int k = 1; k <= n_; k++)
	{
		last_second_term_sumand *= exp_theta_minus_one;
		last_second_term_sumand /= k;
		second_term += last_second_term_sumand;
	}

	//multiply two terms
	normalization_constant = first_term * second_term;
	//cout << endl;
	//cout << "n_of_permus_at_distance_k: " << n_of_permus_at_distance_k << endl;
	//PrintPythonArray(n_of_permus_at_distance_k, n_ + 1);
	//cout << endl;
}

void HammingFunctions::compute_sampling_probabilities(long double theta)
{

	for (int i = 0; i <= n_; i++)
	{
		sampling_probabilities[i] = (double)(n_of_permus_at_distance_k[i] / normalization_constant * expl(-(theta * (long double)i)));
	}

	if (not ALLOW_0_DISTANCE_SAMPLING)
	{
		sampling_probabilities[0] = 0;
	}
	make_probabilities_sum_1(sampling_probabilities);
}

void HammingFunctions::make_probabilities_sum_1(double *sampling_probabilities)
{
	double sum_of_probs = 0;

	for (int i = 0; i <= n_; i++)
	{
		sum_of_probs += sampling_probabilities[i];
	}

	for (int i = 0; i <= n_; i++)
	{
		sampling_probabilities[i] = sampling_probabilities[i] / sum_of_probs;
	}
}

HammingFunctions::HammingFunctions(int n, QAP *qap)
{

	n_ = n;
	model_ = 1;
	esp_red_ = new long double[n_ + 1];
	t_sampling_ = new long double[n_];
	t_ = new long double[n_ + 1];
	aux_esp_ = new long double *[n_ + 1];
	esp_ini_ = new long double[n_ + 1];
	esp_red_yes_a_ = new long double[n_ + 1];
	g_n_ = new long double *[n_ + 1];
	this->qap = qap;

	for (int i = 0; i < n_ + 1; i++)
	{
		aux_esp_[i] = new long double[n_ + 1];
		for (int j = 0; j < n_ + 1; j++)
		{
			aux_esp_[i][j] = 0;
		}
	}

	//initialize factorials
	facts_ = new long double[n + 1];
	facts_[0] = 1;
	for (int i = 1; i <= n; i++)
	{
		facts_[i] = facts_[i - 1] * i;
	}

	deran_num_ = new double[n_ + 1];
	deran_num_[0] = 1;
	deran_num_[1] = 0;
	for (int i = 2; i <= n_; i++)
		deran_num_[i] = deran_num_[i - 1] * (i - 1) + deran_num_[i - 2] * (i - 1);

	for (int i = 0; i < n_ + 1; i++)
	{
		g_n_[i] = new long double[n_ + 1];
		for (int j = 0; j <= i; j++)
		{
			g_n_[i][j] = count_permus_with_at_least_k_unfixed_points(i, j);
		}
		for (int j = 0; j <= n_; j++)
		{
			aux_esp_[i][j] = 0;
		}
	}

	//whi initialization

	esp_ = new long double[n_ + 1];
	esp_no_a_ = new long double *[n_ + 1];
	esp_yes_a_ = new long double *[n_ + 1];
	//aux_esp_  = new long double *[ n_ + 1 ];//esp
	//t_        = new long double  [ n_ + 1 ];//esp
	for (int k = 0; k <= n_; k++)
	{
		esp_no_a_[k] = new long double[n_];
		esp_yes_a_[k] = new long double[n_];
		//aux_esp_[ k ] = new long double[ n_ + 1 ];
		for (int i = 0; i < n_; i++)
		{
			esp_no_a_[k][i] = 0;
			esp_yes_a_[k][i] = 0;
			aux_esp_[k][i] = 0;
		}
	}

	//int iter = 0;
	//double fmin = 0;

	point = new double[n_ + 1];
	dpoint = new double[n_ + 1];

	b_ = 0;

	do_x = new double[n_];
	psi_der = new long double[n_];

	deran_ = new int[n_ - 2];
	conv_ = new int[n_ - 1];

	ran_ = new int[n_];

	h_ = new int[n_];

	e_centered_sample_ = new int[n_];

	freq_ = new int *[n_];
	for (int i = 0; i < n_; i++)
	{
		freq_[i] = new int[n_];
	}
	rows_ = new int[n_];
	cols_ = new int[n_];
	u_ = new int[n_];
	v_ = new int[n_];

	sigma_0_inv_ = new int[n_];
	sigma_0_inv_neig_best_ = new int[n_];
	theta_ = new double[n_];
	for (int i = 0; i < n_; i++)
	{
		for (int j = 0; j < n_; j++)
		{
			freq_[i][j] = 0;
		}
	}

	h_avg_ = new double[n_];
	deranged_positions_ = new int[n_];

	sampling_probabilities = new double[n_ + 1];

	n_of_permus_at_distance_k = new long double[n_ + 1];
	long double first_21_values[21] = {1, 0, 1, 2, 9, 44, 265, 1854, 14833, 133496, 1334961, 14684570, 176214841, 2290792932, 32071101049, 481066515734, 7697064251745, 130850092279664, 2355301661033953, 44750731559645106, 895014631192902121};
	for (int i = 0; i <= n; i++)
	{
		if (i < 21)
		{
			n_of_permus_at_distance_k[i] = facts_[n_] / facts_[n_ - i] / facts_[i] * (long double)first_21_values[i];
		}
		else
		{
			// facts[i] simplified.
			n_of_permus_at_distance_k[i] = facts_[n_] / facts_[n_ - i] / expl((long double)1);
		}
	}

	GMM_consensus = new int[n_];
}

void HammingFunctions::sample_with_kernel(int *kernel_consensus, double theta, int *result_individual)
{
	multistage_sampling_consensus(theta, kernel_consensus, result_individual); // fastfitness is calculated here

}

HammingFunctions::~HammingFunctions()
{

	delete[] deran_num_;
	delete[] esp_red_;
	delete[] esp_ini_;
	delete[] esp_red_yes_a_;
	delete[] t_;
	delete[] t_sampling_;
	for (int i = 0; i < n_ + 1; i++)
	{
		delete[] aux_esp_[i];
		delete[] g_n_[i];
	}
	delete[] aux_esp_;
	delete[] g_n_;
	delete[] do_x;
	//delete [] removed;
	delete[] psi_der;
	delete[] point;
	delete[] dpoint;
	delete[] facts_;
	delete[] esp_;

	for (int k = 0; k <= n_; k++)
	{
		delete[] esp_no_a_[k];
		delete[] esp_yes_a_[k];
		//delete [] aux_esp_[ k ];
	}

	delete[] esp_no_a_;
	delete[] esp_yes_a_;
	delete[] deran_;
	delete[] conv_;
	delete[] ran_;
	delete[] h_;
	delete[] e_centered_sample_;
	delete[] rows_;
	delete[] cols_;
	delete[] u_;
	delete[] v_;
	for (int i = 0; i < n_; i++)
	{
		delete[] freq_[i];
	}
	delete[] freq_;
	delete[] theta_;
	delete[] sigma_0_inv_neig_best_;
	delete[] sigma_0_inv_;
	delete[] h_avg_;
	delete[] deranged_positions_;

	delete[] sampling_probabilities;
	delete[] n_of_permus_at_distance_k;
	delete[] GMM_consensus;
}
