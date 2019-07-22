#include <iostream>
#include "QAP.h"
#include <limits.h>

#ifndef _HAMMINGFUNCTIONS_H__
#define _HAMMINGFUNCTIONS_H__

using namespace std;

class HammingFunctions
{

  public:
	HammingFunctions(int n, QAP *qap);
	~HammingFunctions();

	void compose(int *s1, int *s2, int *res);
	void generate_random_permutation(int *sigma);
	void generate_random_sample(int m, int **sample);
	void multistage_sampling(double theta, int *sample);
	void multistage_sampling_consensus(double theta, int *sigma_0, int *sample);
	void estimate_consensus_exact_mm(int m, int **samples, int *sigma_0);
	void estimate_consensus_approx_gmm(int m, int **samples, int *sigma_0);
	int distance_to_sample(int **samples, int m, int *sigma);
	double estimate_theta_mm(int m, int *sigma_0, int **samples);
	void estimate_theta_gmm(int m, int *sigma_0, int **samples, double *theta);

	template <typename T>
	void print_permutation_sample(int m, T **sample)
	{
		for (int i = 0; i < m; i++)
		{
			cout << "( ";
			for (int j = 0; j < n_; j++)
			{
				if (sample[i][j] < 10)
				{
					cout << "0";
				}
				cout << sample[i][j] << " ";
			}
			cout << ")\n";
		}
	}

	void sample_with_kernel(int *kernel_consensus, double theta, int *result_individual);
	int *GMM_consensus;
	void make_probabilities_sum_1(double *theta_array);
	double expectation(double theta);

	bool fastfitness_can_be_used = false;
	int *swap_sigma1;
	int *swap_sigma2;
	int swap_pos1 = 0;
	int swap_pos2 = 0;

	//protected:
	QAP *qap;
	void random_shuffle_sampling_MM(double theta, int *sample);
	double last_theta = -10.0;
	double *sampling_probabilities;
	long double normalization_constant;
	long double *n_of_permus_at_distance_k;

	int which_distance_to_sample_MM(double theta);
	void compute_normalization_constant(long double theta);
	void compute_sampling_probabilities(long double theta);
	bool is_theta_array_uniform(double *theta_array);
	long double compute_marginal_iterative(int *h, double *theta, int marginal_order);
	void random_derangement(int n, int *sigma);
	void generate_permu_from_list(int *ran, int dist, int *sigma);
	long double count_permus_with_at_least_k_unfixed_points(int n, int k);
	void dist_decomp_vector2perm(int *vec, int *sigma);
	void elementary_symmetric_polynomial(double *theta, int n, long double *theta_exp_aux, long double **esp_aux, long double *esp);
	void split_elementary_symmetric_polynomial(long double *esp, double *theta, int n, long double **esp_no_a, long double **esp_yes_a);
	double psi_whm(double *theta);
	void free_vector(double *v, long nl, long nh);
	double dbrent(double ax, double bx, double cx, double tol, double *xmin);
	double f1dim(double x);
	double df1dim(double x);
	void dlinmin(double p[], double xi[], int n, double *fret);
	void frprmn(double p[], int n, double ftol, int *iter, double *fret);
	double likeli_wmh(double x[]);
	void dlikeli_wmh(double x[], double deriv[]);
	void mle_theta_weighted_mallows_hamming(int m, double *h_avg, double *theta);
	double f(double theta);
	double fdev(double theta);
	void funcd(double theta, double *ff, double *ffdev);
	double rtsafe(double x1, double x2, double xacc);
	double Newton_raphson_method(double dAvg_val, double initialGuess, int model, int j_index, long double *count);
	void sample_to_h_vector(int **samples, int m, int *sigma, double *h_avg);
	void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc);

	int *h_;
	int n_;
	int model_;
	long double *facts_;
	long double *t_sampling_;
	long double *t_;
	long double **aux_esp_;
	long double *esp_ini_;
	long double theta_acum_not_in_A;
	int b_;
	long double *esp_red_;		 //compute marginal
	long double *esp_red_yes_a_; //compute marginal
	long double **g_n_;			 //compute marginal
	double *deran_num_;			 // coutn the number of derangements of n items

	double *point;
	double *dpoint;

	double *do_x;
	long double *psi_der;

	// NR
	double *h_avg_;
	long double *esp_;
	long double **esp_no_a_;
	long double **esp_yes_a_;
	int m_;

	//NR mm theta
	double dist_avg_;
	double UPPER_THETA = 5;

	//estimate theta
	int *deran_;
	int *conv_;

	//  remove all news from code
	int *ran_;

	int *e_centered_sample_;

	int **freq_;
	int *rows_;
	int *cols_;
	int *u_;
	int *v_;

	int *sigma_0_inv_;
	int *sigma_0_inv_neig_best_;
	double *theta_;
	int *deranged_positions_;
};

#endif
