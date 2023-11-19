/*
 * Population.h
 *
 *  Created on: Jun 27, 2018
 *      Author: paran
 */

#ifndef POPULATION_H_
#define POPULATION_H_

extern int MAX_ITERATIONS;

//#define START_MINIMUN_THETA 0.25
//#define FINAL_EXPECTATION 0.11

#include "QAP.h"
#include "SimpleLog.h"
#include "HammingFunctions.h"
#include "Parameters.h"

extern SimpleLog *LOG;

class Population
{

  public:
	Parameters *param;

	Population(Parameters *param);
	virtual ~Population();
	void sort();
	void calculate_fitness();
	virtual void sample() = 0;
	void get_best_permutaion(int *sigma);
	virtual void learn() = 0;

	int get_best_fitness() { return best_fitness; }
	double get_mean_fitness() { return Average(fitness_array, POPSIZE); }
	double get_fitness_variance() { return Variance(fitness_array, POPSIZE); }
	long int get_fitness_evals() { return fitness_evals; }
	int get_best_fit_change_count() { return best_fit_change_count; }
	int get_popsize() { return POPSIZE; }
	void local_search();
	void purge_repeated();

	int **pop;

	QAP *qap;
	HammingFunctions *ham;

	int n;

	void initialize_Population();
	void evaluate_individual_on_position_i(int i);
	void GMM_limit_theta_value();

	long int fitness_evals = 0;

	double pop_mult_factor = 0.0;
	int POPSIZE;
	int avoid_repetition_mode = 0;
	int best_fitness;
	int best_fit_change_count;
	int *best_permutation;
	int *best_permutation_one_position_swapped_local_search;
	int *fitness_array;
	bool *repeated_indexes;
	int *consensus;
	double *theta_array;
	double theta;
	double FINAL_THETA;

	//KEDAMM VARIABLES
	int ker_iteration = 0;
	int ker_pool_mode = 0;
	int ker_theta_increase_mode = 0;
	double ker_percentage_done = 0.0;
	bool ker_finished = false;
	bool *already_evaluated;
	double MIDDLE_THETA;

	//local search
	int *random_i_indexes;
  	int *random_j_indexes;

};

#endif /* POPULATION_H_ */
