/*
 * KEDAMM.h
 *
 *  Created on: Jun 27, 2018
 *      Author: paran
 */

#ifndef KEDAMM_H_
#define KEDAMM_H_
#include "Population.h"
#include "QAP.h"
#include "Parameters.h"
#include <string.h>

extern SimpleLog *RESULT_LOG;

class KEDAMM : public Population
{

  public:
	KEDAMM(Parameters *param); //read from cin
	~KEDAMM();

	void init_ker_memory_pop();
	void learn();
	void sample();

	int pool_mode;

  private:
	int **ker_memory_pop;
	int N_OF_ITERATIONS;

	int t_wait = 0;

	int n_of_tops_without_improvement = 0;

	int *print_theta_consensus;

	void sample_individual_at_position_i(int i);

	double get_KEDAMM_theta(double percentage_done);

	int last_iteration_with_fitness_change = 0;
	int last_finess_change_counts = 0;

	bool theta_previously_increased = false;
};

#endif
