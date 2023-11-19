/*
 * EDAMM.cpp
 *
 *  Created on: Jun 6, 2018
 */
#include "KEDAMM.h"
#include <iostream>
#include "Tools.h"
#include "QAP.h"
#include <stdexcept>
#include <assert.h>
#include <math.h>
#include "Parameters.h"
#include <float.h>

#define USE_EXPO true

using namespace std;

//#define T_ADD 30

//#define FAST_INCREASE_EACH_ITERATION_COUNTS_AS 10

KEDAMM::KEDAMM(Parameters *param) : Population(param)
{
	this->N_OF_ITERATIONS = param->MAX_ITERATIONS;
}

void KEDAMM::learn()
{

    // print_vector("fitness_vector pre get_... calls",this->fitness_array,this->n);

	calculate_fitness();
	sort();

	// If best fitness was improved
	if (last_finess_change_counts < best_fit_change_count)
	{
		last_finess_change_counts = best_fit_change_count;
		if (theta > MIDDLE_THETA)
		{
			t_wait += param->T_ADD;
			t_wait = min(t_wait, param->T_WAIT_MAX);
		}
	}
	else
	{
		if (t_wait == 0)
		{
			if (theta < MIDDLE_THETA)
			{
				ker_iteration += param->FAST_INCREASE_EACH_ITERATION_COUNTS_AS;
			}
			else
			{
				ker_iteration++;
			}
			ker_percentage_done = (double)ker_iteration / (double)N_OF_ITERATIONS;
		}
		else
		{
			t_wait--;
		}
	}

	if (ker_percentage_done >= 1.0)
	{
		ker_finished = true;
		//cout << "-finished eda-" << endl;
	}

	theta = get_KEDAMM_theta(ker_percentage_done);
}

void KEDAMM::sample()
{
	for (int i = POPSIZE / 2; i < POPSIZE; i++)
	{
		sample_individual_at_position_i(i);
	}
}

void KEDAMM::sample_individual_at_position_i(int i)
{
	int r_index;
	r_index = random_integer_uniform(0, POPSIZE / 2 - 2);

	assert(i > r_index);	

	ham->sample_with_kernel(pop[r_index], theta, pop[i]);

	if (ham->fastfitness_can_be_used)
	{
		already_evaluated[i] = true;
		fitness_array[i] = fitness_array[r_index] + qap->update_fitness_on_swap(pop[r_index], pop[i]);
		ham->fastfitness_can_be_used = false;

		// print_vector("pop[r_index]", pop[r_index], n);
		// cout << "r: " << r_index << " , New Fitness " << fitness_array[i] << endl;
		// cout << "delta: " << qap->update_fitness_on_swap(pop[r_index], pop[i]) << endl;
		// print_vector("pop[i]", pop[i], n);
		// cout << "i: " << i << " , old fitness " << fitness_array[r_index] << endl;
	}
	else
	{
		already_evaluated[i] = false;
	}
}


double KEDAMM::get_KEDAMM_theta(double percentage_done)
{
	// plot f(x) = 0.25 + (exp((0.0005)*x) - 1) / (exp(0.0005) - 1) * (6 - 0.25)   from x=0 to x = 1	
	return param->START_MINIMUN_THETA + (exp(-1 * (param->EXPONENT + 0.0005) * percentage_done) - 1) / (exp(-1 * (param->EXPONENT + 0.0005)) - 1) * (FINAL_THETA - param->START_MINIMUN_THETA);
}


// double KEDAMM::get_KEDAMM_theta(double percentage_done)
// {
	
// }

KEDAMM::~KEDAMM()
{
	for (int i = 0; i < POPSIZE; i++)
	{
		delete[] ker_memory_pop[i];
	}
	delete[] ker_memory_pop;
	delete[] print_theta_consensus;
}
