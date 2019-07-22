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


using namespace std;

#define USE_ADAPTIVE_POPSIZE_IN_SAMPLE false

KEDAMM::KEDAMM(Parameters *param) : Population(param)
{

}

void KEDAMM::learn()
{

    // print_vector("fitness_vector pre get_... calls",this->fitness_array,this->n);

	calculate_fitness();
	sort();

	theta = get_KEDAMM_theta();
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
	if(USE_ADAPTIVE_POPSIZE_IN_SAMPLE){
		r_index = random_integer_uniform(0, tools_round((POPSIZE / 2) * ker_percentage_done) + 1);
	}else{
		r_index = random_integer_uniform(0, POPSIZE / 2 - 2);
	}

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


double KEDAMM::get_KEDAMM_theta()
{	
	
	// double progress = toc() / param->MAX_TIME;
	double progress = this->get_fitness_evals() / param->MAX_EVALS;

	progress = min(progress, 1.0);
	progress = max(progress, 0.0);
	
	// plot f(x) = 0 + (exp((EXPO + 0.0005)*x) - 1) / (exp(EXPO + 0.0005) - 1)   from x=0 to x = 1	
	progress = (exp(-1 * (param->EXPONENT + 0.0005) * progress) - 1) / (exp(-1 * (param->EXPONENT + 0.0005)) - 1);

	
	
	double target_expectation = (1 - progress) * (param->START_EXPECTATION - param->FINAL_EXPECTATION) + param->FINAL_EXPECTATION;
	
	assert(target_expectation > 0);
	assert(target_expectation < n);
	//print_variable("target_expc", target_expectation);
	theta = get_the_inverse_of_a_func_with_bisection(0.25, 6.00, target_expectation, [&](double x) { return ham->expectation(x); });
	return theta;
}




KEDAMM::~KEDAMM()
{
	for (int i = 0; i < POPSIZE; i++)
	{
		delete[] ker_memory_pop[i];
	}
	delete[] ker_memory_pop;
	delete[] print_theta_consensus;
}
