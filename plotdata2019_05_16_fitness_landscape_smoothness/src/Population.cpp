/*
 * Population.cpp
 *
 *  Created on: Jun 27, 2018
 *      Author: paran
 */
#include "Population.h"
#include "Parameters.h"
#include "Tools.h"
#include <float.h>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <cstdlib>

using namespace std;
using std::string;
using std::stringstream;


Population::Population(Parameters *param)
{
  this->param = param;
  qap = new QAP();
  qap->read();
  initialize_Population();
}

void Population::initialize_Population()
{
  n = qap->get_dim();
  assert(n > 0);

  param->set_n(n);

  best_fitness = -INT_MAX;
  fitness_evals = 0;
  best_fit_change_count = 0;

  //POPSIZE = tools_round((double)n * (double)n * pop_mult_factor);
  POPSIZE = param->POPSIZE;

  ham = new HammingFunctions(n, qap);
  fitness_array = new int[POPSIZE];
  best_permutation = new int[n];
  best_permutation_one_position_swapped_local_search = new int[n];
  pop = new int *[POPSIZE];
  for (int i = 0; i < POPSIZE; ++i)
  {
    pop[i] = new int[n];
  }

  // initialize random uniform population
  ham->generate_random_sample(POPSIZE, pop);

  consensus = new int[n];
  theta_array = new double[n];
  theta = 0.0;

  repeated_indexes = new bool[POPSIZE];

  // Set max theta based on expectation
  double a = 0.25;
  double b = 10;
  double fab;
  double middle_ab;

  for (int i = 0; i < 1000; i++)
  {
    middle_ab = (a + b) / 2.0;
    fab = ham->expectation(middle_ab, n);
    if (fab > param->FINAL_EXPECTATION)
    {
      a = middle_ab;
    }
    else
    {
      b = middle_ab;
    }
  }



  FINAL_THETA = middle_ab;

  // Set FAST_INCREASE_UNTIL to get perc_exp * n expectation on middle point
  a = 0.25;
  b = 10;


  for (int i = 0; i < 1000; i++)
  {
    middle_ab = (a + b) / 2.0;
    fab = ham->expectation(middle_ab, n);
    if (fab > (double)n * param->TARGET_MID_EXPECTATION_PERCENTAGE)
    {
      a = middle_ab;
    }
    else
    {
      b = middle_ab;
    }
  }


  MIDDLE_THETA = middle_ab;
  //cout << "FINAL_THETA: " << FINAL_THETA << endl;
  //cout << "MIDDLE_THETA:" << MIDDLE_THETA << endl;

  already_evaluated = new bool[POPSIZE];
  for (int i = 0; i < POPSIZE; ++i)
  {
    already_evaluated[i] = false;
    evaluate_individual_on_position_i(i);
  }
  
  random_i_indexes = new int[n];
  random_j_indexes = new int[n];

  
  for(int k = 0; k < n; k++)
  {
    random_i_indexes[k] = k;
    random_j_indexes[k] = k;
  }
  
  shuffle_vector(random_i_indexes, n);
  shuffle_vector(random_j_indexes, n);


  assert(n < POPSIZE);
}


void Population::evaluate_individual_on_position_i(int i)
{
  // evaluate permutation
  fitness_array[i] = qap->evaluate(pop[i]);
  fitness_evals++;

  // check if best solution
  if (fitness_array[i] > best_fitness)
  {
    best_fitness = fitness_array[i];
    copy_vector(pop[i], best_permutation, n);
    best_fit_change_count++;
  };
}

void Population::calculate_fitness()
{

  
  
  if (param -> DEBUG)
  {
    for (int i = 0; i < POPSIZE / 2; ++i)
    {
      assert(fitness_array[i] == qap->evaluate(pop[i]));
    }
  }


  for (int i = POPSIZE / 2; i < POPSIZE; ++i)
  {
    if (!already_evaluated[i])
    {
      evaluate_individual_on_position_i(i);
    }
    else
    {

      if (param->DEBUG)
      {
        int last_fitness = fitness_array[i];
        // print_vector("pop[i]", pop[i], n);
        // print_variable("i", i);
        evaluate_individual_on_position_i(i);
        // cout << last_fitness << " - " << fitness_array[i] << endl;
        assert(last_fitness == fitness_array[i]);
      }
    }
  }
}

// sort the population and the fitness_array
void Population::sort() { QuickSort3Desc(fitness_array, pop, already_evaluated, 0, POPSIZE - 1); }

void Population::get_best_permutaion(int *sigma)
{
  for (int i = 0; i < n; ++i)
  {
    sigma[i] = best_permutation[i];
  }
}

void Population::local_search()
{
  // need goto in order to be able to break nested for loops
  //cout << "-begin local search-" << endl;

  //int log_n = tools_round(log((double) n));
  

  int HOW_MANY_LOCAL_SEARCHES = min(POPSIZE / 4, 1000);
  int fitness_evals_until_local_search = get_fitness_evals(); 


  for(int k = 0; k < HOW_MANY_LOCAL_SEARCHES; k++)
  {

    int new_fitness_delta = 0;
    int best_delta = 0;
    int best_i = 0;
    int best_j = 0;
    bool fitness_was_improved = true;


    copy_vector(pop[k], best_permutation_one_position_swapped_local_search, n);
 


    while (fitness_was_improved)
    {
      fitness_was_improved = false;
      new_fitness_delta = 0;
      best_delta = 0;
      shuffle_vector(random_i_indexes, n);
      shuffle_vector(random_j_indexes, n);

      int i;
      int j;

      for (int i_pos = 0; i_pos < n; i_pos++)
      {
        for (int j_pos = 0; j_pos < n; j_pos++)
        {

          i = random_i_indexes[i_pos];
          j = random_j_indexes[j_pos];
          
          if (i != j)
          {
            swap_two_positions(best_permutation_one_position_swapped_local_search, i, j);
            new_fitness_delta = qap->update_fitness_on_swap(pop[k], best_permutation_one_position_swapped_local_search);
            fitness_evals++;
            //cout << i << ", " << j << " -- " << new_fitness_delta << endl;

            if (new_fitness_delta > best_delta)
            {
              best_i = i;
              best_j = j;
              best_delta = new_fitness_delta;
              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              // the following two lines make the ls first improvement (uncommented) or best improvement (commented)//
              swap_two_positions(best_permutation_one_position_swapped_local_search, i, j);                         
              goto END_LOOP;                                                                                        
              //////////////////////////////////////////////////////////////////////////////////////////////
            }

            swap_two_positions(best_permutation_one_position_swapped_local_search, i, j);
          }
        }
      }
      END_LOOP:

      if (best_delta > 0)
      {
        fitness_was_improved = true;
        swap_two_positions(best_permutation_one_position_swapped_local_search, best_i, best_j);
        swap_two_positions(pop[k], best_i, best_j);
        fitness_array[k] += best_delta;
      }
      else
      {
        fitness_was_improved = false;
      }
    }

    // print_variable("Sol_number: -> ", k, false);
    // cout << " ";
    // print_variable("best_fitness", fitness_array[k]);
    // If fitness evaluations on local exceeds fitness evaluations on EDA, terminate.
    if (get_fitness_evals() > fitness_evals_until_local_search * 2) {
      break;
    }
    
  }
  //cout << "-end local search-" << endl;
  
  sort();

  copy_vector(pop[0], best_permutation, n);

  best_fitness = fitness_array[0]; 

}

void Population::purge_repeated()
{
  int rep_counter = 0;
  int rep_fitness_counter = 0;
  for (int i = 0; i < POPSIZE / 2 - 1; i++)
  {
    if (fitness_array[i] == fitness_array[i + 1])
    {
      if (param->DEBUG)
      {
        rep_fitness_counter++;
      }
      if (compare_vectors(pop[i], pop[i + 1], n))
      {
        rep_counter++;
        shuffle_vector(pop[i], n);
        evaluate_individual_on_position_i(i);
      }
    }
  }
  for(int i = POPSIZE/2 - 2 ; i < POPSIZE/2 + rep_counter + 2; i++)
  {
    evaluate_individual_on_position_i(i);
  }
  if (param->DEBUG) {
  cout << "Repeated individuals: " << rep_counter << "rep_fitness_counter: " << rep_fitness_counter << endl;
  }
  this->sort();
}

Population::~Population()
{
  delete ham;
  //delete qap;
  for (int i = 0; i < POPSIZE; ++i)
  {
    delete[] pop[i];
  }
  delete[] pop;
  delete[] fitness_array;
  delete[] best_permutation;
  delete[] consensus;
  delete[] theta_array;
  delete[] repeated_indexes;
  delete[] already_evaluated;
  delete[] best_permutation_one_position_swapped_local_search;
  delete[] random_i_indexes;
  delete[] random_j_indexes;
}
