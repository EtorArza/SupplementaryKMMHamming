#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "Parameters.h"

#include "KEDAMM.h"
#include "Population.h"
#include "Tools.h"
#include "SimpleLog.h"
#include "Parameters.h"
using namespace std;

Population *eda;
Parameters *param;

#define HOW_MANY_MAX_PRINT_STATUS 100000

int main(int argc, char *argv[])
{

  int best_fitness_eda = 0;
  int n_evals = 0;

  param = new Parameters(argc, argv);
  eda = new KEDAMM(param);

  //param->print_parameters();

  int purge_every = 100;
  int i = 0;
  int evals = 0;

  tic();

  if (param->DEBUG)
  {
    param->print_parameters();
  }
  

  // while (toc() < param->MAX_TIME * param->PERC_TIME_EDA)
  while (n_evals < (int) param->MAX_EVALS)
  {

    i++;
    if (param->DEBUG)
    {
      print_variable("time: ", toc(), false);
      print_variable("| evals: ", eda->get_fitness_evals(), false);
      print_variable("| fitness: ", eda->get_best_fitness(), false);
      print_variable(" | theta:", eda->theta);
    }
    eda->learn();


    n_evals = eda->get_fitness_evals();


    if (i % purge_every == 0) // purging hast to be done only after the learning phase, else, some individuals may have the wrong fitness value
    {
      // cout << i << ": " << eda->ker_percentage_done << " | ";
      eda->purge_repeated();
    }

    eda->sample();

    // if (i % print_every == 0)
    // {
    //   cout << i << ", ";
    //   cout << eda->ker_percentage_done << ": ";
    //   cout << eda->get_best_fitness() << endl;
    //   best_fitness_eda = eda->get_best_fitness();
    // }
  }

  eda->learn();

  //print out the results [n_evals_eda, score_eda, n_evals_local, score_local]
  //cout << "[";
  //cout << n_evals_eda << ", ";
  //cout << best_fitness_eda << ", ";
  //cout << eda->get_fitness_evals() << ", ";
  cout << param->FILENAME << "|";
  cout << param->SEED << "|";
  eda->evaluate_individual_on_position_i(0);
  cout << eda->fitness_array[0];

  // PrintPythonArray(eda->best_permutation, eda->n);

  // // Write result summary on results.txt
  // SimpleLog result_writer("results.txt", false);
  // result_writer.write(best_fitness_eda, false);
  // result_writer.write(" , ", false);
  // result_writer.write(eda->best_fitness, false);
  // result_writer.write(" , ", false);
  // result_writer.write(param->return_parameter_string(), false);
  exit(0);
}
