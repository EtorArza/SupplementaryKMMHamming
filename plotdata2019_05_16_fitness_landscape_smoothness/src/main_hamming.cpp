#include <assert.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <fcntl.h>
#include "Parameters.h"

#include "KEDAMM.h"
#include "Population.h"
#include "Tools.h"
#include "SimpleLog.h"
#include "Parameters.h"

#include <stdio.h>

#include "Cayley.h"
#include "Ulam.h"
#include "Kendall.h"
#include "Generic.h"
#include "Exponential_model.h"

using namespace std;

Population *eda;
Parameters *param;
#define n 45
#define AVERAGE false

#define HOW_MANY_WAIT_TO_PRINT 2
#define EL_PER_LVL 50

#define CAYLEY_DISTANCE 0
#define KENDALL_DISTANCE 1
#define HAMMING_DISTANCE 2
#define ULAM_DISTANCE 3

void print_permutation(int *sigma)
{

  for (int i = 0; i < n; i++)
  {
    cout << sigma[i] << " ";
  }
  cout << endl;
}

void summ_to_all(int *permu, int diff)
{
  for (int i = 0; i < n; i++)
  {
    permu[i] += diff;
  }
}

void generate_at_hamming_distance_d(int *sigma_0, int *sigma, int d)
{
  assert(d != 1);

  isPermutation(sigma_0, n, true);

  int *h;
  int *Id;
  h = new int[n];
  Id = new int[n];

  for (int i = 0; i < n; i++)
  {
    Id[i] = i + 1;
  }

  for (int i = 0; i < d; i++)
  {
    h[i] = 1;
  }

  for (int i = d; i < n; i++)
  {
    h[i] = 0;
  }

  //  Yates shuffle
  for (int i = n - 1; i > 0; --i)
  {
    int r_index = random_range_integer_uniform(i + 1);
    int temp = h[i];
    h[i] = h[r_index];
    h[r_index] = temp;
  }

  eda->ham->dist_decomp_vector2perm(h, Id);
  eda->ham->compose(Id, sigma_0, sigma);

  assert(isPermutation(sigma, n, true));

  delete[] h;
  delete[] Id;
}

// ulam copied to ekhiñe, therefore, unlike other distances (except hamming) perms are considerd as (1 2 3 4 ... n)
void generate_at_cku_distance_d(int *sigma_0, int *sigma, int d, Exponential_model *exp_mod, Generic *gen)
{

  int *central_change = new int[n];

  for (int i = 0; i < n; i++)
  {
    central_change[i] = i + 1;
  }

  exp_mod->distances_sampling(d, central_change);
  gen->compose(n, sigma_0, central_change, sigma);

  delete[] central_change;
}

int evaluate_at_distance_d(int *sigma_0, int *sigma, int distance_type, int d, Exponential_model *exp_mod, Generic *gen)
{
  copy_vector(sigma_0, sigma, n);
  switch (distance_type)
  {
  case HAMMING_DISTANCE:
    generate_at_hamming_distance_d(sigma_0, sigma, d);
    assert(isPermutation(sigma_0, n, true));
    assert(isPermutation(sigma, n, true));
    break;

  case CAYLEY_DISTANCE:
    generate_at_cku_distance_d(sigma_0, sigma, d, exp_mod, gen);
    assert(isPermutation(sigma_0, n, true));
    assert(isPermutation(sigma, n, true));
    break;

  case KENDALL_DISTANCE:
    generate_at_cku_distance_d(sigma_0, sigma, d, exp_mod, gen);
    assert(isPermutation(sigma_0, n, true));
    assert(isPermutation(sigma, n, true));
    break;

  case ULAM_DISTANCE:
    generate_at_cku_distance_d(sigma_0, sigma, d, exp_mod, gen);
    assert(isPermutation(sigma_0, n, true));
    assert(isPermutation(sigma, n, true));

    break;

  default:
    break;
  }

  int score = -eda->qap->evaluate(sigma);
  copy_vector(sigma_0, sigma, n);
  return score;
}

void get_random_local_optima(int *sigma, QAP *qap)
{

  int *random_i_indexes;
  int *random_j_indexes;
  int *new_sigma;

  random_i_indexes = new int[n];
  random_j_indexes = new int[n];
  new_sigma = new int[n];

  for (int i = 0; i < n; i++)
  {
    random_i_indexes[i] = i;
    random_j_indexes[i] = i;
    sigma[i] = i + 1;
  }
  shuffle_vector(sigma, n);
  int new_fitness_delta = 0;
  int best_delta = 0;
  int best_i = 0;
  int best_j = 0;
  bool fitness_was_improved = true;

  copy_vector(sigma, new_sigma, n);

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
          swap_two_positions(new_sigma, i, j);
          new_fitness_delta = qap->update_fitness_on_swap(sigma, new_sigma);
          //cout << i << ", " << j << " -- " << new_fitness_delta << endl;

          if (new_fitness_delta > best_delta)
          {
            best_i = i;
            best_j = j;
            best_delta = new_fitness_delta;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////
            // the following two lines make the ls first improvement (uncommented) or best improvement (commented)//
            swap_two_positions(new_sigma, i, j);
            goto END_LOOP;
            //////////////////////////////////////////////////////////////////////////////////////////////
          }

          swap_two_positions(new_sigma, i, j);
        }
      }
    }
  END_LOOP:

    if (best_delta > 0)
    {
      fitness_was_improved = true;
      swap_two_positions(new_sigma, best_i, best_j);
      swap_two_positions(sigma, best_i, best_j);
    }
    else
    {
      fitness_was_improved = false;
    }
  }
  assert(isPermutation(sigma, n, true));
  delete[] new_sigma;
  delete[] random_i_indexes;
  delete[] random_j_indexes;
}

int main(int argc, char *argv[])
{

  // {
  //   // programmatically redirect stdio
  //   const char *stdin_filename = INSTANCE_LOCATION; //tai125e01.dat.dat"; //nug6.dat.dat ";
  //   assert(dup2(open(stdin_filename, O_RDONLY), 0) != -1);
  //   asm("int3"); // optional breakpoint -- kills program when not debugging
  // }

  // setup
  Generic gen;
  Exponential_model *ulam_exp_mod = gen.new_instance(ULAM_DISTANCE, n);
  Exponential_model *cayley_exp_mod = gen.new_instance(CAYLEY_DISTANCE, n);
  Exponential_model *kendall_exp_mod = gen.new_instance(KENDALL_DISTANCE, n);

  //gen.seed();
  //loop_menu();

  int d = 2;
  int *sigma = new int[n];
  int *sigma_0 = new int[n];

  for (int i = 0; i < n; i++)
  {
    sigma[i] = i + 1;
    sigma_0[i] = i + 1;
  }

  shuffle_vector(sigma_0, n);

  generate_at_cku_distance_d(sigma_0, sigma, d, cayley_exp_mod, &gen);
  print_permutation(sigma_0);
  print_permutation(sigma);
  cout << "---------" << endl;
  generate_at_cku_distance_d(sigma_0, sigma, d, kendall_exp_mod, &gen);
  print_permutation(sigma_0);
  print_permutation(sigma);
  cout << "---------" << endl;
  generate_at_cku_distance_d(sigma_0, sigma, d, ulam_exp_mod, &gen);
  print_permutation(sigma_0);
  print_permutation(sigma);

  param = new Parameters(argc, argv);
  eda = new KEDAMM(param);

  assert(isPermutation(sigma_0, n, true));

  shuffle_vector(sigma_0, n);
  shuffle_vector(sigma_0, n);

  get_random_local_optima(sigma_0, eda->qap);

  assert(isPermutation(sigma_0, n, true));

  copy_vector(sigma_0, sigma, n);

  double *res;
  res = new double[EL_PER_LVL];
  for (int i = 0; i < 100; i++)
  {
    shuffle_vector(sigma_0, n);
    get_random_local_optima(sigma_0, eda->qap);
    assert(isPermutation(sigma_0, n, true));

    int f_sigma = -eda->qap->evaluate(sigma);

    for (int d = 1; d < 15; d++)
    {
      if (d != 1)
      {
        for (int i = 0; i < EL_PER_LVL; i++)
          res[i] = (double)abs(abs(f_sigma) - abs(evaluate_at_distance_d(sigma_0, sigma, HAMMING_DISTANCE, d, cayley_exp_mod, &gen))) / (double)abs(f_sigma);
        if (AVERAGE)
        {
          cout << Average(res, EL_PER_LVL) << ", ";
        }
        else
        {
          cout << Variance(res, EL_PER_LVL) << ", ";
        }
      }
      else
      {
        cout << 0 << ", ";
      }
    }
    cout << "$" << endl;
    for (int d = 1; d < 15; d++)
    {
      for (int i = 0; i < EL_PER_LVL; i++)
      {
        res[i] = (double)abs(abs(f_sigma) - abs(evaluate_at_distance_d(sigma_0, sigma, CAYLEY_DISTANCE, d, cayley_exp_mod, &gen))) / (double)abs(f_sigma);
      }
      if (AVERAGE)
      {
        cout << Average(res, EL_PER_LVL) << ", ";
      }
      else
      {
        cout << Variance(res, EL_PER_LVL) << ", ";
      }
    }
    cout << "$" << endl;
    for (int d = 1; d < 15; d++)
    {
      for (int i = 0; i < EL_PER_LVL; i++)
      {
        res[i] = (double)abs(abs(f_sigma) - abs(evaluate_at_distance_d(sigma_0, sigma, KENDALL_DISTANCE, d, kendall_exp_mod, &gen))) / (double)abs(f_sigma);
      }
      if (AVERAGE)
      {
        cout << Average(res, EL_PER_LVL) << ", ";
      }
      else
      {
        cout << Variance(res, EL_PER_LVL) << ", ";
      }
    }
    cout << "$" << endl;
    for (int d = 1; d < 15; d++)
    {
      for (int i = 0; i < EL_PER_LVL; i++)
      {
        res[i] = (double)abs(abs(f_sigma) - abs(evaluate_at_distance_d(sigma_0, sigma, ULAM_DISTANCE, d, ulam_exp_mod, &gen))) / (double)abs(f_sigma);
      }
      if (AVERAGE)
      {
        cout << Average(res, EL_PER_LVL) << ", ";
      }
      else
      {
        cout << Variance(res, EL_PER_LVL) << ", ";
      }
    }
    cout << "$" << endl;
  }

  cout << "-----"
       << "end" << endl;

  //print out the results [n_evals_eda, score_eda, n_evals_local, score_local]

  // // Write result summary on results.txt
  // SimpleLog result_writer("results.txt", false);
  // result_writer.write(best_fitness_eda, false);
  // result_writer.write(" , ", false);
  // result_writer.write(eda->best_fitness, false);
  // result_writer.write(" , ", false);
  // result_writer.write(param->return_parameter_string(), false);
  exit(0);
}
