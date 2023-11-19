
#include "Parameters.h"
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "Tools.h"

#define N_OF_PARAMS 3


// argv[0] --> name of the binary executable
// argv[1] --> POPSIZE
// argv[2] --> MAX_ITERATIONS
// argv[3] --> SEED
// argv[4] --> START_MINIMUN_THETA
// argv[5] --> TARGET_MID_EXPECTATION_PERCENTAGE
// argv[6] --> FINAL_EXPECTATION
// argv[7] --> T_ADD
// argv[8] --> FAST_INCREASE_EACH_ITERATION_COUNTS_AS
// argv[9] --> EXPONENT
// argv[10] --> INSTANCE AND FILE NAME

Parameters::Parameters(int argc, char *argv[])
{
  assert(N_OF_PARAMS == argc);

  POPSIZE = 20000;

  MAX_ITERATIONS = 1250;

  SEED = atoi(argv[1]);
  srand(SEED);

  START_MINIMUN_THETA = 0.25;

  EXPONENT = 0.01;

  T_ADD = 20;

  FAST_INCREASE_EACH_ITERATION_COUNTS_AS = 10;

  FILENAME = std::string(argv[2]);

  parameters_have_been_set = true;

  T_WAIT_MAX = WAIT_T_ADD_MULTIPLYER * T_ADD;  
}

void Parameters::set_n(int n){
  this->n = n;
  int a;
  int b;
  #define N_OF_KNOWN_SIZES 4
  int first_bigger = 4;
  for(int i = 0; i < N_OF_KNOWN_SIZES; i++)
  { 
    if (known_problem_sizes[i] > n) {
      first_bigger = i;
      break;
    }
  }
  

  if (first_bigger == 0) {
    TARGET_MID_EXPECTATION_PERCENTAGE = known_mid[0];
    FINAL_EXPECTATION = known_final[0];
  }

   

  else if(first_bigger == N_OF_KNOWN_SIZES) {
    TARGET_MID_EXPECTATION_PERCENTAGE = known_mid[N_OF_KNOWN_SIZES - 1];
    FINAL_EXPECTATION = known_final[N_OF_KNOWN_SIZES - 1];
  }
  else {

        TARGET_MID_EXPECTATION_PERCENTAGE = \
    (
      known_mid[first_bigger] * ( n - known_problem_sizes[first_bigger - 1] ) +
      known_mid[first_bigger - 1] * ( known_problem_sizes[first_bigger] - n ) 
    ) / (known_problem_sizes[first_bigger] - known_problem_sizes[first_bigger - 1]);

      FINAL_EXPECTATION = \
    (
      known_final[first_bigger] * ( n - known_problem_sizes[first_bigger - 1] ) +
      known_final[first_bigger - 1] * ( known_problem_sizes[first_bigger] - n ) 
    ) / (known_problem_sizes[first_bigger] - known_problem_sizes[first_bigger - 1]);

  }
    



}



void Parameters::print_parameters(void)
{
  assert(parameters_have_been_set);
  std::cout << return_parameter_string();
  //std::cout << POPSIZE << ", ";
  //std::cout << MAX_ITERATIONS << ", ";
  //std::cout << SEED << ", ";
  //std::cout << START_MINIMUN_THETA << ", ";
  //std::cout << TARGET_MID_EXPECTATION_PERCENTAGE << ", ";
  //std::cout << FINAL_EXPECTATION << ", ";
  //std::cout << T_ADD << ", ";
  //std::cout << FAST_INCREASE_EACH_ITERATION_COUNTS_AS << "] " << std::endl;
}

std::string Parameters::return_parameter_string(void)
{
  assert(parameters_have_been_set);
  std::string result;
  result += "[";
  result += std::to_string(POPSIZE);
  result += ", ";
  result += std::to_string(MAX_ITERATIONS);
  result += ", ";
  result += std::to_string(SEED);
  result += ", ";
  result += std::to_string(START_MINIMUN_THETA);
  result += ", ";
  result += std::to_string(TARGET_MID_EXPECTATION_PERCENTAGE);
  result += ", ";
  result += std::to_string(FINAL_EXPECTATION);
  result += ", ";
  result += std::to_string(T_ADD);
  result += ", ";
  result += std::to_string(FAST_INCREASE_EACH_ITERATION_COUNTS_AS);
  result += ", '";
  result += FILENAME; // FILENAME is already a string.
  result += "']\n";

  return result;
}