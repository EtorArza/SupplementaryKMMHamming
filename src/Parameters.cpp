#include "Parameters.h"
#include <assert.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "Tools.h"

#define N_OF_PARAMS 5



// argv[0] --> name of the binary executable
// argv[1] --> INSTANCE AND/OR FILE NAME
// argv[2] --> SEED
// argv[3] --> POPSIZE
// argv[4] --> EXPONENT


Parameters::Parameters(int argc, char *argv[])
{
  assert(N_OF_PARAMS == argc);


  // Parameters set at compile time
  this->FINAL_EXPECTATION = 0.25;
  this-> START_EXPECTATION_PERC = 0.5;
  // MAX_TIME = 600.0;


  // Parameters inputed from argv
  FILENAME = argv[1];
  
  SEED = atoi(argv[2]);
  srand(SEED);

  POPSIZE = atoi(argv[3]);

  EXPONENT = atof(argv[4]);


}

void Parameters::set_n(int n){
  this->n = n;
  this->START_EXPECTATION = (double) n * START_EXPECTATION_PERC;
  this->MAX_EVALS = 1000*n*n;
}


void Parameters::print_parameters(void)
{
  assert(parameters_have_been_set);
  std::cout << return_parameter_string();
}

std::string Parameters::return_parameter_string(void)
{
  assert(parameters_have_been_set);
  std::string result;
  result += "['";
  result += FILENAME;
  result += "', ";
  result += std::to_string(SEED);
  result += ", ";
  result += std::to_string(POPSIZE);
  result += ", ";
  result += std::to_string(EXPONENT);
  result += "]\n";

  return result;
}

