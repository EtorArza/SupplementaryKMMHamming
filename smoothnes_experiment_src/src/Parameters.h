
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <string>


#define WAIT_T_ADD_MULTIPLYER 3


class Parameters
{
public:
  Parameters(int argc, char *argv[]);

  // argv[0] --> name of the binary executable
  // argv[1] --> POPSIZE
  // argv[2] --> MAX_ITERATIONS
  // argv[3] --> SEED
  // argv[4] --> START_MINIMUN_THETA
  // argv[5] --> TARGET_MID_EXPECTATION_PERCENTAGE
  // argv[6] --> FINAL_EXPECTATION
  // argv[7] --> T_ADD
  // argv[8] --> FAST_INCREASE_EACH_ITERATION_COUNTS_AS

  int POPSIZE;
  int MAX_ITERATIONS;
  int SEED;
  int n;

  double START_MINIMUN_THETA;
  double TARGET_MID_EXPECTATION_PERCENTAGE;
  double FINAL_EXPECTATION;
  double EXPONENT;

  int T_ADD;
  int FAST_INCREASE_EACH_ITERATION_COUNTS_AS;

  int T_WAIT_MAX;

  std::string FILENAME;


  bool DEBUG = false;

  int known_problem_sizes[4] = {27, 75, 125, 175}; 
  double known_mid[4] = {0.75, 0.60, 0.56, 0.24};
  double known_final[4] = {0.45, 8.00, 1.55, 0.67};

  void print_parameters();
  std::string return_parameter_string();
  void set_n(int n);

  


private:
  bool parameters_have_been_set = false;
};

#endif
