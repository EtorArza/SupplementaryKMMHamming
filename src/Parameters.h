
#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include <string>




class Parameters
{
public:

  Parameters(int argc, char *argv[]);

  bool DEBUG = false;




  int POPSIZE;
  int SEED;
  int n;
  double MAX_TIME;
  double FINAL_EXPECTATION;
  double START_EXPECTATION;
  double PERC_TIME_EDA;
  double MAX_EVALS;
  double EXPONENT;


  std::string FILENAME;




  void print_parameters();
  std::string return_parameter_string();
  void set_n(int n);

  


private:
  double START_EXPECTATION_PERC;
  bool parameters_have_been_set = true;
};

#endif
