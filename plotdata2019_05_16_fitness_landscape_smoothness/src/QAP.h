/*
 * QAP.h
 * This class is used to load a QAP problem, and evaluate permutations on
 * the given problem.
 *  Created on: Jun 27, 2018
 *      Author: paran
 */

#ifndef QAP_H_
#define QAP_H_

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>

using std::ifstream;
using std::istream;
using std::ofstream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;

class QAP
{

  public:
	QAP();
	~QAP();

	//read problem from txt file.
	void read(string filename);

	//read problem form cin.
	void read(void);

	// evaluate the fitness of a permutation. It is a negative number.
	int evaluate(int *sigma);

	// calc fitness based on swap
	int update_fitness_on_swap(int *sigma_1, int *sigma_2);

	// returns the dimension of the matrix
	int get_dim()
	{
		return n;
	}

	int n;

	//protected:
	// Inherited properties:
	//int permutation_size;

	int *m_aux;
	int **m_distance_matrix;
	int **m_flow_matrix;
};

#endif /* QAP_H_ */
