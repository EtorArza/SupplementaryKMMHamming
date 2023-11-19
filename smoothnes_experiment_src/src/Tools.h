#ifndef TOOLS_H
#define TOOLS_H
/*
 *  Tools.h
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <vector>

using std::istream;
using std::ostream;
using namespace std;
using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::stringstream;

/*
 * Returs the first position at which value appears in array. If it does not appear, then it returns -1;
 */
int Find(int *array, int size, int value);

/*
 * Calculates Kullback-Leibeler divergence between P and Q distributions.
 */
double KullbackLeibelerDivergence(double *P, double *Q, int size);

/*
 * Calculates Total Variation divergence between P and Q distributions.
 */
double TotalVariationDivergence(double *P, double *Q, int size);

/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int *permutation, int size, bool perm_contains_no_zeroes = false);

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int *permutation, int n);

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr);

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int *array, int length, string text);

/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double *array, int length, string text);

/*
 * Prints the given doubles array in the standard output.
 */
void PrintArray(double *array, int length, string text);


/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys(int *a, double *criteriaValues, int size);

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int *permutationA, int *permutationB, int size);

/*
 * Calculates the Kendall tau distance between 2 permutations.
 */
int Kendall(int *permutationA, int *permutationB, int size, int *m_aux);

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple executions.
 */
int Kendall(int *permutationA, int *permutationB, int size, int *m_aux, int *invertedB, int *composition, int *v);

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int *permutationA, int *permutationB, int size);

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int *permutationA, int *permutationB, int size, int *invertedB, int *composition, int *elemsToCycles, int *maxPosInCycle, int *freeCycle);
int FindNewCycle(int *freeCycle, int size);
int NextUnasignedElem(int *elemsToCycles, int size);
int CalculateDistance(int *sigma, int size);

/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int *sigma, int size);

/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int *s1, int *s2, int *res, int n);

/*
* Calculates V_j-s vector.
*/
void vVector(int *v, int *permutation, int n);

/*
 *  Optimized version by Leti of the V_j-s vector calculation.
 */
void vVector_Fast(int *v, int *permutation, int n, int *m_aux);

/*
 * Inverts a permutation.
 */
void Invert(int *permu, int n, int *inverted);

/*
 * This method moves the value in position i to the position j.
 */
void InsertAt(int *array, int i, int j, int n);

/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val);

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void swap_two_positions(int *array, int i, int j);

// not copied the ones below this


/*
 * Calculate the Hamming distance between two permutations
 */

int Hamming_distance(int* sigma1, int* sigma2, int len);

/*
 * Set timer to 0.
 */
void tic();

/*
 * Return time since last tic().
 */
double toc();

/*
 * Convert to string.
 */
template <class T>
string toString(const T &t, bool *ok = NULL)
{
    ostringstream stream;
    stream << t;
    if (ok != NULL)
        *ok = (stream.fail() == false);
    return stream.str();
}

/*
 * Convert to string.
 * https://stackoverflow.com/questions/3909272/sorting-two-corresponding-arrays
 * //sort 2 arrays simultaneously
 */
template <class A, class B>
void QuickSort2Desc(A a[], B b[], int l, int r)
{
    int i = l;
    int j = r;
    A v = a[(l + r) / 2];
    do
    {
        while (a[i] > v)
            i++;
        while (v > a[j])
            j--;
        if (i <= j)
        {
            std::swap(a[i], a[j]);
            std::swap(b[i], b[j]);
            i++;
            j--;
        };
    } while (i <= j);
    if (l < j)
        QuickSort2Desc(a, b, l, j);
    if (i < r)
        QuickSort2Desc(a, b, i, r);
}





//sort 3 arrays simultaneously
template <class A, class B, class C>
void QuickSort3Desc(A a[], B b[], C c[], int l, int r)
{
	int i = l;
	int j = r;
	A v = a[(l + r) / 2];
	do
	{
		while (a[i] > v)
			i++;
		while (v > a[j])
			j--;
		if (i <= j)
		{
			std::swap(a[i], a[j]);
			std::swap(b[i], b[j]);
			std::swap(c[i], c[j]);

			i++;
			j--;
		};
	} while (i <= j);
	if (l < j)
		QuickSort3Desc(a, b, c, l, j);
	if (i < r)
		QuickSort3Desc(a, b, c, i, r);
}


/*
 * Return wether two vectors are equal or not.
 */
template <class T>
bool compare_vectors(T *vec1, T *vec2, int len)
{
    for (int i = 0; i < len; i++)
    {
        if (vec1[i] != vec2[i])
        {
            return false;
        }
    }
    return true;
}

/*
Function to find all the repeated rows on a matrix
Writes in  bool *is_ith_position_repeated (true --> vector is a repetition, false--> vector is not a repetition)
*/
template <class T>
void which_indexes_correspond_to_repeated_vectors(T **vec_array, int vec_len, int n_of_vecs, bool *is_ith_position_repeated, bool is_known_last_repeated_indexes)
{
    if (n_of_vecs == 1)
    {
        is_ith_position_repeated[0] = false;
        return;
    }
    else if (n_of_vecs == 2)
    {
        is_ith_position_repeated[0] = false;
        is_ith_position_repeated[1] = compare_vectors(vec_array[0], vec_array[1], vec_len);
        return;
    }

    is_ith_position_repeated[0] = false;
    is_ith_position_repeated[1] = compare_vectors(vec_array[0], vec_array[1], vec_len);

    for (int i = 2; i < n_of_vecs; i++)
    {
        if (is_known_last_repeated_indexes && not is_ith_position_repeated[i])
        {
            continue;
        }
        for (int j = i - 1; j >= 0; j--)
        {
            is_ith_position_repeated[i] = false;
            if (compare_vectors(vec_array[i], vec_array[j], vec_len))
            {
                is_ith_position_repeated[i] = true;
                break;
            }
        }
    }
}

/*
Shuffle vector given its length.
*/
void shuffle_vector(int *vec, int len);

/*
Get random integer on interval [min, max - 1], faster but slightly biased
*/
int random_integer_fast(int min, int max);

// Get random integer on interval [min, max - 1]
// https://ericlippert.com/2013/12/16/how-much-bias-is-introduced-by-the-remainder-technique/
int random_integer_uniform(int min, int max = 0);

// chooses a random integer from {0,1,2, range_max - 1}, at uniform (may be a little slower)
int random_range_integer_uniform(int range_max);

// return a random uniform float on the interval [0,1]
double random_0_1_float();

// apply sigmoid function
double sigmoid(double x);

// Choose an index given the probabilities
int choose_index_given_probabilities(double *probabilities_array, int max_index);


// Choose an index given positive weights
int choose_index_given_wheights(double *weights_array, int max_index);

// Sample from a bernouilli distribution.
bool coin_toss(double p_of_true);

// round a double into the nearest integer
int tools_round(double x);

// compute the average value of the elements on the array
template <class T>
double Average(T *array, int len)
{

    double sum = 0;

    for (int i = 0; i < len; i++)
    {
        sum += array[i];
    }

    return (double)sum / len;
}

// compute the variance of the elements on the array
template <class T>
double Variance(T *array, int len)
{

    double mean = Average(array, len);

    double var = 0;
    for (int i = 0; i < len; i++)
    {
        var += (array[i] - mean) * (array[i] - mean);
    }

    return (double)var / len;
}

// Normalize a vector so that the sum of all the elements on it is 1
template <class T>
double normalize_vector(T *array, int len)
{
    int sum = 0;
    for (int i = 0; i < len; i++)
    {
        sum += array[i];
    }
    for (int i = 0; i < len; i++)
    {
        array[i] = array[i] / sum;
    }
}

void select_indexes_to_insert_towards(int* genome_to_be_changed, int* reference_genome, int* i, int* j);



template <class T>
void PrintMatrix(T **M, int m, int n)
{

	cout << "\n";
	for (int i = 0; i < m; i++)
	{
		cout << "| i = " << i << " ( ";
		for (int j = 0; j < n; j++)
		{
			cout << M[i][j] << " ";
		}
		cout << ")\n";
	}
}


template <class T>
void copy_vector(T *v1,T *v2, int n){
    for (int i = 0; i < n; i++)
    {
        v2[i] = v1[i];
    }
    
}




#endif /* TOOLS_H */
