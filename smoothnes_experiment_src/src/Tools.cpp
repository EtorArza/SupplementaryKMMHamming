/*
 *  Tools.cpp
 *  RankingEDAsCEC
 *
 *  Created by Josu Ceberio Uribe on 11/21/11.
 *  Copyright 2011 University of the Basque Country. All rights reserved.
 *
 */

#include "Tools.h"
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h> /* srand, rand */
#include <assert.h>
#include "Parameters.h"

#define MAX_INTEGER 100000000

/*
 * Returs the first position at which value appears in array. If it does not appear, then it returns -1;
 */
int Find(int *array, int size, int value)
{
    int i = 0;
    while (i < size)
    {
        if (array[i] == value)
            return i;
        i++;
    }
    return -1;
}

/*
 * Calculates Kullback-Leibeler divergence between P and Q distributions.
 */
double KullbackLeibelerDivergence(double *P, double *Q, int SearchSpaceSize)
{
    double divergence = 0;
    double auxi = 0;
    for (int i = 0; i < SearchSpaceSize; i++)
    {
        auxi = P[i] * (log(P[i] / Q[i]));
        divergence = divergence + auxi;
    }
    return divergence;
}

/*
 * Calculates Total Variation divergence between P and Q distributions.
 */
double TotalVariationDivergence(double *P, double *Q, int size)
{
    double divergence = 0;
    for (int i = 0; i < size; i++)
    {
        divergence += fabs(P[i] - Q[i]);
    }
    return (0.5 * divergence);
}

/*
 * It determines if the given int sequecen if it is indeed a permutation or not.
 */
bool isPermutation(int *permutation, int size, bool perm_contains_no_zeroes)
{
    if (!perm_contains_no_zeroes)
    {
        int flags[size];
        //int * flags=new int[size];
        for (int i = 0; i < size; i++)
        {
            if (permutation[i] >= size)
            {
                return false;
            }
            flags[i] = 1;
        }

        for (int i = 0; i < size; i++)
        {
            int value = permutation[i];
            flags[value] = 0;
        }

        int result, sum = 0;
        for (int i = 0; i < size; i++)
            sum += flags[i];
        if (sum == 0)
            result = true;
        else
            result = false;
        //delete [] flags;
        return result;
    }
    else
    {
        for (int i = 0; i < size; i++)
        {
            permutation[i]--;
        }
        bool result = isPermutation(permutation, size, false);
        for (int i = 0; i < size; i++)
        {
            permutation[i]++;
        }
        return result;
    }
}

/*
 * Generates a random permutation of size 'n' in the given array.
 */
void GenerateRandomPermutation(int *permutation, int n)
{
    for (int i = 0; i < n; ++i)
    {
        int j = rand() % (i + 1);
        permutation[i] = permutation[j];
        permutation[j] = i;
    }
}

/*
 * Determines if a given string contains a certain substring.
 */
bool strContains(const string inputStr, const string searchStr)
{
    size_t contains;

    contains = inputStr.find(searchStr);

    if (contains != string::npos)
        return true;
    else
        return false;
}

/*
 * Prints in standard output 'length' integer elements of a given array.
 */
void PrintArray(int *array, int length, string text)
{
    cout << text;
    for (int i = 0; i < length; i++)
    {
        cout << array[i] << " ";
    }
    cout << " " << endl;
}

/*
 * Prints in standard output 'length' long double elements of a given array.
 */
void PrintArray(long double *array, int length, string text)
{
    cout << text;
    for (int i = 0; i < length; i++)
    {
        cout << array[i] << " ";
    }
    cout << " " << endl;
}

/*
 * Prints in the standard output the 'size' elements of the given array..
 */
void PrintArray(int *array, int size)
{
    for (int i = 0; i < size; i++)
    {
        cout << array[i] << " ";
    }
    cout << " " << endl;
}

/*
 * Prints in standard output 'length' double elements of a given array.
 */
void PrintArray(double *array, int length, string text)
{
    int i;
    cout << text;
    for (i = 0; i < length; i++)
        printf("%3.5f,", array[i]);
    printf("\n ");
}

/*
 * Prints in the standard output given matrix.
 */
void PrintMatrix(int **matrix, int length, int length2, string text)
{
    int i, j;
    cout << text;
    for (i = 0; i < length; i++)
    {
        cout << "" << endl;
        for (j = 0; j < length2; j++)
        {
            cout << matrix[i][j] << " ";
        }
    }
    cout << " " << endl;
}

void PrintMatrix(double **matrix, int length, int length2, string text)
{
    int i, j;
    cout << text;
    for (i = 0; i < length; i++)
    {
        cout << "" << endl;
        for (j = 0; j < length2; j++)
        {
            cout << matrix[i][j] << ", ";
        }
    }
    cout << " " << endl;
}

/*
 * Prints in standard output 'lengthxlength' double elements of a given matrix.
 */
void PrintMatrixDouble(double **matrix, int length, int length2, string text)
{
    int i, j;
    cout << text << endl;
    for (i = 0; i < length; i++)
    {
        for (j = 0; j < length2; j++)
            printf("%3.9f, ", matrix[i][j]);
        printf("\n");
    }
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int *permutationA, int *permutationB, int size)
{
    int i, dist;
    int *v, *composition, *invertedB, *m_aux;

    dist = 0;
    v = new int[size - 1];
    composition = new int[size];
    invertedB = new int[size];
    // m_aux = new int[size];

    Invert(permutationB, size, invertedB);
    Compose(invertedB, permutationA, composition, size);
    vVector(v, composition, size);

    for (i = 0; i < size - 1; i++)
        dist += v[i];

    delete[] composition;
    delete[] invertedB;
    delete[] v;
    //delete[] m_aux;

    return dist;
}

/*
 * Calculates the tau Kendall distance between 2 permutations.
 */
int Kendall(int *permutationA, int *permutationB, int size, int *m_aux)
{
    int i, dist;
    int *v, *composition, *invertedB;

    dist = 0;
    v = new int[size - 1];
    composition = new int[size];
    invertedB = new int[size];

    Invert(permutationB, size, invertedB);
    Compose(permutationA, invertedB, composition, size);
    vVector_Fast(v, composition, size, m_aux);

    for (i = 0; i < size - 1; i++)
        dist += v[i];

    delete[] composition;
    delete[] invertedB;
    delete[] v;

    return dist;
}

/*
 * Calculates the Kendall tau distance between 2 permutations.
 * Auxiliary parameters are used for multiple continuous executions.
 */
int Kendall(int *permutationA, int *permutationB, int size, int *m_aux, int *invertedB, int *composition, int *v)
{
    int i, dist;
    dist = 0;
    Invert(permutationB, size, invertedB);
    Compose(permutationA, invertedB, composition, size);
    vVector_Fast(v, composition, size, m_aux);
    for (i = 0; i < size - 1; i++)
        dist += v[i];

    return dist;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int *permutationA, int *permutationB, int size)
{
    int *invertedB = new int[size];
    int *composition = new int[size];
    int *elemsToCycles = new int[size];
    int *maxPosInCycle = new int[size];
    int *freeCycle = new int[size];

    Invert(permutationB, size, invertedB);
    Compose(permutationA, invertedB, composition, size);

    int index, cycle, distance;

    for (int i = 0; i < size; i++)
    {
        elemsToCycles[i] = -1;
        maxPosInCycle[i] = -1;
        freeCycle[i] = 1;
    }

    while ((index = NextUnasignedElem(elemsToCycles, size)) != -1)
    {
        cycle = FindNewCycle(freeCycle, size);
        freeCycle[cycle] = 0;
        do
        {
            elemsToCycles[index] = cycle;
            index = composition[index]; //para permus de 1..n =>index = sigma[index]-1;
        } while (elemsToCycles[index] == -1);
    }
    distance = size - FindNewCycle(freeCycle, size);

    delete[] invertedB;
    delete[] composition;
    delete[] elemsToCycles;
    delete[] maxPosInCycle;
    delete[] freeCycle;

    return distance;
}

/*
 * Calculates the Cayley distance between 2 permutations.
 */
int Cayley(int *permutationA, int *permutationB, int size, int *invertedB, int *composition, int *elemsToCycles, int *maxPosInCycle, int *freeCycle)
{

    Invert(permutationB, size, invertedB);
    Compose(permutationA, invertedB, composition, size);

    int index, cycle, distance;

    for (int i = 0; i < size; i++)
    {
        elemsToCycles[i] = -1;
        maxPosInCycle[i] = -1;
        freeCycle[i] = 1;
    }

    while ((index = NextUnasignedElem(elemsToCycles, size)) != -1)
    {
        cycle = FindNewCycle(freeCycle, size);
        freeCycle[cycle] = 0;
        do
        {
            elemsToCycles[index] = cycle;
            index = composition[index]; //para permus de 1..n =>index = sigma[index]-1;
        } while (elemsToCycles[index] == -1);
    }
    distance = size - FindNewCycle(freeCycle, size);

    return distance;
}

int Hamming_distance(int *sigma1, int *sigma2, int len)
{
    int res = 0;
    for (int i = 0; i < len; i++)
    {
        if (sigma1[i] == sigma2[i])
        {
            res++;
        }
    }
    return res;
}

int FindNewCycle(int *freeCycle, int size)
{
    int i;

    for (i = 0; i < size; i++)
        if (freeCycle[i])
            return i;
    return size;
}

int NextUnasignedElem(int *elemsToCycles, int size)
{
    int i;
    for (i = 0; i < size; i++)
        if (elemsToCycles[i] == -1)
            return i;
    return -1;
}

/*
 * Calculates the length of the longest increasing subsequence in the given array of ints.
 */
int getLISLength(int *sigma, int size)
{

    // O(n log k)

    int i;
    vector<int> vc(1, sigma[0]);
    vector<int>::iterator vk;

    for (i = 1; i < size; i++)
    {
        for (vk = vc.begin(); vk != vc.end(); vk++)
            if (*vk >= sigma[i])
                break;
        if (vk == vc.end())
            vc.push_back(sigma[i]);
        else
            *vk = sigma[i];
    }

    return (int)vc.size();
}
/*
 * Implements the compose of 2 permutations of size n.
 */
void Compose(int *s1, int *s2, int *res, int n)
{
    int i;
    for (i = 0; i < n; i++)
        res[i] = s1[s2[i]];
}

/*
 * Calculates V_j-s vector.
 */
void vVector(int *v, int *permutation, int n)
{

    int i, j;
    for (i = 0; i < n - 1; i++)
        v[i] = 0;

    for (i = n - 2; i >= 0; i--)
        for (j = i + 1; j < n; j++)
            if (permutation[i] > permutation[j])
                v[i]++;
}

/*
 *  Optimized version proposed by Leti for the calculation of the V_j-s vector.
 */
void vVector_Fast(int *v, int *permutation, int n, int *m_aux)
{
    int i, j, index;

    for (i = 0; i < n - 1; i++)
    {
        v[i] = 0;
        m_aux[i] = 0;
    }
    m_aux[n - 1] = 0;
    for (j = 0; j < n - 1; j++)
    {
        index = permutation[j];
        v[j] = index - m_aux[index];
        for (i = index; i < n; i++)
            m_aux[i]++;
    }
}

/*
 * Inverts a permutation.
 */
void Invert(int *permu, int n, int *inverted)
{
    int i;
    for (i = 0; i < n; i++)
        inverted[permu[i]] = i;
}

/*
 * Applies the random keys sorting strategy to the vector of doubles
 */
void RandomKeys(int *a, double *criteriaValues, int size)
{
    bool *fixedValues = new bool[size];
    double criteria, min;
    int i, j;
    for (i = 0; i < size; i++)
    {
        fixedValues[i] = false;
        a[i] = 0;
    }
    int minPos = 0;
    for (i = 0; i < size; i++)
    {
        min = MAX_INTEGER;
        for (j = 0; j < size; j++)
        {
            criteria = criteriaValues[j];
            if (!fixedValues[j] && min > criteria)
            {
                min = criteria;
                minPos = j;
            }
        }

        fixedValues[minPos] = true;
        //a[i]=minPos;// modification por el asunto ordering /ranking
        a[minPos] = i; // original.
    }
    delete[] fixedValues;
}

/*
 * This method moves the value in position i to the position j.
 */
void InsertAt(int *array, int i, int j, int n)
{
    assert(0 <= i);
    assert(0 <= j);
    assert(i < n);
    assert(j < n);
    if (i != j)
    {
        int res[n];
        int val = array[i];
        if (i < j)
        {
            memcpy(res, array, sizeof(int) * i);

            for (int k = i + 1; k <= j; k++)
                res[k - 1] = array[k];

            res[j] = val;

            for (int k = j + 1; k < n; k++)
                res[k] = array[k];
        }
        else if (i > j)
        {
            memcpy(res, array, sizeof(int) * j);

            res[j] = val;

            for (int k = j; k < i; k++)
                res[k + 1] = array[k];

            for (int k = i + 1; k < n; k++)
                res[k] = array[k];
        }
        memcpy(array, res, sizeof(int) * n);
    }
}

/*
 * Calculates the factorial of a solution.
 */
long double factorial(int val)
{
    if (val <= 0)
        return 1;
    //long  N, b, c, p; // use int for fast calculation and small range of calculation..
    long b, c;
    long double p, N;
    N = (long double)val;
    c = (long)N - 1;
    p = 1;
    while (c > 0)
    {
        p = 0;
        b = c;
        while (b > 0)
        {
            if (b & 1)
            {
                p += N; // p = p + N;
            }
            // if you would like to use double choose the alternative forms instead shifts
            // the code is fast even!
            // you can use the same tips on double or 64 bit int etc.... but you must... ;-)
            //b >>= 1; // b/=2; (b = b / 2;) ( b >> 1; a.s.r. is more efficent for int or long..!)
            b /= 2;
            //N <<= 1; // N += N; N = N + N; N = N * 2; (N <<=1; a.s.l. is more efficent for int or long..!)
            N += N;
        } // end of: while(b>0)
        N = p;
        c--; // c = c - 1;
    }        // end of: while(c > 0)
    //printf("[%d] is the factorial! \n", p);
    return p;
}

/*
 * This method applies a swap of the given i,j positions in the array.
 */
void swap_two_positions(int *array, int i, int j)
{
    int aux = array[i];
    array[i] = array[j];
    array[j] = aux;
}

static double tt_tic = 0;

double getTick()
{
    struct timespec ts;
    double theTick;
    clock_gettime(CLOCK_REALTIME, &ts);
    theTick = (double)ts.tv_nsec / 1000000000.0;
    theTick += (double)ts.tv_sec;
    return theTick;
}

void tic() { tt_tic = getTick(); }

double toc() { return (getTick() - tt_tic); }

int random_integer_fast(int min, int max)
{
    return min + (rand() % static_cast<int>(max - min + 1));
}

// https://ericlippert.com/2013/12/16/how-much-bias-is-introduced-by-the-remainder-technique/
int random_integer_uniform(int min, int max)
{

    //LOG->write("min: ", false);
    //LOG->write(min);
    //LOG->write("max: ", false);
    //LOG->write(max);

    if (max == 0)
    {
        //LOG->write("MAX WAS 0");
        int range = min;
        while (true)
        {
            int value = rand();
            if (value < RAND_MAX - (RAND_MAX % range))
            {

                //LOG->write("range: ", false);
                //LOG->write(range);
                //LOG->write("value mod range: ", false);
                //LOG->write(value % range);

                return value % range;
            }
        }
    }
    else
    {
        int r = random_integer_uniform(max - min) + min;
        return r;
    }
}

// chooses a random integer from {0,1,2, range_max - 1}
int random_range_integer_uniform(int range_max)
{
    return random_integer_uniform(range_max);
}

double random_0_1_float()
{
    return (double)rand() / RAND_MAX;
}

double sigmoid(double x)
{
    return 1.0 / (1.0 + exp(-x));
}

int choose_index_given_probabilities(double *probabilities_array, int len)
{
    double r = random_0_1_float();
    double cum_prob = 0;

    for (int i = 0; i < len; i++)
    {
        cum_prob += probabilities_array[i];
        if (r < cum_prob)
        {
            return i;
        }
    }

    return choose_index_given_probabilities(probabilities_array, len);

    //cout << endl;
    //cout << "cum_prob = " << cum_prob << endl;
    //assert(cum_prob > 0.99999);
}

int choose_index_given_weights(double *weights_array, int len)
{
    double r = random_0_1_float();
    double cum_sum = 0;
    double total = 0;

    for (size_t i = 0; i < len; i++)
    {
        total += weights_array[i];
    }

    r = total * r;

    for (int i = 0; i < len; i++)
    {
        cum_sum += weights_array[i];
        if (r < cum_sum)
        {
            return i;
        }
    }

    return choose_index_given_weights(weights_array, len);

    //cout << endl;
    //cout << "cum_prob = " << cum_prob << endl;
    //assert(cum_prob > 0.99999);
}

bool coin_toss(double p_of_true)
{
    if (random_0_1_float() < p_of_true)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int tools_round(double x)
{
    if (x <= 0.0)
    {
        return (int)(x - 0.5);
    }
    else
    {
        return (int)(x + 0.5);
    }
}

void shuffle_vector(int *vec, int len)
{
    for (int i = 0; i < len - 1; i++)
    {
        int pos = rand() % (len - i) + i;
        //int pos = (int) (unif_rand() * (len-i) + i);
        int aux = vec[i];
        vec[i] = vec[pos];
        vec[pos] = aux;
    }
}

int *_random_permu;
bool _permu_was_set = false;

void select_indexes_to_insert_towards(int *genome_to_be_changed, int *reference_genome, int *i, int *j, int n)
{

    if (!_permu_was_set)
    {
        _random_permu = new int[n];
        GenerateRandomPermutation(_random_permu, n);
        _permu_was_set = true;
    }
    else
    {
        shuffle_vector(_random_permu, n);
    }

    for (int k = 0; k < n; i++)
    {
        int index = _random_permu[k];
        if (genome_to_be_changed[index] != reference_genome[index])
        {
            *i = index;
            *j = Find(reference_genome, n, genome_to_be_changed[index]);
            return;
        }
    }

    *i = 0;
    *j = 0;
    return;
}