#include <stdlib.h> /* srand, rand */
#include "Tools.h"
#include <float.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

//extern SimpleLog *LOG;

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
            if (value < RAND_MAX - RAND_MAX % range)
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
        assert(max > min);
        int range = static_cast<int>(max - min + 1);
        while (true)
        {
            int value = rand();
            if (value < RAND_MAX - RAND_MAX % range)
            {
                return min + value % range;
            }
        }
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

int chose_index_given_probabilities(double *probabilities_array, int len)
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

    return chose_index_given_probabilities(probabilities_array, len);

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






                                                



