import os
import subprocess
from multiprocessing import Pool

from statistics import mean, median
import matplotlib.pyplot as plt
import numpy as np

from skopt.space import Real, Integer
from skopt.utils import use_named_args

from skopt import gp_minimize

from numpy import argmax


from skopt.plots import plot_objective, plot_evaluations


INSTANCE_PATH = 'QAPinstances/tai31a.dat.dat'
BINARY_PATH = 'binary'


def _prep_cmd(arg_list):

    SEED, POPSIZE, EXPO = arg_list

    cmd = "./"
    cmd += BINARY_PATH + ' '
    cmd += str(INSTANCE_PATH.split("/")[-1]) + ' '
    cmd += str(SEED) + ' '
    cmd += str(POPSIZE) + ' '
    cmd += str(EXPO) + ' < '
    cmd += INSTANCE_PATH
    return cmd



def black_box_f_individual(arg_list):



    cmd = _prep_cmd(arg_list)

    # subprocess.check_output("./binary  tai 2 500 0.5 0.8 < QAPinstances/tai125e01.dat.dat ", shell=True)
    best = 1945072
    return (abs(int(subprocess.check_output(cmd, shell=True))) - best) / best * 100

    




def black_box_function(POPSIZE, EXPO): 

    results = list()
    N_OF_EVALS = 600

    POPSIZE = int(POPSIZE)
    EXPO = float(EXPO)

    args = [[np.random.randint(10000), POPSIZE, EXPO] for i in range(N_OF_EVALS)]


    pool = Pool(processes=6)            
    results = pool.map(black_box_f_individual, args)
    pool.terminate()


    print("best: ",max(results), " | average: ", mean(results), args[0])


    return  median(results)




if __name__ == '__main__':


    # The list of hyper-parameters we want to optimize.

    space  =  [
              Real(50, 1000, name='POPSIZE'),
              Real(0.05, 6.0, name='EXPO')
              ]

    # tai75 [mid=0.35, FIN=8.0, EXPO=0.45]


    # estimate variance
    # space  = [Integer(15000, 15001, name='POPSIZE'),
    #           Integer(700, 701, name='MAX_ITERATIONS'),
    #           Real(0.7, 0.71, name='TARGET_MID_EXPECTATION_PERCENTAGE'),
    #           Real(5.0, 5.1, name='FINAL_EXPECTATION'),
    #           Integer(10, 11, name='T_ADD'),
    #           Integer(25, 26, name='FAST_INCREASE_EACH_ITERATION_COUNTS_AS'),
    #           ]


    @use_named_args(space)
    def objective(**params):

        return black_box_function(**params)


   

    # this decorator allows your objective function to receive a the parameters as
    # keyword arguments. This is particularly convenient when you want to set scikit-learn
    # estimator parameters

    res_gp = gp_minimize(objective,
                         space,
                         n_random_starts=6,
                         n_calls=12,
                         random_state=0,
                         verbose=True,
                         noise=1,) 

    _ = plot_objective(res_gp)
    plt.tight_layout(pad=0.2)
    plt.savefig("bo_heatmap.png", figsize=(60,80))
    _ = plot_evaluations(res_gp)
    plt.tight_layout(pad=0.2)
    plt.savefig("bo_sampled_points.png", figsize=(60,80))
    
    print(res_gp.x)


# Parameters estimated on instance tai31a.dat
# [972.0110483144621, 5.143023755631278]




    
