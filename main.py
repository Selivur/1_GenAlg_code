from config import env, NR
from model.fitness_functions import *
from selection.rws import *
from selection.sus import *
from model.encoding import *
from model.gen_operators import *
from runner import run_experiment
from datetime import datetime
from output import excel
import time

fitness_functions = [
    #(FconstALL(100), 'FconstALL'),
    #(Fx2(FloatEncoder(0.0, 10.23, 10)), 'Fx2'),
    #(Fx2(FloatEncoder(0.0, 10.23, 10, is_gray=True)), 'Fx2_gray'),
    #(F5122subx2(FloatEncoder(-5.12, 5.11, 10)), 'F5122subx2'),
    #(F5122subx2(FloatEncoder(-5.12, 5.11, 10, is_gray=True)), 'F5122subx2_gray'),
    #(Fexp(0.25, FloatEncoder(0.0, 10.23, 10)), 'Fexp0.25'),
    #(Fexp(0.25, FloatEncoder(0.0, 10.23, 10, is_gray=True)), 'Fexp0.25_gray'),
    #(Fexp(1, FloatEncoder(0.0, 10.23, 10)), 'Fexp1'),
    #(Fexp(1, FloatEncoder(0.0, 10.23, 10, is_gray=True)), 'Fexp1_gray'),
    (Fexp(2, FloatEncoder(0.0, 10.23, 10)), 'Fexp2'),
    (Fexp(2, FloatEncoder(0.0, 10.23, 10, is_gray=True)), 'Fexp2_gray'),

    #(FH(100), 'FH'),
    (Rastrigin(FloatEncoder(-5.12, 5.11, 10)), 'Rastrigin_binary'),
    #(Rastrigin(FloatEncoder(-5.12, 5.11, 10, is_gray=True)), 'Rastrigin_gray'),
    #(DebsTestFunction2(FloatEncoder(0.0, 1.023, 10)), 'Debs2_binary'),
    #(DebsTestFunction2(FloatEncoder(0.0, 1.023, 10, is_gray=True)), 'Debs2_gray'),
    #(DebsTestFunction4(FloatEncoder(0.0, 1.023, 10)), 'Debs4_binary'),
    #(DebsTestFunction4(FloatEncoder(0.0, 1.023, 10, is_gray=True)), 'Debs4_gray'),
]
selection_methods = [
    ##(RWS, 'RWS', None),
    ##(TruncationRWS, 'TruncationRWS_2', 2),
    #(TruncationRWS, 'TruncationRWS_4', 4),
    (PowerScalingRWS, 'PowerScalingRWS_0.8', 0.8),
    #(PowerScalingRWS, 'PowerScalingRWS_0.9', 0.9),
    #(PowerScalingRWS, 'PowerScalingRWS_1.0', 1.0),
    #(PowerScalingRWS, 'PowerScalingRWS_1.1', 1.1),
    #(PowerScalingRWS, 'PowerScalingRWS_1.2', 1.2),
    ##(SUS, 'SUS', None),
    ##(TruncationSUS, 'TruncationSUS_2', 2),
    #(TruncationSUS, 'TruncationSUS_4', 4),
    (PowerScalingSUS, 'PowerScalingSUS_0.8', 0.8),
    #(PowerScalingSUS, 'PowerScalingSUS_0.9', 0.9),
    #(PowerScalingSUS, 'PowerScalingSUS_1.0', 1.0),
    #(PowerScalingSUS, 'PowerScalingSUS_1.1', 1.1),
    #(PowerScalingSUS, 'PowerScalingSUS_1.2', 1.2),

    ##(ExpRankingRWS, 'ExpRankingRWS_0', 0),
    #(ExpRankingRWS, 'ExpRankingRWS_1', 1),
    #(ExpRankingRWS, 'ExpRankingRWS_2', 2),
    #(ExpRankingRWS, 'ExpRankingRWS_3', 3),
    ##(ExpRankingSUS, 'ExpRankingSUS_0', 0),
    #(ExpRankingSUS, 'ExpRankingSUS_1', 1),
    #(ExpRankingSUS, 'ExpRankingSUS_2', 2),
    #(ExpRankingSUS, 'ExpRankingSUS_3', 3)
]
gen_operators = [
    (BlankGenOperator, 'no_operators'),
    #(Crossover, 'crossover'),
    #(Mutation, 'mutation'),
    #(CrossoverAndMutation, 'crossover_and_mutation')
]

# a list of tuples of parameters for each run that involves a certain fitness function 
# {fitness_func_name: [(tuples with run parameters), (), ..., ()], other_func: [], ...}
experiment_params = {
    ff: [
        (sm, sm_arg, go, (ff_name, sm_name, go_name))
        for (sm, sm_name, sm_arg) in selection_methods
        for (go, go_name) in gen_operators
    ] for (ff, ff_name) in fitness_functions
}

# only keeping one list of populations in memory at a time (for one fitness function)
def generate_all_populations_for_fitness_function(ff):
    return [ff.generate_population_for_run(run_i) for run_i in range(NR)]

def log(x):
    datetime_prefix = str(datetime.now())[:-4]
    print(f'{datetime_prefix} | {x}')

if __name__ == '__main__':
    log('Program start')
    print('----------------------------------------------------------------------')
    start_time = time.time()
    results = []

    for ff in experiment_params:
        ff_start_time = time.time()
        populations = generate_all_populations_for_fitness_function(ff)
        params = [params + (populations,) for params in experiment_params[ff]]
        experiment_stats_list = [run_experiment(*p) for p in params]

        excel.write_ff_stats(experiment_stats_list)
        for experiment_stats in experiment_stats_list:
            del experiment_stats.runs
            results.append(experiment_stats)

        ff_end_time = time.time()
        ff_name = experiment_params[ff][0][3][0]
        log(f'{ff_name} experiments finished in {(ff_end_time - ff_start_time):.2f}s')

    #excel.write_aggregated_stats(results)

    print('----------------------------------------------------------------------')
    end_time = time.time()
    log(f'Program end. Total runtime: {end_time - start_time:.2f}s')
