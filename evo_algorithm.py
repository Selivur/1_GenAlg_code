from config import *
from model.fitness_functions import *
from selection.selection_method import SelectionMethod
from model.gen_operators import GeneticOperator
from stats.run_stats import RunStats
from stats.generation_stats import GenerationStats
from output import plotting


class EvoAlgorithm:
    def __init__(self,
                 initial_population: Population,
                 selection_method: SelectionMethod,
                 genetic_operator: GeneticOperator,
                 param_names: tuple[str]):
        self.population: Population = initial_population
        self.selection_method = selection_method
        self.genetic_operator = genetic_operator
        self.param_names = param_names

        self.gen_i = 0
        self.run_stats = RunStats(self.param_names)
        self.prev_gen_stats = None
        self.gen_stats_list = None
        self.has_converged = False
        
    def run(self, run_i):
        if run_i < RUNS_TO_PLOT:
            self.gen_stats_list = []

        f_avgs = []
        while not self.has_converged and self.gen_i < G:
            gen_stats = self.__calculate_stats_and_evolve(run_i)

            f_avgs.append(gen_stats.f_avg)
            if len(f_avgs) > N_LAST_GENS:
                f_avgs.pop(0)

            self.has_converged = (
                self.population.has_converged(self.param_names))
            self.prev_gen_stats = gen_stats
            self.gen_i += 1

        gen_stats = self.__calculate_final_stats(run_i)
        self.run_stats.NI = self.gen_i
        self.run_stats.is_successful = self.__check_success(gen_stats)

        if run_i < RUNS_TO_PLOT:
            plotting.plot_generation_stats(self.population, self.param_names, run_i, self.gen_i)
            plotting.plot_run_stats(self.gen_stats_list, self.param_names, run_i)

        return self.run_stats

    def __calculate_stats_and_evolve(self, run_i):
        if run_i < RUNS_TO_PLOT and self.gen_i < DISTRIBUTIONS_TO_PLOT:
            plotting.plot_generation_stats(self.population, self.param_names, run_i, self.gen_i)
        
        gen_stats = GenerationStats(self.population, self.param_names)
        if run_i < RUNS_TO_PLOT:
            self.gen_stats_list.append(gen_stats)

        gen_stats.calculate_stats_before_selection(self.prev_gen_stats)
        self.selection_method.select(self.population)
        gen_stats.calculate_stats_after_selection()
        self.run_stats.update_stats_for_generation(gen_stats, self.gen_i)
        self.genetic_operator.apply(self.population)

        return gen_stats

    def __calculate_final_stats(self, run_i):
        if run_i < RUNS_TO_PLOT and self.gen_i < DISTRIBUTIONS_TO_PLOT:
            plotting.plot_generation_stats(self.population, self.param_names, run_i, self.gen_i)

        gen_stats = GenerationStats(self.population, self.param_names)
        if run_i < RUNS_TO_PLOT:
            self.gen_stats_list.append(gen_stats)

        gen_stats.calculate_stats_before_selection(self.prev_gen_stats)
        self.run_stats.update_final_stats(gen_stats, self.gen_i)

        return gen_stats

    def __check_success(self, gen_stats: GenerationStats):
        has_gen_op = self.param_names[2] != 'no_operators'
        is_FconstALL = self.param_names[0] == 'FconstALL'
        is_binary_chain = self.param_names[0] == 'FH'
        is_real_arg = not is_binary_chain

        if is_FconstALL and not has_gen_op:
            return self.has_converged
        elif is_FconstALL and has_gen_op:
            return self.has_converged and self.population.is_homogenous_90()
        elif is_binary_chain and not has_gen_op:
            return self.has_converged and self.population.count_optimal_genotype() == N
        elif is_binary_chain and has_gen_op:
            return self.has_converged and self.population.count_optimal_genotype() >= N * 0.9
        elif is_real_arg:
            return self.has_converged and self.population.found_close_to_optimal()

