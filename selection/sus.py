from config import G, N, c_values
from model.population import Population
import numpy as np
from selection.selection_method import SelectionMethod
from copy import copy


class SUS(SelectionMethod):
    def select(self, population: Population):
        fitness_sum = 0
        fitness_scale = []

        for index, chromosome in enumerate(population.chromosomes):
            fitness_sum += chromosome.fitness
            if index == 0:
                fitness_scale.append(chromosome.fitness)
            else:
                fitness_scale.append(chromosome.fitness + fitness_scale[index - 1])

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            fitness_scale = [0.0001 * (i+1) for i in range(N)]

        mating_pool = self.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)
    
    @staticmethod
    def basic_sus(population: Population, fitness_sum, fitness_scale):
        mating_pool = np.empty(N, dtype=object)
        fitness_step = fitness_sum / N
        random_offset = np.random.uniform(0, fitness_step)
        current_fitness_pointer = random_offset
        last_fitness_scale_position = 0

        for i in range(N):
            for fitness_scale_position in range(last_fitness_scale_position, len(fitness_scale)):
                if fitness_scale[fitness_scale_position] >= current_fitness_pointer:
                    mating_pool[i] = copy(population.chromosomes[fitness_scale_position])
                    last_fitness_scale_position = fitness_scale_position
                    break
            current_fitness_pointer += fitness_step

        return mating_pool


class DisruptiveSUS(SelectionMethod):
    def select(self, population: Population):
        fitness_sum = 0
        fitness_scale = []
        f_avg = population.get_fitness_avg()

        for index, chromosome in enumerate(population.chromosomes):
            f_scaled = abs(chromosome.fitness - f_avg)
            fitness_sum += f_scaled
            if index == 0:
                fitness_scale.append(f_scaled)
            else:
                fitness_scale.append(f_scaled + fitness_scale[index - 1])

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            fitness_scale = [0.0001 * (i+1) for i in range(N)]

        mating_pool = SUS.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)


class BlendedSUS(SelectionMethod):
    def __init__(self):
        self.i = 0

    def select(self, population: Population):
        fitness_sum = 0
        fitness_scale = []

        for index, chromosome in enumerate(population.chromosomes):
            f_scaled = chromosome.fitness / (G + 1 - self.i)
            fitness_sum += f_scaled
            if index == 0:
                fitness_scale.append(f_scaled)
            else:
                fitness_scale.append(f_scaled + fitness_scale[index - 1])

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            fitness_scale = [0.0001 * (i+1) for i in range(N)]

        mating_pool = SUS.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)

        self.i += 1


class WindowSUS(SelectionMethod):
    def __init__(self, h=2):
        self.h = h
        self.f_h_worst = []

    def select(self, population: Population):
        if len(self.f_h_worst) < self.h:
            self.f_h_worst.append(min(population.fitnesses))
        else:
            self.f_h_worst[0] = self.f_h_worst[1]
            self.f_h_worst[1] = min(population.fitnesses)
        f_worst = min(self.f_h_worst)

        fitness_sum = 0
        fitness_scale = []

        for index, chromosome in enumerate(population.chromosomes):
            f_scaled = chromosome.fitness - f_worst
            fitness_sum += f_scaled
            if index == 0:
                fitness_scale.append(f_scaled)
            else:
                fitness_scale.append(f_scaled + fitness_scale[index - 1])

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            fitness_scale = [0.0001 * (i+1) for i in range(N)]

        mating_pool = SUS.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)


class ExpRankingSUS(SelectionMethod):
    def __init__(self, c_variant):
        self.c_variant = c_variant
        self.c = c_values[min(i for i in c_values.keys() if N <= i)][self.c_variant]

    def select(self, population: Population):
        population.chromosomes = sorted(population.chromosomes, key=lambda chr: chr.fitness)

        fitness_sum = 0
        fitness_scale = []

        for index, chromosome in enumerate(population.chromosomes):
            f_scaled = ((self.c - 1) / (self.c ** N - 1)) * self.c ** (N - (index + 1))
            fitness_sum += f_scaled
            if index == 0:
                fitness_scale.append(f_scaled)
            else:
                fitness_scale.append(f_scaled + fitness_scale[index - 1])

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            fitness_scale = [0.0001 * (i + 1) / N for i in range(N)]

        mating_pool = SUS.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)

class PowerScalingSUS(SelectionMethod):
    def __init__(self, k_start=0.8, k_end=1.2):
        self.k_start = k_start
        self.k_end = k_end
        self.k = k_start

    def select(self, population: Population):
        if np.median(population.fitnesses) >= np.mean(population.fitnesses):
            self.k = self.k_end

        fitnesses = [f ** self.k for f in population.fitnesses]
        fitness_sum = sum(fitnesses)
        fitness_scale = np.cumsum(fitnesses)

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            fitnesses = [0.0001 * (i + 1) for i in range(N)]
            fitness_scale = np.cumsum(fitnesses)

        mating_pool = SUS.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)

class TruncationSUS(SelectionMethod):
    def __init__(self, c, sigma=2):
        self.c = c
        self.sigma = sigma

    def select(self, population: Population):
        fitnesses = population.fitnesses
        fitness_mean = np.mean(fitnesses)
        fitness_std = np.std(fitnesses)

        f_scaled = [f - (fitness_mean - self.c * fitness_std) for f in fitnesses]
        f_scaled = [f if f >= 0 else 0 for f in f_scaled]
        fitness_sum = sum(f_scaled)
        fitness_scale = np.cumsum(f_scaled)

        if fitness_sum == 0:
            fitness_sum = 0.0001 * N
            f_scaled = [0.0001 * (i+1) for i in range(N)]
            fitness_scale = np.cumsum(f_scaled)

        mating_pool = SUS.basic_sus(population, fitness_sum, fitness_scale)
        population.update_chromosomes(mating_pool)
