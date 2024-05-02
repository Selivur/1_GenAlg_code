import numpy as np
from config import N, EPS, N_LAST_GENS
from model.chromosome import Chromosome
from copy import deepcopy, copy


class Population:
    def __init__(self, fitness_function, seed=0,
                 chromosomes=None,
                 optimal_fraction=0,
                 is_single_optimal=True,
                 use_normal_distribution=False):
        self.fitness_function = fitness_function

        if chromosomes is not None:
            self.chromosomes = chromosomes
        else:
            self.chromosomes = np.empty(N, dtype=object)
            rng = np.random.default_rng(seed=seed)

            if is_single_optimal:
                num_optimal = 1
            else:
                num_optimal = max(1, int(N * optimal_fraction))

            if use_normal_distribution:
                for chr_i in range(N):
                    genotype = rng.binomial(1, 0.5, fitness_function.chr_length).astype('b1')
                    while genotype != fitness_function.get_optimal():
                        genotype = rng.binomial(1, 0.5, fitness_function.chr_length).astype('b1')
                    self.chromosomes[chr_i] = Chromosome(chr_i, genotype, fitness_function)
            else:
                for chr_i in range(num_optimal):
                    self.chromosomes[chr_i] = copy(fitness_function.get_optimal())
                for chr_i in range(num_optimal, N):
                    genotype = rng.choice([b'0', b'1'], fitness_function.chr_length)
                    self.chromosomes[chr_i] = Chromosome(chr_i, genotype, fitness_function)

        # Shuffle population.
        np.random.shuffle(self.chromosomes)
        self.update()

    def has_converged(self, param_names):
        has_gen_op = param_names[2] != 'no_operators'

        if not has_gen_op:
            return self.is_homogenous_100()

        return self.is_homogenous_99()
    
    def has_f_avg_converged(self, f_avgs):
        if len(f_avgs) < N_LAST_GENS:
            return False

        diffs = []
        for i in range(1, len(f_avgs)):
            curr = f_avgs[i]
            prev = f_avgs[i-1]
            diffs.append(abs(curr - prev))

        return all(x <= EPS for x in diffs)

    def count_unique_chromosomes(self):
        # Отримуємо список усіх генотипів у популяції
        all_genotypes = [chr.genotype for chr in self.chromosomes]

        # Використовуємо набір (set) для фільтрації унікальних генотипів
        unique_genotypes = set(map(tuple, all_genotypes))  # Використовуємо map, щоб перетворити генотипи у кортежі

        # Повертаємо кількість унікальних генотипів
        return len(unique_genotypes)

    def is_homogenous_90(self):
        l = self.fitness_function.chr_length
        for i in range(l):
            n_zeros = len([True for g in self.genotypes if g[i] == b'0'])
            percentage = n_zeros / N
            if 0.01 < percentage < 0.9:
                return False
        return True

    def is_homogeneous_frac(self, frac):
        """
        check if the population is homogenous by at least (frac*100) percent
        :param frac: 0.5 < frac < 1 - fraction of the population that should have the same value for every gene
        :return: True, if the population satisfies the statement above, False, otherwise
        """
        l = self.fitness_function.chr_length
        for i in range(l):
            n_zeros = len([True for g in self.genotypes if g[i] == b'0'])
            percentage = n_zeros / N
            if percentage > (1 - frac) and percentage < frac:
                return False
        return True

    def get_unique_X(self):
        '''
        Find the number of different chromosomes in the population
        '''
        unique_genotypes = set(tuple(genotype) for genotype in self.genotypes)
        return len(unique_genotypes)

    def is_homogenous_99(self):
        l = self.fitness_function.chr_length
        for i in range(l):
            n_zeros = len([True for g in self.genotypes if g[i] == b'0'])
            percentage = n_zeros / N
            if 0.01 < percentage < 0.99:
                return False
        return True

    def is_homogenous_100(self):
            return all([np.array_equal(geno, self.genotypes[0]) for geno in self.genotypes[1:]])

    def found_close_to_optimal(self):
        for chr in self.chromosomes:
            if self.fitness_function.check_chromosome_success(chr):
                return True
        return False

    def get_fitness_max(self):
        res = np.max(self.fitnesses)
        return res

    def get_fitness_avg(self):
        return np.mean(self.fitnesses)

    def get_fitness_std(self):
        return np.std(self.fitnesses)
    
    def count_fitness_at_least(self, min_fitness):
        return len([True for f in self.fitnesses if f >= min_fitness])

    def count_optimal_genotype(self):
        optimal = self.fitness_function.get_optimal().genotype
        return len([True for g in self.genotypes if np.array_equal(g, optimal)])

    def get_ids(self):
        return [chr.id for chr in self.chromosomes]

    def update(self):
        if any(chr is None for chr in self.chromosomes):
            raise ValueError("One or more chromosomes is None")
        self.fitnesses = np.array([chr.fitness for chr in self.chromosomes])
        self.genotypes = np.array([chr.genotype for chr in self.chromosomes])

    def update_chromosomes(self, chromosomes):
        self.chromosomes = chromosomes
        self.update()
    
    def __deepcopy__(self, memo):
        return Population(self.fitness_function, chromosomes=deepcopy(self.chromosomes))
    
    def __str__(self):
        return str(np.array([str(chr) for chr in self.chromosomes]))

    def get_optimal_individual_id(self):
        optimal_genotype = self.fitness_function.get_optimal().genotype
        for chr in self.chromosomes:
            if np.array_equal(chr.genotype, optimal_genotype):
                return chr.id
        return None
