from config import N
from stats.generation_stats import GenerationStats
import numpy as np

class RunStats:
    def __init__(self, param_names: tuple[str]):
        self.param_names = param_names

        # Criteria after first selection
        self.I_start = None
        self.GR_start = None
        self.Pr_start = None

        # Criteria over iterations
        self.Pr_min = None
        self.NI_Pr_min = None
        self.Pr_max = None
        self.NI_Pr_max = None
        self.Pr_avg = None
        self.NI_Pr_avg = None

        self.NI = None
        self.F_found = None
        self.F_avg = None
        self.is_successful = False

        # Reproduction Rate
        self.RR_min = None
        self.NI_RR_min = None
        self.RR_max = None
        self.NI_RR_max = None
        self.RR_avg = None

        # Loss of Diversity
        self.Teta_min = None
        self.NI_Teta_min = None
        self.Teta_max = None
        self.NI_Teta_max = None
        self.Teta_avg = None

        # Selection Intensity
        self.I_min = None
        self.NI_I_min = None
        self.I_max = None
        self.NI_I_max = None
        self.I_avg = None

        # Selection Difference
        self.s_min = None
        self.NI_s_min = None
        self.s_max = None
        self.NI_s_max = None
        self.s_avg = None

        # Growth Rate
        self.GR_early = None
        self.GR_late = None
        self.NI_GR_late = None
        self.GR_avg = None

        # New val
        self.NI_loose = 0
        self.Num_loose = 0

        # Нові поля для втраченої оптимальної особини
        self.Avg_NI_loose = None
        self.Sigma_NI_loose = None
        self.Avg_Num_loose = None
        self.Sigma_Num_loose = None

    def update_stats_for_generation(self, gen_stats: GenerationStats, gen_i):
        # Reproduction Rate
        if self.RR_min is None or gen_stats.reproduction_rate < self.RR_min:
            self.RR_min = gen_stats.reproduction_rate
            self.NI_RR_min = gen_i
        if self.RR_max is None or gen_stats.reproduction_rate > self.RR_max:
            self.RR_max = gen_stats.reproduction_rate
            self.NI_RR_max = gen_i
        if self.RR_avg is None:
            self.RR_avg = gen_stats.reproduction_rate
        else:
            self.RR_avg = (self.RR_avg * (gen_i - 1) + gen_stats.reproduction_rate) / gen_i

        # Loss of Diversity
        if self.Teta_min is None or gen_stats.loss_of_diversity < self.Teta_min:
            self.Teta_min = gen_stats.loss_of_diversity
            self.NI_Teta_min = gen_i
        if self.Teta_max is None or gen_stats.loss_of_diversity > self.Teta_max:
            self.Teta_max = gen_stats.loss_of_diversity
            self.NI_Teta_max = gen_i
        if self.Teta_avg is None:
            self.Teta_avg = gen_stats.loss_of_diversity
        else:
            self.Teta_avg = (self.Teta_avg * (gen_i - 1) + gen_stats.loss_of_diversity) / gen_i

        if self.param_names[0] != 'FconstALL':
            # Selection Intensity
            if self.I_min is None or gen_stats.intensity < self.I_min:
                self.I_min = gen_stats.intensity
                self.NI_I_min = gen_i
            if self.I_max is None or gen_stats.intensity > self.I_max:
                self.I_max = gen_stats.intensity
                self.NI_I_max = gen_i
            if self.I_avg is None:
                self.I_avg = gen_stats.intensity
            else:
                self.I_avg = (self.I_avg * (gen_i - 1) + gen_stats.intensity) / gen_i

            # Selection Difference
            if self.s_min is None or gen_stats.difference < self.s_min:
                self.s_min = gen_stats.difference
                self.NI_s_min = gen_i
            if self.s_max is None or gen_stats.difference > self.s_max:
                self.s_max = gen_stats.difference
                self.NI_s_max = gen_i
            if self.s_avg is None:
                self.s_avg = gen_stats.difference
            else:
                self.s_avg = (self.s_avg * (gen_i - 1) + gen_stats.difference) / gen_i

            # Update criteria after first selection
            if gen_i == 1:
                self.I_start = gen_stats.intensity
                self.GR_start = gen_stats.growth_rate
                self.Pr_start = gen_stats.f_best / gen_stats.f_avg

            # Update Pr criteria
            if self.Pr_min is None or gen_stats.f_best / gen_stats.f_avg < self.Pr_min:
                self.Pr_min = gen_stats.f_best / gen_stats.f_avg
                self.NI_Pr_min = gen_i
            if self.Pr_max is None or gen_stats.f_best / gen_stats.f_avg > self.Pr_max:
                self.Pr_max = gen_stats.f_best / gen_stats.f_avg
                self.NI_Pr_max = gen_i
            if self.Pr_avg is None:
                self.Pr_avg = gen_stats.f_best / gen_stats.f_avg
            else:
                self.Pr_avg = (self.Pr_avg * (gen_i - 1) + gen_stats.f_best / gen_stats.f_avg) / gen_i
            self.NI_Pr_avg = gen_i

            # Growth Rate
            if gen_i == 2:
                self.GR_early = gen_stats.growth_rate
            if self.GR_late is None and gen_stats.num_of_best >= N / 2:
                self.GR_late = gen_stats.growth_rate
                self.NI_GR_late = gen_i
            if self.GR_avg is None:
                self.GR_avg = gen_stats.growth_rate
            else:
                self.GR_avg = (self.GR_avg * (gen_i - 1) + gen_stats.growth_rate) / gen_i
        if gen_stats.optimal_individual_lost:
            self.Num_loose += 1
            self.NI_loose = gen_i
            # Оновлення статистики для втраченої оптимальної особини
            if self.Num_loose > 0:
                if self.Avg_NI_loose is None:
                    self.Avg_NI_loose = gen_i
                else:
                    self.Avg_NI_loose = (self.Avg_NI_loose * (self.Num_loose - 1) + gen_i) / self.Num_loose

                if self.Sigma_NI_loose is None:
                    self.Sigma_NI_loose = 0
                else:
                    self.Sigma_NI_loose = np.sqrt(((gen_i - self.Avg_NI_loose) ** 2 + self.Sigma_NI_loose ** 2 * (
                                self.Num_loose - 1)) / self.Num_loose)

                if self.Avg_Num_loose is None:
                    self.Avg_Num_loose = self.Num_loose
                else:
                    self.Avg_Num_loose = (self.Avg_Num_loose * (self.Num_loose - 1) + self.Num_loose) / self.Num_loose

                if self.Sigma_Num_loose is None:
                    self.Sigma_Num_loose = 0
                else:
                    self.Sigma_Num_loose = np.sqrt(((
                                                                self.Num_loose - self.Avg_Num_loose) ** 2 + self.Sigma_Num_loose ** 2 * (
                                                                self.Num_loose - 1)) / self.Num_loose)

    def update_final_stats(self, gen_stats: GenerationStats, gen_i):
        if self.param_names[0] != 'FconstALL':
            self.F_found = gen_stats.f_best
            self.F_avg = gen_stats.f_avg

            if gen_i == 2:
                self.GR_early = gen_stats.growth_rate
            if self.GR_late is None and gen_stats.num_of_best >= N / 2:
                self.GR_late = gen_stats.growth_rate
                self.NI_GR_late = gen_i
            if self.GR_avg is None:
                self.GR_avg = gen_stats.growth_rate
            else:
                self.GR_avg = (self.GR_avg * (gen_i - 1) + gen_stats.growth_rate) / gen_i
