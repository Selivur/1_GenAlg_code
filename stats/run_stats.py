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

        # Нові поля для зберігання втрат оптимальної особини
        self.NI_loose = 0
        self.Num_loose = 0

        # Нові поля для втраченої оптимальної особини
        self.Avg_NI_loose = None
        self.Sigma_NI_loose = None
        self.Avg_Num_loose = None
        self.Sigma_Num_loose = None

        self.RR_start = None  # Швидкість репродукції на 1-ій ітерації
        self.RR_fin = None  # Швидкість репродукції на останній ітерації

        self.Teta_start = None  # Втрата різноманітності на 1-ій ітерації
        self.Teta_fin = None  # Втрата різноманітності на останній ітерації

        self.unique_X_start = None  # Кількість унікальних хромосом в початковій популяції
        self.unique_X_fin = None  # Кількість унікальних хромосом в фінальній популяції

        # Поле для збереження різниці відбору після першого відбору
        self.s_start = None

        # Нові поля для критеріїв PFET (Fish)
        self.Fish_start = None
        self.Fish_min = None
        self.NI_Fish_min = None
        self.Fish_max = None
        self.NI_Fish_max = None
        self.Fish_avg = None

        # Нові поля для критеріїв Pτ (Kend)
        self.Kend_start = None
        self.Kend_min = None
        self.NI_Kend_min = None
        self.Kend_max = None
        self.NI_Kend_max = None
        self.Kend_avg = None

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

        # Оновлення RR_start, RR_fin
        if gen_i == 1:
            self.RR_start = gen_stats.reproduction_rate
        self.RR_fin = gen_stats.reproduction_rate

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

        # Оновлення Teta_start, Teta_fin
        if gen_i == 1:
            self.Teta_start = gen_stats.loss_of_diversity
        self.Teta_fin = gen_stats.loss_of_diversity

            # Оновлення NI_loose та Num_loose
        if gen_stats.optimal_individual_lost:
            # Якщо оптимальна особина була втрачена в цій генерації
            self.Num_loose += 1  # Збільшуємо кількість втрат оптимальної особини
            self.NI_loose = gen_i  # Оновлюємо номер ітерації, коли втрачена оптимальна особина

            # Оновлення unique_X_start, unique_X_fin
        if gen_i == 1:
            self.unique_X_start = gen_stats.population.count_unique_chromosomes()
        self.unique_X_fin = gen_stats.population.count_unique_chromosomes()

        # Різниця відбору після першого відбору (s_start)
        if gen_i == 1:
            # Обчислюємо різницю відбору s_start
            self.s_start = gen_stats.difference

        # Оновлення Fish_start, Fish_min, Fish_max, Fish_avg
        if gen_i == 1:
            self.Fish_start = gen_stats.fish_value  # Початкове значення Fish
        if self.Fish_min is None or gen_stats.fish_value < self.Fish_min:
            self.Fish_min = gen_stats.fish_value
            self.NI_Fish_min = gen_i
        if self.Fish_max is None or gen_stats.fish_value > self.Fish_max:
            self.Fish_max = gen_stats.fish_value
            self.NI_Fish_max = gen_i
        if self.Fish_avg is None:
            self.Fish_avg = gen_stats.fish_value
        else:
            self.Fish_avg = (self.Fish_avg * (gen_i - 1) + gen_stats.fish_value) / gen_i

        # Оновлення Kend_start, Kend_min, Kend_max, Kend_avg
        if gen_i == 1:
            self.Kend_start = gen_stats.kend_value  # Початкове значення Kend
        if self.Kend_min is None or gen_stats.kend_value < self.Kend_min:
            self.Kend_min = gen_stats.kend_value
            self.NI_Kend_min = gen_i
        if self.Kend_max is None or gen_stats.kend_value > self.Kend_max:
            self.Kend_max = gen_stats.kend_value
            self.NI_Kend_max = gen_i
        if self.Kend_avg is None:
            self.Kend_avg = gen_stats.kend_value
        else:
            self.Kend_avg = (self.Kend_avg * (gen_i - 1) + gen_stats.kend_value) / gen_i

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

    def update_final_stats(self, gen_stats: GenerationStats, gen_i):
        # Оновлення unique_X_fin
        self.unique_X_fin = gen_stats.population.count_unique_chromosomes()  # Оновлення кількості унікальних хромосом

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
