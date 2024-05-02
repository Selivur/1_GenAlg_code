from config import NR
from stats.run_stats import RunStats
import numpy as np

class ExperimentStats:
    def __init__(self, experiment_params: tuple[str]):
        self.params = experiment_params
        self.runs = np.empty(NR, dtype=object)

        self.Suc = 0
        self.N_Suc = 0
        self.Min_NI = None
        self.Max_NI = None
        self.Avg_NI = None
        self.Sigma_NI = None

        # Reproduction Rate
        self.Min_RR_min = None
        self.NI_RR_min = None
        self.Max_RR_max = None
        self.NI_RR_max = None
        self.Avg_RR_min = None
        self.Avg_RR_max = None
        self.Avg_RR_avg = None
        self.Sigma_RR_min = None
        self.Sigma_RR_max = None
        self.Sigma_RR_avg = None

        # Loss of Diversity
        self.Min_Teta_min = None
        self.NI_Teta_min = None
        self.Max_Teta_max = None
        self.NI_Teta_max = None
        self.Avg_Teta_min = None
        self.Avg_Teta_max = None
        self.Avg_Teta_avg = None
        self.Sigma_Teta_min = None
        self.Sigma_Teta_max = None
        self.Sigma_Teta_avg = None

        # Selection Intensity
        self.Min_I_min = None
        self.NI_I_min = None
        self.Max_I_max = None
        self.NI_I_max = None
        self.Avg_I_min = None
        self.Avg_I_max = None
        self.Avg_I_avg = None
        self.Sigma_I_min = None
        self.Sigma_I_max = None
        self.Sigma_I_avg = None

        # Selection Difference
        self.Min_s_min = None
        self.NI_s_min = None
        self.Max_s_max = None
        self.NI_s_max = None
        self.Avg_s_min = None
        self.Avg_s_max = None
        self.Avg_s_avg = None

        # Growth Rate
        self.Min_GR_early = None
        self.Max_GR_early = None
        self.Avg_GR_early = None
        self.Min_GR_late = None
        self.Max_GR_late = None
        self.Avg_GR_late = None
        self.Min_GR_avg = None
        self.Max_GR_avg = None
        self.Avg_GR_avg = None

        # New
        self.nonSuc = None
        self.nonMin_NI = None
        self.nonMax_NI = None
        self.nonAvg_NI = None
        self.nonSigma_NI = None
        self.nonAvg_F_found = None
        self.nonSigma_F_found = None
        self.nonMax_F_found = None

        # Нові поля для зберігання критеріїв
        self.Min_I_start = None
        self.Max_I_start = None
        self.Avg_I_start = None
        self.Sigma_I_start = None

        self.Min_GR_start = None
        self.Max_GR_start = None
        self.Avg_GR_start = None
        self.Sigma_GR_start = None

        self.Min_Pr_min = None
        self.NI_Pr_min = None
        self.Max_Pr_max = None
        self.NI_Pr_max = None
        self.Avg_Pr_min = None
        self.Avg_Pr_max = None
        self.Avg_Pr_avg = None
        self.Sigma_Pr_min = None
        self.Sigma_Pr_max = None
        self.Sigma_Pr_avg = None

        self.Min_Pr_start = None
        self.Max_Pr_start = None
        self.Avg_Pr_start = None
        self.Sigma_Pr_start = None

        # Нові поля для зберігання критеріїв
        self.Avg_NI_loose = None
        self.Sigma_NI_loose = None
        self.Avg_Num_loose = None
        self.Sigma_Num_loose = None

        # Reproduction Rate
        self.Avg_RR_start = None
        self.Avg_RR_fin = None
        self.Min_RR_start = None
        self.Max_RR_start = None
        self.Sigma_RR_start = None
        self.Sigma_RR_fin = None

        # Loss of Diversity
        self.Avg_Teta_start = None
        self.Avg_Teta_fin = None
        self.Min_Teta_start = None
        self.Max_Teta_start = None
        self.Sigma_Teta_start = None
        self.Sigma_Teta_fin = None

        # Unique Chromosomes
        self.Avg_unique_X_start = None
        self.Avg_unique_X_fin = None
        self.Min_unique_X_start = None
        self.Max_unique_X_start = None
        self.Sigma_unique_X_start = None
        self.Sigma_unique_X_fin = None

        self.Min_s_start = None
        self.Max_s_start = None
        self.Avg_s_start = None
        self.Sigma_s_start = None

        # Додайте нові поля для критеріїв Фішера (Fish)
        self.Min_Fish_min = None
        self.NI_Fish_min = None
        self.Max_Fish_max = None
        self.NI_Fish_max = None
        self.Avg_Fish_min = None
        self.Avg_Fish_max = None
        self.Avg_Fish_avg = None
        self.Sigma_Fish_max_ = None
        self.Sigma_Fish_min_ = None
        self.Sigma_Fish_avg_ = None
        self.Min_Fish_start = None
        self.Max_Fish_start = None
        self.Avg_Fish_start = None
        self.Sigma_Fish_start = None

        # Додайте нові поля для критеріїв Кендела (Kend)
        self.Min_Kend_min = None
        self.NI_Kend_min = None
        self.Max_Kend_max = None
        self.NI_Kend_max = None
        self.Avg_Kend_min = None
        self.Avg_Kend_max = None
        self.Avg_Kend_avg = None
        self.Sigma_Kend_max_ = None
        self.Sigma_Kend_min_ = None
        self.Sigma_Kend_avg_ = None
        self.Min_Kend_start = None
        self.Max_Kend_start = None
        self.Avg_Kend_start = None
        self.Sigma_Kend_start = None

    def add_run(self, run: RunStats, run_i):
        self.runs[run_i] = run

    def calculate(self):
        successful_runs = [run for run in self.runs if run.is_successful]
        unsuccessful_runs = [run for run in self.runs if not run.is_successful]


        self.N_Suc = len(successful_runs)
        self.Suc = self.N_Suc / NR

        self.__calculate_convergence_stats(successful_runs)
        self.__calculate_rr_stats(successful_runs)
        self.__calculate_teta_stats(successful_runs)
        self.calculate_criteria(successful_runs)
        self.__calculate_unsuccessful_run_stats(unsuccessful_runs)
        self.__calculate_loose_stats(unsuccessful_runs)

        if self.params[0] != 'FconstALL':
            self.__calculate_s_stats(successful_runs)
            self.__calculate_i_stats(successful_runs)
            self.__calculate_gr_stats(successful_runs)

    def calculate_criteria(self, runs: list[RunStats]):
        successful_runs = [run for run in self.runs if run.is_successful]

        # Обчислення критеріїв для I_start
        I_start_values = [run.I_start for run in successful_runs if run.I_start is not None]
        if I_start_values:
            self.Min_I_start = min(I_start_values)
            self.Max_I_start = max(I_start_values)
            self.Avg_I_start = np.mean(I_start_values)
            self.Sigma_I_start = np.std(I_start_values)

        # Обчислення критеріїв для GR_start
        GR_start_values = [run.GR_start for run in successful_runs if run.GR_start is not None]
        if GR_start_values:
            self.Min_GR_start = min(GR_start_values)
            self.Max_GR_start = max(GR_start_values)
            self.Avg_GR_start = np.mean(GR_start_values)
            self.Sigma_GR_start = np.std(GR_start_values)

        # Обчислення критеріїв для Pr_min, Pr_max, Pr_avg
        Pr_min_values = [run.Pr_min for run in successful_runs if run.Pr_min is not None]
        Pr_max_values = [run.Pr_max for run in successful_runs if run.Pr_max is not None]
        Pr_avg_values = [run.Pr_avg for run in successful_runs if run.Pr_avg is not None]
        if Pr_min_values:
            self.Min_Pr_min = min(Pr_min_values)
            self.NI_Pr_min = successful_runs[Pr_min_values.index(self.Min_Pr_min)].NI_Pr_min
            self.Max_Pr_max = max(Pr_max_values)
            self.NI_Pr_max = successful_runs[Pr_max_values.index(self.Max_Pr_max)].NI_Pr_max
            self.Avg_Pr_min = np.mean(Pr_min_values)
            self.Avg_Pr_max = np.mean(Pr_max_values)
            self.Avg_Pr_avg = np.mean(Pr_avg_values)
            self.Sigma_Pr_min = np.std(Pr_min_values)
            self.Sigma_Pr_max = np.std(Pr_max_values)
            self.Sigma_Pr_avg = np.std(Pr_avg_values)

        # Обчислення критеріїв для Pr_start
        Pr_start_values = [run.Pr_start for run in successful_runs if run.Pr_start is not None]
        if Pr_start_values:
            self.Min_Pr_start = min(Pr_start_values)
            self.Max_Pr_start = max(Pr_start_values)
            self.Avg_Pr_start = np.mean(Pr_start_values)
            self.Sigma_Pr_start = np.std(Pr_start_values)

        # Обчислення RR_start, RR_fin, Sigma_RR_start, Sigma_RR_fin, Min_RR_start, Max_RR_start
        RR_start_values = [run.RR_start for run in successful_runs if run.RR_start is not None]
        if RR_start_values:
            self.Min_RR_start = min(RR_start_values)
            self.Max_RR_start = max(RR_start_values)
            self.Avg_RR_start = np.mean(RR_start_values)
            self.Sigma_RR_start = np.std(RR_start_values)

        RR_fin_values = [run.RR_fin for run in successful_runs if run.RR_fin is not None]
        if RR_fin_values:
            self.Avg_RR_fin = np.mean(RR_fin_values)
            self.Sigma_RR_fin = np.std(RR_fin_values)

        # Обчислення Teta_start, Teta_fin, Sigma_Teta_start, Sigma_Teta_fin, Min_Teta_start, Max_Teta_start
        Teta_start_values = [run.Teta_start for run in successful_runs if run.Teta_start is not None]
        if Teta_start_values:
            self.Min_Teta_start = min(Teta_start_values)
            self.Max_Teta_start = max(Teta_start_values)
            self.Avg_Teta_start = np.mean(Teta_start_values)
            self.Sigma_Teta_start = np.std(Teta_start_values)

        Teta_fin_values = [run.Teta_fin for run in successful_runs if run.Teta_fin is not None]
        if Teta_fin_values:
            self.Avg_Teta_fin = np.mean(Teta_fin_values)
            self.Sigma_Teta_fin = np.std(Teta_fin_values)

        # Обчислення Avg_unique_X_start, Avg_unique_X_fin, Sigma_unique_X_start, Sigma_unique_X_fin, Min_unique_X_start, Max_unique_X_start
        unique_X_start_values = [run.unique_X_start for run in successful_runs if run.unique_X_start is not None]
        if unique_X_start_values:
            self.Min_unique_X_start = min(unique_X_start_values)
            self.Max_unique_X_start = max(unique_X_start_values)
            self.Avg_unique_X_start = np.mean(unique_X_start_values)
            self.Sigma_unique_X_start = np.std(unique_X_start_values)

        unique_X_fin_values = [run.unique_X_fin for run in successful_runs if run.unique_X_fin is not None]
        if unique_X_fin_values:
            self.Avg_unique_X_fin = np.mean(unique_X_fin_values)
            self.Sigma_unique_X_fin = np.std(unique_X_fin_values)

        # Обчислення критеріїв для s_start
        s_start_values = [run.s_start for run in successful_runs if run.s_start is not None]
        if s_start_values:
            # Мінімальне значення s_start
            self.Min_s_start = min(s_start_values)

            # Максимальне значення s_start
            self.Max_s_start = max(s_start_values)

            # Середнє значення s_start
            self.Avg_s_start = np.mean(s_start_values)

            # Стандартне відхилення s_start
            self.Sigma_s_start = np.std(s_start_values)

        # Обчислення статистичних характеристик для Фішера (Fish)
        fish_min_values = [run.Fish_min for run in successful_runs if run.Fish_min is not None]
        fish_max_values = [run.Fish_max for run in successful_runs if run.Fish_max is not None]
        fish_start_values = [run.Fish_start for run in successful_runs if run.Fish_start is not None]
        fish_avg_values = [run.Fish_avg for run in successful_runs if run.Fish_avg is not None]

        if fish_min_values:
            self.Min_Fish_min = min(fish_min_values)
            self.NI_Fish_min = successful_runs[fish_min_values.index(self.Min_Fish_min)].NI_Fish_min

        if fish_max_values:
            self.Max_Fish_max = max(fish_max_values)
            self.NI_Fish_max = successful_runs[fish_max_values.index(self.Max_Fish_max)].NI_Fish_max

        if fish_start_values:
            self.Min_Fish_start = min(fish_start_values)
            self.Max_Fish_start = max(fish_start_values)
            self.Avg_Fish_start = np.mean(fish_start_values)
            self.Sigma_Fish_start = np.std(fish_start_values)

        if fish_avg_values:
            self.Avg_Fish_avg = np.mean(fish_avg_values)
            self.Sigma_Fish_avg_ = np.std(fish_avg_values)

            # Обчислення статистичних характеристик для Кендела (Kend)
        kend_min_values = [run.Kend_min for run in successful_runs if run.Kend_min is not None]
        kend_max_values = [run.Kend_max for run in successful_runs if run.Kend_max is not None]
        kend_start_values = [run.Kend_start for run in successful_runs if run.Kend_start is not None]
        kend_avg_values = [run.Kend_avg for run in successful_runs if run.Kend_avg is not None]

        if kend_min_values:
            self.Min_Kend_min = min(kend_min_values)
            self.NI_Kend_min = successful_runs[kend_min_values.index(self.Min_Kend_min)].NI_Kend_min

        if kend_max_values:
            self.Max_Kend_max = max(kend_max_values)
            self.NI_Kend_max = successful_runs[kend_max_values.index(self.Max_Kend_max)].NI_Kend_max

        if kend_start_values:
            self.Min_Kend_start = min(kend_start_values)
            self.Max_Kend_start = max(kend_start_values)
            self.Avg_Kend_start = np.mean(kend_start_values)
            self.Sigma_Kend_start = np.std(kend_start_values)

        if kend_avg_values:
            self.Avg_Kend_avg = np.mean(kend_avg_values)
            self.Sigma_Kend_avg_ = np.std(kend_avg_values)

    def __calculate_unsuccessful_run_stats(self, runs: list[RunStats]):
        nonSucs = len(runs) / NR
        self.nonSuc = nonSucs

        if runs:
            nonNIs = [run.NI for run in runs]
            self.nonMin_NI = min(nonNIs)
            self.nonMax_NI = max(nonNIs)
            self.nonAvg_NI = np.mean(nonNIs)
            self.nonSigma_NI = np.std(nonNIs)

            nonFounds = [run.F_found for run in runs if run.F_found is not None]
            if nonFounds:
                self.nonAvg_F_found = np.mean(nonFounds)
                self.nonSigma_F_found = np.std(nonFounds)
                self.nonMax_F_found = max(nonFounds)

    def __calculate_convergence_stats(self, runs: list[RunStats]):
        NIs = [run.NI for run in runs]
        if NIs:
            self.Min_NI = min(NIs)
            self.Max_NI = max(NIs)
            self.Avg_NI = np.mean(NIs)
            self.Sigma_NI = np.std(NIs)

    def __calculate_rr_stats(self, runs: list[RunStats]):
        RR_min_list = [run.RR_min for run in runs]
        if RR_min_list:
            run_i_RR_min = np.argmin(RR_min_list)
            self.NI_RR_min = runs[run_i_RR_min].NI_RR_min
            self.Min_RR_min = RR_min_list[run_i_RR_min]
            self.Avg_RR_min = np.mean(RR_min_list)
            self.Sigma_RR_min = np.std(RR_min_list)
        RR_max_list = [run.RR_max for run in runs]
        if RR_max_list:
            run_i_RR_max = np.argmax(RR_max_list)
            self.NI_RR_max = runs[run_i_RR_max].NI_RR_max
            self.Max_RR_max = RR_max_list[run_i_RR_max]
            self.Avg_RR_max = np.mean(RR_max_list)
            self.Sigma_RR_max = np.std(RR_max_list)
        RR_avg_list = [run.RR_avg for run in runs]
        if RR_avg_list:
            self.Avg_RR_avg = np.mean(RR_avg_list)
            self.Sigma_RR_avg = np.std(RR_avg_list)

    def __calculate_teta_stats(self, runs: list[RunStats]):
        Teta_min_list = [run.Teta_min for run in runs]
        if Teta_min_list:
            run_i_Teta_min = np.argmin(Teta_min_list)
            self.NI_Teta_min = runs[run_i_Teta_min].NI_Teta_min
            self.Min_Teta_min = Teta_min_list[run_i_Teta_min]
            self.Avg_Teta_min = np.mean(Teta_min_list)
            self.Sigma_Teta_min = np.std(Teta_min_list)
        Teta_max_list = [run.Teta_max for run in runs]
        if Teta_max_list:
            run_i_Teta_max = np.argmax(Teta_max_list)
            self.NI_Teta_max = runs[run_i_Teta_max].NI_Teta_max
            self.Max_Teta_max = Teta_max_list[run_i_Teta_max]
            self.Avg_Teta_max = np.mean(Teta_max_list)
            self.Sigma_Teta_max = np.std(Teta_max_list)
        Teta_avg_list = [run.Teta_avg for run in runs]
        if Teta_avg_list:
            self.Avg_Teta_avg = np.mean(Teta_avg_list)
            self.Sigma_Teta_avg = np.std(Teta_avg_list)

    def __calculate_s_stats(self, runs: list[RunStats]):
        s_min_list = [run.s_min for run in runs]
        if s_min_list:
            run_i_s_min = np.argmin(s_min_list)
            self.NI_s_min = runs[run_i_s_min].NI_s_min
            self.Min_s_min = s_min_list[run_i_s_min]
            self.Avg_s_min = np.mean(s_min_list)
        s_max_list = [run.s_max for run in runs]
        if s_max_list:
            run_i_s_max = np.argmax(s_max_list)
            self.NI_s_max = runs[run_i_s_max].NI_s_max
            self.Max_s_max = s_max_list[run_i_s_max]
            self.Avg_s_max = np.mean(s_max_list)
        s_avg_list = [run.s_avg for run in runs]
        if s_avg_list:
            self.Avg_s_avg = np.mean(s_avg_list)

    def __calculate_i_stats(self, runs: list[RunStats]):
        I_min_list = [run.I_min for run in runs]
        if I_min_list:
            run_i_I_min = np.argmin(I_min_list)
            self.NI_I_min = runs[run_i_I_min].NI_I_min
            self.Min_I_min = I_min_list[run_i_I_min]
            self.Avg_I_min = np.mean(I_min_list)
            self.Sigma_I_min = np.std(I_min_list)
        I_max_list = [run.I_max for run in runs]
        if I_max_list:
            run_i_I_max = np.argmax(I_max_list)
            self.NI_I_max = runs[run_i_I_max].NI_I_max
            self.Max_I_max = I_max_list[run_i_I_max]
            self.Avg_I_max = np.mean(I_max_list)
            self.Sigma_I_max = np.std(I_max_list)
        I_avg_list = [run.I_avg for run in runs]
        if I_avg_list:
            self.Avg_I_avg = np.mean(I_avg_list)
            self.Sigma_I_avg = np.std(I_avg_list)

    def __calculate_gr_stats(self, runs: list[RunStats]):
        gre_list = [run.GR_early for run in runs if run.GR_early is not None]
        grl_list = [run.GR_late for run in runs if run.GR_late is not None]
        gra_list = [run.GR_avg for run in runs if run.GR_avg is not None]
        if gre_list:
            self.Avg_GR_early = np.mean(gre_list)
            self.Min_GR_early = min(gre_list)
            self.Max_GR_early = max(gre_list)
        if grl_list:
            self.Avg_GR_late = np.mean(grl_list)
            self.Min_GR_late = min(grl_list)
            self.Max_GR_late = max(grl_list)
        if gra_list:
            self.Avg_GR_avg = np.mean(gra_list)
            self.Min_GR_avg = min(gra_list)
            self.Max_GR_avg = max(gra_list)

    def __calculate_loose_stats(self, runs: list[RunStats]):
        # Фільтрація прогонів з втратами оптимальної особини
        loose_runs = [run for run in runs if run.Num_loose > 0]

        if loose_runs:
            # Обчислення NI_loose_list - списку номерів ітерацій, на яких була втрачена оптимальна особина
            NI_loose_list = [run.NI_loose for run in loose_runs]

            # Обчислення середнього та стандартного відхилення номеру ітерації, на якій була втрачена оптимальна особина
            self.Avg_NI_loose = np.mean(NI_loose_list)
            self.Sigma_NI_loose = np.std(NI_loose_list)

            # Обчислення Num_loose_list - списку загальної кількості втрат оптимальної особини
            Num_loose_list = [run.Num_loose for run in loose_runs]

            # Обчислення середнього та стандартного відхилення загальної кількості втрат оптимальної особини
            self.Avg_Num_loose = np.mean(Num_loose_list)
            self.Sigma_Num_loose = np.std(Num_loose_list)
        else:
            # Якщо немає прогонів з втратами, встановлюємо параметри в None
            self.Avg_NI_loose = None
            self.Sigma_NI_loose = None
            self.Avg_Num_loose = None
            self.Sigma_Num_loose = None

    def __str__(self):
        return ("Suc: " + str(self.Suc) + "%" +
                "\nMin: " + str(self.Min_NI) + "\nMax: " + str(self.Max_NI) + "\nAvg: " + str(self.Avg_NI))
