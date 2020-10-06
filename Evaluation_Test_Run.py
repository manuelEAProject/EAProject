from Tape_EA_Wenzel import load_settings, delete_old_population_and_start_chromo, copy_loaded_settings_to_main_folder, GUI_select_folder_directory
import matplotlib.pyplot as plt
import os
import numpy as np

def get_results_of_test_run(sub_dir):
    list_setup_results = []

    for setup in os.listdir(sub_dir):
        setup_subdir = sub_dir + setup + "/"
        list_setup_executions = []

        for name in os.listdir(setup_subdir):
            if os.path.isdir(setup_subdir + name):
                file = setup_subdir + name + "/fitness_over_generation.txt"
                file2 = open(file, "r")
                generations_list_execution, fitness_list_execution, distance_fit_list_execution, length_fit_list_execution, border_fit_start_list_execution, border_fit_end_list_execution, mutation_rate_list_execution, penalty_alpha_list_execution, penalty_beta_list_execution, penalty_length_list_execution, penalty_negativ_list_execution, penalty_minimise_beta_list_execution = [], [], [], [], [], [], [], [], [], [], [], []
                lists_execution = [generations_list_execution, fitness_list_execution, distance_fit_list_execution,
                                   length_fit_list_execution, border_fit_start_list_execution,
                                   border_fit_end_list_execution, mutation_rate_list_execution,
                                   penalty_alpha_list_execution, penalty_beta_list_execution,
                                   penalty_length_list_execution,
                                   penalty_negativ_list_execution, penalty_minimise_beta_list_execution]

                for i in range(12):

                    list1 = file2.readline()[1:]
                    while list1[-2] != "]":
                        list1[:-2]
                        list1 = list1 + file2.readline()

                    list1 = list1[:-2]

                    list1 = list1.split(" ")
                    list2 = []
                    for j in range(len(list1)):
                        list1[j].rsplit("\n")
                        if list1[j] == "":
                            pass
                        else:
                            list2.append(list1[j])
                    list2 = list(map(float, list2))

                    lists_execution[i] = list2
                list_setup_executions.append(lists_execution)
        list_setup_results.append(list_setup_executions)
        # Load/Open Fitness Over Generation
    return list_setup_results
def plot_setup_results(j, list_setup_results,color1,color2):
    plt.figure(1)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][1], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][1], linewidth=1, color=color2)

    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.figure(2)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][2], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][2], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Distance-Fit')
    plt.figure(3)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][3], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][3], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Length-Fit')
    plt.figure(4)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][4], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][4], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Border-Fit-Start')
    plt.figure(5)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][5], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][5], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Border-Fit-End')
    plt.figure(6)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][6], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][6], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Mutation Rate')
    plt.figure(7)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][7], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][7], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Penalty - Alpha')
    plt.figure(8)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][8], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][8], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Penalty - Beta')
    plt.figure(9)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][9], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][9], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Penalty - Length')
    plt.figure(10)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][10], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][10], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Penalty Negative')
    plt.figure(11)
    for i in range(10):
        plt.plot(list_setup_results[j][i][0], list_setup_results[j][i][11], linewidth=1, color=color1)
        plt.plot(list_setup_results[j + 1][i][0], list_setup_results[j + 1][i][11], linewidth=1, color=color2)
    plt.xlabel('Generation')
    plt.ylabel('Penalty - minimize Beta')


############################## GUI and read in start chromosoms ####################################################################

# Select Folder of test run
sub_dir = GUI_select_folder_directory()
sub_dir2 = GUI_select_folder_directory()



list_setup_results = get_results_of_test_run(sub_dir)
list_setup_results2 = get_results_of_test_run(sub_dir2)

j = 0
plot_setup_results(j, list_setup_results,"r","b")
plt.show(block=False)
k = 0
plot_setup_results(k, list_setup_results2,"y","green")
plt.show()

max_fitness_end = []
max_fitness_end1 = []
max_fitness_end2 = []
max_fitness_end12 = []
for i in range(len(list_setup_results[j])):
    max_fitness_end.append(list_setup_results[j][i][1][-1])
    max_fitness_end1.append(list_setup_results[j+1][i][1][-1])
    max_fitness_end2.append(list_setup_results2[k][i][1][-1])
    max_fitness_end12.append(list_setup_results2[k + 1][i][1][-1])




plt.figure(12)
plt.boxplot([max_fitness_end,max_fitness_end1])#,max_fitness_end2,max_fitness_end12],positions=[1,3,2,4])

plt.show()
