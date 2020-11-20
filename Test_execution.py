from Tape_EA_Wenzel import load_settings, delete_old_population_and_start_chromo, copy_loaded_settings_to_main_folder, GUI_select_folder_directory

def load_setup(setup_subdir):
    delete_old_population_and_start_chromo()
    copy_loaded_settings_to_main_folder(setup_subdir)  # Copy the loaded settingssheet and startparameter into the current folder to work on it, without changing the initial settings
    try:
        load_settings("./")  # Since we copy the settings into our main folder, we can load the settings directly from here.
    except:
        print("Settings are missing - Check folder")
        exit()
def append_best_fitness_of_run(adap_mutation, border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, mutation_rate_list, num_gen_list,penalty_alpha_list, penalty_beta_list, penalty_length_list,penalty_minimise_beta_list, penalty_negativ_list):
    global fitness_list_execution, distance_fit_list_execution, length_fit_list_execution, border_fit_start_list_execution, border_fit_end_list_execution, mutation_rate_list_execution, penalty_alpha_list_execution, penalty_beta_list_execution, penalty_length_list_execution, penalty_negativ_list_execution, penalty_minimise_beta_list_execution
    fitness_list_execution.append(fitness_list[-1])  # Comment_DB: stack gen_list with list of fitnesses
    distance_fit_list_execution.append(distance_fit_list[-1])
    length_fit_list_execution.append(length_fit_list[-1])
    border_fit_start_list_execution.append(border_fit_start_list[-1])
    border_fit_end_list_execution.append(border_fit_end_list[-1])
    if adap_mutation == 1:
        mutation_rate_list_execution.append(mutation_rate_list[-1])
    penalty_alpha_list_execution.append(penalty_alpha_list[-1])
    penalty_beta_list_execution.append(penalty_beta_list[-1])
    penalty_length_list_execution.append(penalty_length_list[-1])
    penalty_negativ_list_execution.append(penalty_negativ_list[-1])
    penalty_minimise_beta_list_execution.append(penalty_minimise_beta_list[-1])
def save_fitness_per_setup(setup_subdir,setup_nr):
    file2 = open(setup_subdir+"fitness_for_setup"+str(setup_nr)+".txt", "w")
    list_gen_index = [fitness_list_execution, distance_fit_list_execution, length_fit_list_execution, border_fit_start_list_execution, border_fit_end_list_execution, mutation_rate_list_execution, penalty_alpha_list_execution, penalty_beta_list_execution, penalty_length_list_execution, penalty_negativ_list_execution, penalty_minimise_beta_list_execution]

    for list_1 in list_gen_index:
        file2.write(str(list_1)+ "\n")

    file2.close()
def save_results_of_run():
    create_subdir_with_results(setup_subdir)
    from Tape_EA_Wenzel import adap_mutation, border_fit_end_list, border_fit_start_list, distance_fit_list, \
        fitness_list, length_fit_list, mutation_rate_list, num_gen_list, penalty_alpha_list, penalty_beta_list, \
        penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list
    append_best_fitness_of_run(adap_mutation, border_fit_end_list, border_fit_start_list, distance_fit_list,
                               fitness_list, length_fit_list, mutation_rate_list, num_gen_list,
                               penalty_alpha_list, penalty_beta_list, penalty_length_list,
                               penalty_minimise_beta_list, penalty_negativ_list)
############################## GUI and read in start chromosoms ####################################################################

# Select Folder of test run
sub_dir = GUI_select_folder_directory()


different_setups = 2
amount_of_execution_per_setup = 2
execute_main_EA = True
execute_individual = False
execute_Population_EA = False
use_best_Solution = False

for i in range(different_setups):
    setup_subdir = sub_dir + "Setup" + str(i) + "/"
    load_setup(setup_subdir)
    from Tape_EA_Wenzel import load_preprocessor_start_chromosomes, calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution,calc_best_Solution  # import loaded settings

    # Read in Startchromosoms from sub_dir
    load_preprocessor_start_chromosomes(calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution, calc_best_Solution, setup_subdir)
    from Tape_EA_Wenzel import *

    # Initialize lists to save result of setup run
    fitness_list_execution, distance_fit_list_execution, length_fit_list_execution, border_fit_start_list_execution, border_fit_end_list_execution, mutation_rate_list_execution, penalty_alpha_list_execution, penalty_beta_list_execution, penalty_length_list_execution, penalty_negativ_list_execution, penalty_minimise_beta_list_execution = [], [], [], [], [], [], [], [], [], [], []

    for j in range(amount_of_execution_per_setup):


        if execute_individual:
            print("Individual Optimization")
            individual_optimization_and_removing_unnecessary_bends(use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution, use_best_Solution)
            save_fitness_over_gen()
            save_results_of_run()

        if execute_main_EA:
            print("Main EA Optimization")
            repair_start_chromosomes_to_same_amount_of_bendpoints(use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution, use_best_Solution)
            main_EA(adap_mutation, amount_of_bends, num_gen, p_mutation, pop_size, startchromo2D, startchromo2D_edge,
            startchromo3D,startchromo_current_best, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution, use_best_Solution,
                    False)

            save_results_of_run()

        if execute_Population_EA:
            input_pop1,input_pop2,input_pop3 = [],[],[]
            input_pop1 = setup_subdir+"population_main_Start.txt"

            if input_pop1 != [] or input_pop2 != [] or input_pop3 != []:
                main_EA_with_loaded_pop(adap_mutation, amount_of_bends, input_pop1, input_pop2, input_pop3, num_gen,
                                        p_mutation, pop_size, False)
                save_results_of_run()
    save_fitness_per_setup(setup_subdir, i)



