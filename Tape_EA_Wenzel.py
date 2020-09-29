import numpy as np
import inspect
import os
from pyquaternion import Quaternion
import math
import trimesh
from trimesh import proximity
from galileo_EA_Wenzel import Population
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import stl_preprocessing_Wenzel as stlprep3_6
from tkinter import *
from tkinter import filedialog
from timeit import default_timer as timer
import datetime
import time
import shutil


# GUI-Settings
# Preprocessor GUI
def get_Preprocessor_Vars_from_GUI():
    master = Tk()
    master.protocol("WM_DELETE_WINDOW", sys.exit)
    Label(master, text="Settings for Preprocessor").grid(row=10, sticky=W)
    Label(master, justify=LEFT, text=" ").grid(row=11, sticky=W)

    input_file = Entry(master)
    if os.path.isfile('./settingssheet.txt'):
        settingssheet = open('./settingssheet.txt')
        input_file.insert(0, settingssheet.readline())
        settingssheet.close()
    else:
        input_file.insert(0, "Select file...")
    input_file.grid(row=20, columnspan=3, sticky=W + E + N + S)

    Button(text='Browse..', command=lambda: select_stl_file(input_file)).grid(row=20, column=3, sticky=W)
    Label(master, justify=LEFT, text=" ").grid(row=21, sticky=W)
    Label(master, text="Tape-Typ :").grid(row=22, sticky=W)
    tape_type = Entry(master)
    tape_type.insert(0, "BASF_CFK_Tape")
    tape_type.grid(row=22, column=1, sticky=W)
    Label(master, text="Tape width [mm]:").grid(row=23, sticky=W)
    width = Entry(master)
    width.insert(0, 15)
    width.grid(row=23, column=1, sticky=W)

    #calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution
    Label(master, text="Preprocessor start solutions:").grid(row=25, sticky=W)

    calc_2D_Solution = BooleanVar()
    calc_2D_Solution.set(True)
    calc_2D_Solution_cb = Checkbutton(master, text="2D",
                                                 variable=calc_2D_Solution)
    calc_2D_Solution_cb.grid(row=25, column=1, sticky=W)

    calc_2D_with_edge_detection = BooleanVar()
    calc_2D_with_edge_detection.set(True)
    calc_2D_with_edge_detection_cb = Checkbutton(master, text="2D with edge-detection", variable=calc_2D_with_edge_detection)
    calc_2D_with_edge_detection_cb.grid(row=26, column=1, sticky=W)

    calc_3D_Solution = BooleanVar()
    calc_3D_Solution.set(True)
    calc_3D_Solution_cb = Checkbutton(master, text="3D with edge-detection",
                                                 variable=calc_3D_Solution)
    calc_3D_Solution_cb.grid(row=27, column=1, sticky=W)


    Label(master, text=" ").grid(row=29, sticky=W)

    Label(master, text="Resolution for surface interpolation (grid)").grid(row=43, sticky=W)
    grid_resolution = Entry(master)
    grid_resolution.insert(0, 1000)
    grid_resolution.grid(row=43, column=1, sticky=W)
    Label(master, text="Width for edge-detection").grid(row=44, sticky=W)
    width_for_edge_detection = Entry(master)
    width_for_edge_detection.insert(0, 2)
    width_for_edge_detection.grid(row=44, column=1, sticky=W)


    Label(master, text='Maximum distance from lin. fitted curve to surface').grid(row=45, sticky=W)
    max_distance = Entry(master)
    max_distance.insert(0, 1)  # Comment_DB: insert 10 at front of list
    max_distance.grid(row=45, column=1, sticky=W)
    Label(master, justify=LEFT, text=" ").grid(row=46, sticky=W)



    Button(master, text='Abort', command=sys.exit).grid(row=1000, column=0, pady=4)
    Button(master, text='Start', command=master.quit).grid(row=1000, column=3, pady=4)


    use_last_setting = BooleanVar()
    use_last_setting_cb = Checkbutton(master, text="continue with current results", command= lambda :load_setting.set(False),
                                       variable=use_last_setting)
    use_last_setting_cb.grid(row=1000, column=2, pady=4)
    use_last_setting.set(False)

    load_setting = BooleanVar()
    load_setting_cb = Checkbutton(master, text="load settings", command= lambda :use_last_setting.set(False),
                                      variable=load_setting)
    load_setting_cb.grid(row=1000, column=1, pady=4)
    load_setting.set(False)



    if_settingssheet_exists_fill_values(input_file, max_distance, grid_resolution, width, tape_type, width_for_edge_detection, calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution)

    mainloop()  # Executes GUI

    use_last_setting = bool(use_last_setting.get())
    load_setting = bool(load_setting.get())


    input_file = input_file.get()
    tape_type = tape_type.get()
    width = float(width.get())  # Tape width

    # Settings for Preprocessor
    calc_2D_Solution = bool(calc_2D_Solution.get())
    calc_2D_with_edge_detection = bool(calc_2D_with_edge_detection.get())
    calc_3D_Solution = bool(calc_3D_Solution.get())

    grid_resolution = int(grid_resolution.get())
    width_for_edge_detection = float(width_for_edge_detection.get())
    max_distance = float(max_distance.get())



    if not use_last_setting and not load_setting:
        save_settings([input_file, width, tape_type,
                    grid_resolution,
                   width_for_edge_detection, max_distance, calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution],"settingssheet.txt")  # Saves the selected settings

    master.destroy()  # closing of settings window
    return use_last_setting, load_setting, tape_type, width, input_file, max_distance, grid_resolution, width_for_edge_detection, \
           calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution
def change_amount_of_bendpoints_GUI(max_distance_value,width_for_edge_detection_value):
    master = Tk()
    master.protocol("WM_DELETE_WINDOW", sys.exit)
    Label(master, text="Settings for Preprocessor").grid(row=10, sticky=W)
    Label(master, justify=LEFT, text=" ").grid(row=11, sticky=W)


    Label(master, text="Width for edge-detection").grid(row=44, sticky=W)
    width_for_edge_detection = Entry(master)
    width_for_edge_detection.insert(0,str(width_for_edge_detection_value))
    width_for_edge_detection.grid(row=44, column=1, sticky=W)


    Label(master, text='Maximum distance from lin. fitted curve to surface').grid(row=45, sticky=W)
    max_distance = Entry(master)
    max_distance.insert(0,str(max_distance_value))  # Comment_DB: insert 10 at front of list
    max_distance.grid(row=45, column=1, sticky=W)
    Label(master, justify=LEFT, text=" ").grid(row=46, sticky=W)



    Button(master, text='Abort', command=sys.exit).grid(row=1000, column=0, pady=4)
    Button(master, text='Start', command=master.quit).grid(row=1000, column=3, pady=4)


    continue_check = BooleanVar()
    continue_check_cb = Checkbutton(master, text="Try different setting:",
                                       variable=continue_check)
    continue_check_cb.grid(row=1000, column=2, pady=4)
    continue_check.set(False)


    mainloop()  # Executes GUI


    width_for_edge_detection = float(width_for_edge_detection.get())
    max_distance = float(max_distance.get())
    continue_check = bool(continue_check.get())

    master.destroy()  # closing of settings window
    return continue_check, max_distance, width_for_edge_detection
# EA GUI
def get_EA_Vars_from_GUI(calc_3D_Solution,calc_2D_Solution,calc_2D_with_edge_detection):
    master = Tk()
    master.protocol("WM_DELETE_WINDOW", sys.exit)
    Label(master, text="Settings for EA").grid(row=10, sticky=W)
    Label(master, text="").grid(row=11, sticky=W)
    Label(master, justify=LEFT, text="Use start chromosoms: ").grid(row=12)

    use_2D_Solution = BooleanVar()
    use_2D_Solution.set(calc_2D_Solution)
    if calc_2D_Solution:
        Label(master, text="Use 2D preprocessor start chromosom:").grid(row=14, sticky=W)
        use_2D_Solution_cb = Checkbutton(master, text="", variable=use_2D_Solution)
        use_2D_Solution_cb.grid(row=14, column=1, sticky=W)


    use_2D_with_edge_detection = BooleanVar()
    use_2D_with_edge_detection.set(calc_2D_with_edge_detection)
    if calc_2D_with_edge_detection:
        Label(master, text="Use 2D(edge detection)preprocessor start chromosom:").grid(row=16, sticky=W)
        use_2D_with_edge_detection_cb = Checkbutton(master, text="", variable=use_2D_with_edge_detection)
        use_2D_with_edge_detection_cb.grid(row=16, column=1, sticky=W)

    use_3D_Solution = BooleanVar()
    use_3D_Solution.set(calc_3D_Solution)
    if calc_3D_Solution:
        Label(master, text="Use 3D preprocessor start chromosom:").grid(row=18, sticky=W)
        use_3D_Solution_cb = Checkbutton(master, text="",variable=use_3D_Solution)
        use_3D_Solution_cb.grid(row=18, column=1, sticky=W)

    # Loading Option for whole population
    if (os.path.isfile('./population_main_Start.txt') or
        os.path.isfile('./population_3D_Start.txt') or
        os.path.isfile('./population_2D_Start.txt')or
        os.path.isfile('./population_2DE_Start.txt')):

        Label(master, text="OR").grid(row=12, column=1)
        Label(master, text="Use Population:").grid(row=12, column=2)
        Label(master, justify=LEFT, text=" ").grid(row=12, column=2, sticky=W)

        Label(master, text="Population 1:").grid(row=14, column=2, sticky=W)
        input_pop1 = Entry(master)
        input_pop1.insert(0, "Select file...")
        input_pop1.grid(row=14, column=3, sticky=W + E + N + S)
        Button(text='Browse..', command=lambda: select_populaton_file(input_pop1)).grid(row=14, column=4, sticky=W)

        Label(master, text="Population 2:").grid(row=16, column=2, sticky=W)
        input_pop2 = Entry(master)
        input_pop2.insert(0, "Select file...")
        input_pop2.grid(row=16, column=3, sticky=W + E + N + S)
        Button(text='Browse..', command=lambda: select_populaton_file(input_pop2)).grid(row=16, column=4, sticky=W)

        Label(master, text="Population 3:").grid(row=18, column=2, sticky=W)
        input_pop3 = Entry(master)
        input_pop3.insert(0, "Select file...")
        input_pop3.grid(row=18, column=3, sticky=W + E + N + S)
        Button(text='Browse..', command=lambda: select_populaton_file(input_pop3)).grid(row=18, column=4, sticky=W)

    Button(text='Show Solutions', command=lambda: show_startchromos_results()).grid(row=20, column=0,  pady=1)


    Label(master, justify=LEFT, text=" ").grid(row=20, sticky=W)

    Label(master, text="Individual optimization:").grid(row=25, sticky=W)
    individual_optimization = BooleanVar()
    individual_optimization.set(True)
    individual_optimization_cb = Checkbutton(master, text="",
                                      variable=individual_optimization)
    individual_optimization_cb.grid(row=25, column=1, sticky=W)



    Label(master, text="Turn around stl-normal:").grid(row=25, column=2, sticky=W)
    turn_around_normal = BooleanVar()
    turn_around_normal.set(False)
    turn_around_normal_cb = Checkbutton(master, text="",
                                             variable=turn_around_normal)
    turn_around_normal_cb.grid(row=25, column=3, sticky=W)


    Label(master, justify=LEFT, text=" ").grid(row=71, sticky=W)
    #####Comment_DB: Different weightings for 4 equally-divided Population Sets!#####
    Label(master, text="Fine tuning fitness function (4 Population Sets)").grid(row=73, sticky=W)
    Label(master, text="Weighting distance fitness [gamma_d]: ").grid(row=74, column=0, sticky=W)
    Label(master, text="Set 1").grid(row=73, column=1, sticky=W)
    gamma_d = Entry(master)
    gamma_d.insert(0, 5.8)
    gamma_d.grid(row=74, column=1, sticky=W)
    Label(master, text="Set 2").grid(row=73, column=2, sticky=W)
    gamma_d2 = Entry(master)
    gamma_d2.insert(0, 5.8)
    gamma_d2.grid(row=74, column=2, sticky=W)
    Label(master, text="Set 3").grid(row=73, column=3, sticky=W)
    gamma_d3 = Entry(master)
    gamma_d3.insert(0, 5.8)
    gamma_d3.grid(row=74, column=3, sticky=W)
    Label(master, text="Set 4").grid(row=73, column=4, sticky=W)
    gamma_d4 = Entry(master)
    gamma_d4.insert(0, 5.8)
    gamma_d4.grid(row=74, column=4, sticky=W)
    Label(master, text="Weighting length fitness [gamma_l] : ").grid(row=75, column=0, sticky=W)
    gamma_l = Entry(master)
    gamma_l.insert(0, 1.4)
    gamma_l.grid(row=75, column=1, sticky=W)
    gamma_l2 = Entry(master)
    gamma_l2.insert(0, 1.4)
    gamma_l2.grid(row=75, column=2, sticky=W)
    gamma_l3 = Entry(master)
    gamma_l3.insert(0, 1.4)
    gamma_l3.grid(row=75, column=3, sticky=W)
    gamma_l4 = Entry(master)
    gamma_l4.insert(0, 1.4)
    gamma_l4.grid(row=75, column=4, sticky=W)
    Label(master, text="Weighting start fitness [gamma_ps] : ").grid(row=76, column=0, sticky=W)
    gamma_ps = Entry(master)
    gamma_ps.insert(0, 1)
    gamma_ps.grid(row=76, column=1, sticky=W)
    gamma_ps2 = Entry(master)
    gamma_ps2.insert(0, 1.4)
    gamma_ps2.grid(row=76, column=2, sticky=W)
    gamma_ps3 = Entry(master)
    gamma_ps3.insert(0, 1.4)
    gamma_ps3.grid(row=76, column=3, sticky=W)
    gamma_ps4 = Entry(master)
    gamma_ps4.insert(0, 1.4)
    gamma_ps4.grid(row=76, column=4, sticky=W)
    Label(master, text="Weighting end fitness [gamma_pe] : ").grid(row=77, column=0, sticky=W)
    gamma_pe = Entry(master)
    gamma_pe.insert(0, 1.4)
    gamma_pe.grid(row=77, column=1, sticky=W)
    gamma_pe2 = Entry(master)
    gamma_pe2.insert(0, 1.4)
    gamma_pe2.grid(row=77, column=2, sticky=W)
    gamma_pe3 = Entry(master)
    gamma_pe3.insert(0, 1.4)
    gamma_pe3.grid(row=77, column=3, sticky=W)
    gamma_pe4 = Entry(master)
    gamma_pe4.insert(0, 1.4)
    gamma_pe4.grid(row=77, column=4, sticky=W)
    Label(master, text="Point of Change (# of Generations) : ").grid(row=78, column=0, sticky=W)
    num_gen_set2 = Entry(master)
    num_gen_set2.insert(0, 5)
    num_gen_set2.grid(row=78, column=2)
    num_gen_set3 = Entry(master)
    num_gen_set3.insert(0, 15)
    num_gen_set3.grid(row=78, column=3)
    num_gen_set4 = Entry(master)
    num_gen_set4.insert(0, 30)
    num_gen_set4.grid(row=78, column=4)
    Label(master, justify=LEFT, text=" ").grid(row=79, sticky=W)

    Label(master, text="Resolution of bending/allel parameter:").grid(row=85, column=0, sticky=W)
    Label(master, justify=LEFT, text="Population size:").grid(row=90, column=0, sticky=W)
    pop_size = Entry(master)
    pop_size.insert(0, 12)
    pop_size.grid(row=90, column=1, sticky=W)


    Label(master, text="Number of generations:").grid(row=100, column=0, sticky=W)
    num_gen = Entry(master)
    num_gen.insert(0, 80)
    num_gen.grid(row=100, column=1, sticky=W)

    Label(master, text="Mutation Rate:").grid(row=120, sticky=W)
    p_mutation = Entry(master)
    p_mutation.insert(0, 0.1)
    p_mutation.grid(row=120, column=1, sticky=W)
    Label(master, text="Mutation range:").grid(row=120, column=2, sticky=W)
    p_mutate_range = Entry(master)
    p_mutate_range.insert(0, 0.05)
    p_mutate_range.grid(row=120, column=3, sticky=W)
    adap_mutation = BooleanVar()
    adap_mutation.set(False)
    adap_mutation_cb = Checkbutton(master, text="Adaptive Mutation", variable=adap_mutation)
    adap_mutation_cb.grid(row=130, column=3, sticky=W)
    Label(master, text="Crossoverrate:").grid(row=130, sticky=W)
    p_crossover = Entry(master)
    p_crossover.insert(0, 0.7)
    p_crossover.grid(row=130, column=1, sticky=W)

    Label(master, text="").grid(row=140, sticky=W)
    Label(master, text="Point resolution of calculated patch:").grid(row=150, sticky=W)
    fix_number_of_pts = BooleanVar()
    fix_number_of_pts_cb = Checkbutton(master, text="Fix number of points between bending points",
                                       variable=fix_number_of_pts,
                                       command=lambda: equidistant_pts_between_bendpts.set(False))
    fix_number_of_pts_cb.grid(row=160, column=0, sticky=W)
    fix_number_of_pts.set(True)
    Label(master, text="Number of Points:").grid(row=160, column=1, sticky=W)
    pointspersection = Entry(master)
    pointspersection.insert(0, 7)
    pointspersection.grid(row=160, column=2, sticky=W)
    equidistant_pts_between_bendpts = BooleanVar()
    equidistant_pts_between_bendpts_cb = Checkbutton(master, text="Equidistant Points",
                                                     variable=equidistant_pts_between_bendpts,
                                                     command=lambda: fix_number_of_pts.set(False),
                                                     disabledforeground="gray")
    equidistant_pts_between_bendpts_cb.grid(row=170, column=0, sticky=W)
    Label(master, text="Distance Points [mm]:").grid(row=170, column=1, sticky=W)
    step_size = Entry(master)
    step_size.insert(0, 6)
    step_size.grid(row=170, column=2, sticky=W)

    Button(master, text='Abort', command=sys.exit).grid(row=1000, column=0, pady=4)
    Button(master, text='Start', command=master.quit).grid(row=1000, column=1, pady=4)

    if_settingssheet_EA_exists_fill_values(individual_optimization,turn_around_normal, adap_mutation,
                                        gamma_d, gamma_d2, gamma_d3, gamma_d4, gamma_l, gamma_l2,
                                        gamma_l3, gamma_l4, gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, gamma_ps,
                                        gamma_ps2, gamma_ps3, gamma_ps4,
                                        num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_crossover,
                                        p_mutate_range, p_mutation, pop_size,
                                        pointspersection, equidistant_pts_between_bendpts, step_size)

    mainloop()  # Executes GUI

    use_2D_Solution = bool(use_2D_Solution.get())
    use_2D_with_edge_detection = bool(use_2D_with_edge_detection.get())
    use_3D_Solution = bool(use_3D_Solution.get())
    individual_optimization = bool(individual_optimization.get())
    turn_around_normal = bool(turn_around_normal.get())

    # Load population if wanted
    if (os.path.isfile('./population_main_Start.txt') or
        os.path.isfile('./population_3D_Start.txt') or
        os.path.isfile('./population_2D_Start.txt')or
        os.path.isfile('./population_2DE_Start.txt')):
        input_pop1 = str(input_pop1.get())
        input_pop2 = str(input_pop2.get())
        input_pop3 = str(input_pop3.get())
    try:
        if input_pop1 == "Select file...": input_pop1 = []
    except:  input_pop1 = []
    try:
        if input_pop2 == "Select file...": input_pop2 = []
    except:  input_pop2 = []
    try:
        if input_pop3 == "Select file...": input_pop3 = []
    except:  input_pop3 = []


    # Different gammas
    gamma_d = float(gamma_d.get())
    gamma_d2 = float(gamma_d2.get())
    gamma_d3 = float(gamma_d3.get())
    gamma_d4 = float(gamma_d4.get())
    gamma_l = float(gamma_l.get())
    gamma_l2 = float(gamma_l2.get())
    gamma_l3 = float(gamma_l3.get())
    gamma_l4 = float(gamma_l4.get())
    gamma_ps = float(gamma_ps.get())
    gamma_ps2 = float(gamma_ps2.get())
    gamma_ps3 = float(gamma_ps3.get())
    gamma_ps4 = float(gamma_ps4.get())
    gamma_pe = float(gamma_pe.get())
    gamma_pe2 = float(gamma_pe2.get())
    gamma_pe3 = float(gamma_pe3.get())
    gamma_pe4 = float(gamma_pe4.get())

    num_gen_set2 = float(num_gen_set2.get())
    num_gen_set3 = float(num_gen_set3.get())
    num_gen_set4 = float(num_gen_set4.get())


    # Settings for Evolutionary Algorithm
    num_gen = int(num_gen.get())  # Number of generations        Comment_DB: User inputs
    pop_size = int(pop_size.get())  # Population size


    # Crossover Probability     Comment_DB: Probability nature not defined yet
    p_crossover = float(p_crossover.get())
    # Mutations Probability
    p_mutation = float(p_mutation.get())
    # Adaptive Mutation
    adap_mutation = bool(adap_mutation.get())
    # Mutation range (in which percentage range may the allele mutate?)
    p_mutate_range = float(p_mutate_range.get())

    # Should the points between the bending points be evenly distributed (True,-> step_size) or a fixed number of points
    # Choose between the bending points (False,-> pointspersection)

    equidistant_pts_between_bendpts = bool(equidistant_pts_between_bendpts.get())
    fix_number_of_pts = bool(fix_number_of_pts.get())
    # distance between patch points - ONLY IF equidist_pts_between_bendpts=True
    step_size = float(step_size.get())
    # Number of points between the bending points
    pointspersection = int(pointspersection.get())


    save_settings([use_2D_with_edge_detection, use_2D_Solution, use_3D_Solution, individual_optimization, turn_around_normal, adap_mutation,
     gamma_d, gamma_d2,
     gamma_d3, gamma_d4, gamma_l, gamma_l2, gamma_l3, gamma_l4, gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, gamma_ps,
     gamma_ps2, gamma_ps3, gamma_ps4, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_crossover, p_mutate_range,
     p_mutation, pop_size, pointspersection, equidistant_pts_between_bendpts, step_size], "settingssheet_EA.txt")  # Saves the selected settings

    master.destroy()  # closing of settings window
    return use_2D_with_edge_detection, use_2D_Solution, use_3D_Solution, individual_optimization, turn_around_normal, adap_mutation, gamma_d, gamma_d2, gamma_d3, gamma_d4, gamma_l, gamma_l2, gamma_l3, gamma_l4, gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, gamma_ps,gamma_ps2, gamma_ps3, gamma_ps4, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_crossover, p_mutate_range,p_mutation, pop_size,pointspersection, equidistant_pts_between_bendpts, step_size, input_pop1, input_pop2, input_pop3

# Functions in GUI-Settings
def select_stl_file(input_file):
    if os.path.isfile('./settingssheet.txt'):
        settingssheet = open('./settingssheet.txt')
        name_of_stl_file = filedialog.askopenfilename(initialdir=settingssheet.readline(), title="Select file",
                                                      filetypes=(("stl files", "*.stl"), ("all files", "*.*")))
    else:
        name_of_stl_file = filedialog.askopenfilename(initialdir="/", title="Select file",
                                                      filetypes=(("stl files", "*.stl"), ("all files", "*.*")))
    input_file.delete(0, 'end')
    input_file.insert(0, name_of_stl_file)
    if len(input_file.get()) < 1:
        input_file.insert(0, "Select stl-file...")

    settingssheet = open('./settingssheet.txt', 'w+')
    settingssheet.write(input_file.get())
    settingssheet.close()
def select_populaton_file(input_file):
    name_of_pop_file = filedialog.askopenfilename(initialdir="./", title="Select file",
                                                      filetypes=(("txt files", "*.txt"), ("all files", "*.*")))
    input_file.delete(0, 'end')
    input_file.insert(0, name_of_pop_file)
    if len(input_file.get()) < 1:
        input_file.insert(0, "Select file...")
def save_settings(settings_list,settingssheet = "settingssheet.txt"):
    settingssheet = open(settingssheet, 'w+')

    for listitem in settings_list:
        settingssheet.write('%s\n' % listitem)
    settingssheet.close()
def if_settingssheet_exists_fill_values(input_file, max_distance, grid_resolution, width, tape_type, width_for_edge_detection, calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution):
    if os.path.isfile('./settingssheet.txt'):

        t = "True" + '\n'
        try:
            settingssheet = open('./settingssheet.txt')
            if not settingssheet.readline() == "Select file....":
                settingssheet.seek(0)

                input_file.delete(0,"end")
                input = settingssheet.readline()
                input_file.insert(0, input[0:-1])

                width.delete(0, 'end')
                width.insert(0, float(settingssheet.readline()))

                tape_type.delete(0, 'end')
                tape_type.insert(0, str(settingssheet.readline())[:-1])

                grid_resolution.delete(0, 'end')
                grid_resolution.insert(0, int(settingssheet.readline()))

                width_for_edge_detection.delete(0, 'end')
                width_for_edge_detection.insert(0, float(settingssheet.readline()))

                max_distance.delete(0, 'end')
                max_distance.insert(0, float(settingssheet.readline()))

                if settingssheet.readline() == t:
                    calc_2D_Solution.set(True)
                else:
                    calc_2D_Solution.set(False)
                if settingssheet.readline() == t:
                    calc_2D_with_edge_detection.set(True)
                else:
                    calc_2D_with_edge_detection.set(False)

                if settingssheet.readline() == t:
                    calc_3D_Solution.set(True)
                else:
                    calc_3D_Solution.set(False)


                settingssheet.close()
        except:
            print("Please delete settingssheet.txt")
            settingssheet.close()
def if_settingssheet_EA_exists_fill_values(individual_optimization,turn_around_normal,adap_mutation,
                                        gamma_d, gamma_d2, gamma_d3, gamma_d4, gamma_l, gamma_l2,
                                        gamma_l3, gamma_l4, gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, gamma_ps,
                                        gamma_ps2, gamma_ps3, gamma_ps4,
                                        num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_crossover,
                                        p_mutate_range, p_mutation, pop_size,
                                        pointspersection, equidistant_pts_between_bendpts, step_size):

    if os.path.isfile('./settingssheet_EA.txt'):

        t = "True" + '\n'
        #try:
        settingssheet = open('./settingssheet_EA.txt')

        # First 4 lines we don´t have to read in. use2D/2DE/3D dependent on preprocessor GUI settings

        settingssheet.readline() #use_2D
        settingssheet.readline() #use_2DE
        settingssheet.readline() #use_3D

        if settingssheet.readline() == t:
            individual_optimization.set(True)
        else:
            individual_optimization.set(False)

        if settingssheet.readline() == t:
            turn_around_normal.set(True)
        else:
            turn_around_normal.set(False)


        if settingssheet.readline() == t:
            adap_mutation.set(True)
        else:
            adap_mutation.set(False)


        gamma_d.delete(0, 'end')
        gamma_d.insert(0, float(settingssheet.readline()))
        gamma_d2.delete(0, 'end')
        gamma_d2.insert(0, float(settingssheet.readline()))
        gamma_d3.delete(0, 'end')
        gamma_d3.insert(0, float(settingssheet.readline()))
        gamma_d4.delete(0, 'end')
        gamma_d4.insert(0, float(settingssheet.readline()))

        gamma_l.delete(0, 'end')
        gamma_l.insert(0, float(settingssheet.readline()))
        gamma_l2.delete(0, 'end')
        gamma_l2.insert(0, float(settingssheet.readline()))
        gamma_l3.delete(0, 'end')
        gamma_l3.insert(0, float(settingssheet.readline()))
        gamma_l4.delete(0, 'end')
        gamma_l4.insert(0, float(settingssheet.readline()))

        gamma_ps.delete(0, 'end')
        gamma_ps.insert(0, float(settingssheet.readline()))
        gamma_ps2.delete(0, 'end')
        gamma_ps2.insert(0, float(settingssheet.readline()))
        gamma_ps3.delete(0, 'end')
        gamma_ps3.insert(0, float(settingssheet.readline()))
        gamma_ps4.delete(0, 'end')
        gamma_ps4.insert(0, float(settingssheet.readline()))

        gamma_pe.delete(0, 'end')
        gamma_pe.insert(0, float(settingssheet.readline()))
        gamma_pe2.delete(0, 'end')
        gamma_pe2.insert(0, float(settingssheet.readline()))
        gamma_pe3.delete(0, 'end')
        gamma_pe3.insert(0, float(settingssheet.readline()))
        gamma_pe4.delete(0, 'end')
        gamma_pe4.insert(0, float(settingssheet.readline()))

        num_gen.delete(0, 'end')
        num_gen.insert(0, int(settingssheet.readline()))

        num_gen_set2.delete(0, 'end')
        num_gen_set2.insert(0, float(settingssheet.readline()))
        num_gen_set3.delete(0, 'end')
        num_gen_set3.insert(0, float(settingssheet.readline()))
        num_gen_set4.delete(0, 'end')
        num_gen_set4.insert(0, float(settingssheet.readline()))




        p_crossover.delete(0, 'end')
        p_crossover.insert(0, float(settingssheet.readline()))

        p_mutate_range.delete(0, 'end')
        p_mutate_range.insert(0, float(settingssheet.readline()))

        p_mutation.delete(0, 'end')
        p_mutation.insert(0, float(settingssheet.readline()))

        pop_size.delete(0, 'end')
        pop_size.insert(0, int(settingssheet.readline()))

        pointspersection.delete(0, 'end')
        pointspersection.insert(0, int(settingssheet.readline()))

        if settingssheet.readline() == t:
            equidistant_pts_between_bendpts.set(True)
        else:
            equidistant_pts_between_bendpts.set(False)

        step_size.delete(0, 'end')
        step_size.insert(0, float(settingssheet.readline()))
        settingssheet.close()
        #except:
          #  print("Please delete settingssheet_EA.txt")
          #  settingssheet.close()

# Load settings
def GUI_select_folder_directory():
    global folder_directory
    master = Tk()
    master.protocol("WM_DELETE_WINDOW", sys.exit)
    Label(master, justify=LEFT, text=" ").grid(row=11, sticky=W)
    folder_directory = Entry(master)

    folder_directory.insert(0, "Select folder with settings and population...")
    folder_directory.grid(row=20, columnspan=3, sticky=W + E + N + S)
    Button(text='Browse..', command=lambda: select_folder(folder_directory)).grid(row=20, column=3, sticky=W)
    Button(master, text='Abort', command=sys.exit).grid(row=1000, column=0, pady=4)
    Label(master, justify=LEFT, text="                                                                      ").grid(row=1000, column=1,sticky=W)
    Button(master, text='Start', command=master.quit).grid(row=1000, column=3, pady=4)
    mainloop()  # Executes GUI

    folder_directory = str(folder_directory.get()+"/")
    master.destroy()
    return folder_directory
def select_folder(folder_directory):

    name_of_folder = filedialog.askdirectory(initialdir="./Testruns/", title="Select file")

    folder_directory.delete(0, 'end')
    folder_directory.insert(0, name_of_folder)
    if len(folder_directory.get()) < 1:
        folder_directory.insert(0, "Select stl-file...")
def copy_loaded_settings_to_main_folder(sub_dir):
    shutil.copy(sub_dir + 'settingssheet.txt', "./")
    shutil.copy(sub_dir + 'settingssheet_EA.txt', "./")
    shutil.copy(sub_dir + 'startparameter.txt', "./")
    try:
        shutil.copy(sub_dir + 'start_chromo_2D.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'start_chromo_3D.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'start_chromo_2DE.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_main_after_EA.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_2D_after_EA.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_3D_after_EA.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_2DE_after_EA.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_main_Start.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_2D_Start.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_3D_Start.txt', "./")
    except:
        pass
    try:
        shutil.copy(sub_dir + 'population_2DE_Start.txt', "./")
    except:
        pass
def delete_old_population_and_start_chromo():
    # Delete 2D
    if os.path.exists("start_chromo_2D.txt"):
        os.remove("start_chromo_2D.txt")
    if os.path.exists("population_2D_Start.txt"):
        os.remove("population_2D_Start.txt")
    if os.path.exists("population_2D_after_EA.txt"):
        os.remove("population_2D_after_EA.txt")

    # Delete 2DE
    if os.path.exists("start_chromo_2DE.txt"):
        os.remove("start_chromo_2DE.txt")
    if os.path.exists("population_2DE_Start.txt"):
        os.remove("population_2DE_Start.txt")
    if os.path.exists("population_2DE_after_EA.txt"):
        os.remove("population_2DE_after_EA.txt")

    # Delete 3D
    if os.path.exists("start_chromo_3D.txt"):
        os.remove("start_chromo_3D.txt")
    if os.path.exists("population_3D_Start.txt"):
        os.remove("population_3D_Start.txt")
    if os.path.exists("population_3D_after_EA.txt"):
        os.remove("population_3D_after_EA.txt")

    # Delete main
    if os.path.exists("fitness_over_generation.txt"):
        os.remove("fitness_over_generation.txt")
    if os.path.exists("population_main_Start.txt"):
        os.remove("population_main_Start.txt")
    if os.path.exists("population_main_after_EA.txt"):
        os.remove("population_main_after_EA.txt")
# Save startparameter
def save_startparameter(L_aim, l_factor, Start_direction_prep_fromstart, Start_normal_prep_fromstart, patch_end,patch_start,amount_of_bends):
    file1 = open("startparameter.txt", "w")
    file1.write('%s\n' % L_aim)
    file1.write('%s\n' % l_factor)
    file1.write('%s\n' % patch_start[0])
    file1.write('%s\n' % patch_start[1])
    file1.write('%s\n' % patch_start[2])
    file1.write('%s\n' % patch_end[0])
    file1.write('%s\n' % patch_end[1])
    file1.write('%s\n' % patch_end[2])
    file1.write('%s\n' % Start_normal_prep_fromstart[0])
    file1.write('%s\n' % Start_normal_prep_fromstart[1])
    file1.write('%s\n' % Start_normal_prep_fromstart[2])
    file1.write('%s\n' % Start_direction_prep_fromstart[0])
    file1.write('%s\n' % Start_direction_prep_fromstart[1])
    file1.write('%s\n' % Start_direction_prep_fromstart[2])
    file1.write('%s\n' % amount_of_bends)
    file1.close()

# Load parameters
def load_settings(sub_dir):
    load_settingssheet(sub_dir)
    load_startparameter(sub_dir)
    load_EA_settings(sub_dir)
# load settingsheets
def load_settingssheet(sub_dir):
    global z_grid_values_trendline_KOS, max_x, max_y, min_x, min_y, fix_number_of_pts, pointspersection, equidistant_pts_between_bendpts, step_size,\
        input_file, testpatch, width, tape_type, grid_resolution, width_for_edge_detection, max_distance,calc_2D_Solution,calc_3D_Solution,calc_2D_with_edge_detection

    settingssheet = open(sub_dir+'settingssheet.txt')
    input = settingssheet.readline()
    input_file = input[0:-1]
    testpatch = trimesh.load(input_file)
    stlprep3_6.calc_trendline_of_geometry_from_stl_file(input_file)
    width = float(settingssheet.readline())
    tape_type = str(settingssheet.readline())
    grid_resolution = int(settingssheet.readline())
    width_for_edge_detection = float(settingssheet.readline())
    max_distance = float(settingssheet.readline())
    if settingssheet.readline() == "True" + '\n':
        calc_2D_Solution = True
    else:
        calc_2D_Solution = False
    if settingssheet.readline() == "True" + '\n':
        calc_2D_with_edge_detection = True
    else:
        calc_2D_with_edge_detection = False
    if settingssheet.readline() == "True" + '\n':
        calc_3D_Solution = True
    else:
        calc_3D_Solution = False


    settingssheet.close()

    # For calculating the signd distance, we need the interpolation. If loaded the parameters,
    [grid_x, max_x, max_y, min_x, min_y, z_grid_values_trendline_KOS] = stlprep3_6.calc_trendline_and_interpolate(grid_resolution,input_file)
def load_startparameter(sub_dir):
    global L_aim, l_factor, patch_start, patch_end, Start_normal_prep_fromstart, Start_direction_prep_fromstart,amount_of_bends
    file1 = open(sub_dir+"startparameter.txt", "r")
    L_aim = float(file1.readline())
    l_factor = float(file1.readline())
    patch_start_x = float(file1.readline())
    patch_start_y = float(file1.readline())
    patch_start_z = float(file1.readline())
    patch_end_x = float(file1.readline())
    patch_end_y = float(file1.readline())
    patch_end_z = float(file1.readline())
    Start_normal_prep_fromstart_x = float(file1.readline())
    Start_normal_prep_fromstart_y = float(file1.readline())
    Start_normal_prep_fromstart_z = float(file1.readline())
    Start_direction_prep_fromstart_x = float(file1.readline())
    Start_direction_prep_fromstart_y = float(file1.readline())
    Start_direction_prep_fromstart_z = float(file1.readline())
    amount_of_bends = int(file1.readline())
    file1.close()
    patch_start = np.asarray([patch_start_x, patch_start_y, patch_start_z])
    patch_end = np.asarray([patch_end_x, patch_end_y, patch_end_z])
    Start_normal_prep_fromstart = np.asarray(
        [Start_normal_prep_fromstart_x, Start_normal_prep_fromstart_y, Start_normal_prep_fromstart_z])
    Start_direction_prep_fromstart = np.asarray(
        [Start_direction_prep_fromstart_x, Start_direction_prep_fromstart_y, Start_direction_prep_fromstart_z])
def load_EA_settings(sub_dir):
    settingssheet = open(sub_dir+"settingssheet_EA.txt", "r")
    global use_2D_with_edge_detection, use_2D_Solution, use_3D_Solution, individual_optimization, turn_around_normal,gamma_d, gamma_d2, gamma_d3, gamma_d4, gamma_l, gamma_l2, gamma_l3, gamma_l4, gamma_ps, gamma_ps2, gamma_ps3, gamma_ps4, \
        gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, num_gen_set2, num_gen_set3, num_gen_set4, pop_size, num_gen, \
        chromo_resolution, p_mutation, adap_mutation, p_mutate_range, p_crossover,fix_number_of_pts, pointspersection,equidistant_pts_between_bendpts,step_size


    if settingssheet.readline() == "True" + '\n':
        use_2D_with_edge_detection = True
    else:
        use_2D_with_edge_detection = False
    if settingssheet.readline() == "True" + '\n':
        use_2D_Solution = True
    else:
        use_2D_Solution = False
    if settingssheet.readline() == "True" + '\n':
        use_3D_Solution = True
    else:
        use_3D_Solution = False

    if settingssheet.readline() == "True" + '\n':
        individual_optimization = True
    else:
        individual_optimization = False

    if settingssheet.readline() == "True" + '\n':
        turn_around_normal = True
    else:
        turn_around_normal = False

    if settingssheet.readline() == "True" + '\n':
        adap_mutation = True
    else:
        adap_mutation = False

    gamma_d = float(settingssheet.readline())
    gamma_d2 = float(settingssheet.readline())
    gamma_d3 = float(settingssheet.readline())
    gamma_d4 = float(settingssheet.readline())
    gamma_l = float(settingssheet.readline())
    gamma_l2 = float(settingssheet.readline())
    gamma_l3 = float(settingssheet.readline())
    gamma_l4 = float(settingssheet.readline())
    gamma_ps = float(settingssheet.readline())
    gamma_ps2 = float(settingssheet.readline())
    gamma_ps3 = float(settingssheet.readline())
    gamma_ps4 = float(settingssheet.readline())
    gamma_pe = float(settingssheet.readline())
    gamma_pe2 = float(settingssheet.readline())
    gamma_pe3 = float(settingssheet.readline())
    gamma_pe4 = float(settingssheet.readline())
    num_gen = int(settingssheet.readline())
    num_gen_set2 = float(settingssheet.readline())
    num_gen_set3 = float(settingssheet.readline())
    num_gen_set4 = float(settingssheet.readline())




    p_crossover = float(settingssheet.readline())
    p_mutate_range = float(settingssheet.readline())
    p_mutation = float(settingssheet.readline())
    pop_size = int(settingssheet.readline())

    pointspersection = int(settingssheet.readline())

    if settingssheet.readline() == "True" + '\n':
        equidistant_pts_between_bendpts = True
    else:
        equidistant_pts_between_bendpts = False
    step_size = float(settingssheet.readline())


    settingssheet.close()


# Save and load Populations:
def save_current_population(p, string_name_ending):
    file2 = open("population" + string_name_ending + ".txt", "w")
    for j in range(p.numChromosomes):
        file2.write('%s\n' % p.currentGeneration[j].genes)
    file2.close()
def load_population(pop_size, file):
    generation = []
    file2 = open(file, "r")
    for j in range(pop_size):
        chromo_read = file2.readline()[1:-2]
        chromo_read = chromo_read.split(",")
        try: chromo_read = list(map(int, chromo_read))
        except:
            print("You have to use same population size then the one you load.")
            exit()

        generation.append(chromo_read)
    file2.close()
    return generation


# Save and load start chromo:
def save_start_chromo(chromo, string_name_ending):
    file2 = open("start_chromo" + string_name_ending + ".txt", "w")
    file2.write(chromo)
    file2.close()
def load_start_chromo(file):
    start_chromo = []
    file2 = open(file, "r")
    start_chromo = file2.readline()[1:-1]
    start_chromo = start_chromo.split(",")
    start_chromo = list(map(int, start_chromo))
    file2.close()
    return start_chromo

""" 
Not used at the moment since we saved start_chromosoms.
Eventually usefull if we want to load the INDIVIDUAL Otimized start chromosoms. 

# In population_main_Start the individual optimized start chromosoms are used, when they are calculated
def load_last_start_chromosoms_used_for_optimization_all(pop_size, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution, subdir):
    global startchromo3D, startchromo2D_edge, startchromo2D
    individuals_all = load_population(pop_size, subdir + "population_main_Start.txt")
    if use_3D_Solution and use_2D_with_edge_detection and use_2D_Solution:
        startchromo3D = individuals_all[0]  # startchromo ist first in start population
        startchromo2D_edge = individuals_all[1]
        startchromo2D = individuals_all[2]
    elif use_3D_Solution and use_2D_with_edge_detection:
        startchromo3D = individuals_all[0]  # startchromo ist first in start population
        startchromo2D_edge = individuals_all[1]
    elif use_2D_with_edge_detection and use_2D_Solution:
        startchromo2D_edge = individuals_all[0]  # startchromo ist first in start population
        startchromo2D = individuals_all[1]
    elif use_3D_Solution and use_2D_Solution:
        startchromo3D = individuals_all[0]  # startchromo ist first in start population
        startchromo2D = individuals_all[1]
    elif use_3D_Solution:
        startchromo3D = individuals_all[0]
    elif use_2D_with_edge_detection:
        startchromo2D_edge = individuals_all[0]
    elif use_2D_Solution:
        startchromo2D = individuals_all[0]
    use_3D_Solution_last = use_3D_Solution
    use_2D_Solution_last = use_2D_Solution
    use_2D_with_edge_detection_last = use_2D_with_edge_detection
    return use_2D_Solution_last, use_2D_with_edge_detection_last, use_3D_Solution_last # # ##
"""

# Load preprocessor start chromosoms
def load_preprocessor_start_chromosoms(calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution, subdir):
    global startchromo3D, startchromo2D_edge, startchromo2D
    startchromo3D, startchromo2D_edge, startchromo2D = [],[],[]


    if calc_3D_Solution:
        startchromo3D = load_start_chromo(subdir+"start_chromo_3D.txt")
    if calc_2D_with_edge_detection:
        startchromo2D_edge = load_start_chromo(subdir + "start_chromo_2DE.txt")
    if calc_2D_Solution:
        startchromo2D = load_start_chromo(subdir + "start_chromo_2D.txt")

# Kinematic description of the patch   COMMENT_DB: This is the translation of values suitable for the evolutionary algorithm!
def ListOfPoints(chromo):  # Comment_DB: chromo not defined elsewhere. chromo here is a parameter. Function definition.
    # Alpha, beta and length from chromosome and translated
    [alpha_list,
     beta_list,
     length_list] = translate_alpha_beta_length_from_chromo(chromo)

    # Start point and direction
    [Start_direction,
     Start_normale_gamma,
     Start_point] = calc_start_point_direction_normal_vector(chromo)

    # Direction vector at each beding point
    direction_vector_list, normal_vector_list = calc_direction_vectors(Start_direction, Start_normale_gamma, alpha_list,
                                                                       beta_list,
                                                                       length_list)

    # Lengths of left and right side from tape,
    [delta_length_start_bend,
     length_left_list,
     length_right_list,
     delta_length_at_bendpoint] = calc_delta_length_start_and_side_lengths(alpha_list, length_list)

    # Eckpunkte
    Start_point_left = Start_point - np.cross(Start_direction,
                                              Start_normale_gamma) * width / 2 + delta_length_start_bend * Start_direction
    Start_point_right = Start_point + np.cross(Start_direction,
                                               Start_normale_gamma) * width / 2 - delta_length_start_bend * Start_direction

    ## mid line after startpoint
    mid_point_list = calc_points_from_start_directions_lengths(Start_point, direction_vector_list, length_list)
    left_point_list = calc_points_from_start_directions_lengths(Start_point_left, direction_vector_list,
                                                                length_left_list)
    right_point_list = calc_points_from_start_directions_lengths(Start_point_right, direction_vector_list,
                                                                 length_right_list)

    # Fill in points between the bend points on the center line. Either fixed number of points or equidistant
    mid_point_filled_up_list = calc_filled_up_points(direction_vector_list, length_list, mid_point_list)
    left_points_filled_up_list = calc_filled_up_points(direction_vector_list, length_left_list, left_point_list)
    right_points_filled_up_list = calc_filled_up_points(direction_vector_list, length_right_list, right_point_list)

    tape_section_index = []
    tape_sections = []
    for j in range(1, len(mid_point_list)):
        try:
            index = [i for i, x in enumerate(mid_point_filled_up_list[:, 0]) if x == mid_point_list[j][0]][0] + 1
        except:
            index = -1

        if tape_section_index == []:
            tape_sections.append(mid_point_filled_up_list[0:index])
        else:
            tape_sections.append(mid_point_filled_up_list[tape_section_index[-1]:index])

        tape_section_index.append(index)

    all_patch_points_filled_up = np.concatenate(
        (mid_point_filled_up_list, left_points_filled_up_list, right_points_filled_up_list), axis=0)

    start = mid_point_filled_up_list[0]  # Comment_DB: Change to beginning of list
    end = mid_point_filled_up_list[-1]

    # Only corner and bend points of the tape for easy visualization
    patch_visualisation_points = []

    for i in range(len(left_point_list)):
        patch_visualisation_points.append(left_point_list[i])
        patch_visualisation_points.append(right_point_list[i])

    patch_visualisation_points = np.stack(patch_visualisation_points, axis=0)

    return all_patch_points_filled_up, start, end, patch_visualisation_points, length_list, alpha_list, beta_list, Start_point, Start_direction, direction_vector_list, normal_vector_list, mid_point_list, tape_sections, delta_length_at_bendpoint
# Functions in ListofPoints
def calc_delta_length_start_and_side_lengths(alpha_list, length_list):
    delta_length_start_bend = 0
    delta_length_at_bendpoint = [delta_length_start_bend]
    length_left_list = []
    length_right_list = []

    i = 1
    while i <= len(length_list):
        if alpha_list[i - 1] > math.pi / 2:
            # Delta at bend i
            delta_length_at_bendpoint_i = +(width / 2) * math.tan(math.pi - alpha_list[i - 1])
        else:  # alpha_list[i] < math.pi / 2:
            # Delta at bend i
            delta_length_at_bendpoint_i = -(width / 2) * math.tan(alpha_list[i - 1])

        # Length
        length_left_new = length_list[i - 1] - delta_length_at_bendpoint_i + delta_length_at_bendpoint[i - 1]
        length_right_new = length_list[i - 1] + delta_length_at_bendpoint_i - delta_length_at_bendpoint[i - 1]

        length_left_list.append(length_left_new)
        length_right_list.append(length_right_new)

        delta_length_at_bendpoint.append(delta_length_at_bendpoint_i)

        i += 1
    return delta_length_start_bend, length_left_list, length_right_list, delta_length_at_bendpoint
def calc_filled_up_points(direction_vector_list, length_list, point_list):
    filled_up_list = []
    if equidistant_pts_between_bendpts:

        for i, length in enumerate(length_list):
            if math.floor(length / step_size) < 0:
                num_of_points_to_insert = 0
            else:
                num_of_points_to_insert = math.floor(length / step_size)
            a_new = point_list[i][np.newaxis, :] + np.outer(
                np.linspace(0, length, num_of_points_to_insert, endpoint=True), direction_vector_list[i])

            filled_up_list.append(a_new)
    else:

        for i, length in enumerate(length_list):
            a_new = point_list[i][np.newaxis, :] + np.outer(np.linspace(0, length, pointspersection, endpoint=True),
                                                            direction_vector_list[i])
            filled_up_list.append(a_new)
    filled_up_list = np.concatenate(filled_up_list)
    return filled_up_list
def calc_points_from_start_directions_lengths(Start_point, direction_vector_list, length_list):
    midpoint_list = [Start_point]
    old_p = Start_point
    for i, length in enumerate(length_list):
        direction_vector_length = length * direction_vector_list[i, :]
        new_p = old_p + direction_vector_length
        midpoint_list.append(new_p)
        old_p = new_p
    return midpoint_list
def calc_direction_vectors(Start_direction, Start_normale_gamma, alpha_list, beta_list, length_list):
    direction_vector_list = [Start_direction]
    normal_vector_list = [Start_normale_gamma]
    #### Vector calculation for tape page after the starting point
    for i in range(1, len(length_list)):

        # Rotate Direction Vector around alpha
        if alpha_list[i - 1] < math.pi / 2:
            direction_rotation_alpha = Quaternion(axis=normal_vector_list[i - 1],
                                                  angle=(-alpha_list[i - 1] - (math.pi) / 2)).rotate(
                direction_vector_list[i - 1])

        else:
            direction_rotation_alpha = Quaternion(axis=normal_vector_list[i - 1],
                                                  angle=(-alpha_list[i - 1] + (math.pi) / 2)).rotate(
                direction_vector_list[i - 1])

        # Rotate new Direction Vector around beta
        direction_rotation_alpha_beta = Quaternion(axis=direction_rotation_alpha, angle=-beta_list[i - 1]).rotate(
            direction_vector_list[i - 1])
        # Save Direction vector
        direction_vector_list.append(direction_rotation_alpha_beta)
        n_new = Quaternion(axis=direction_rotation_alpha, angle=-beta_list[i - 1]).rotate(normal_vector_list[i - 1])
        normal_vector_list.append(n_new)

    normal_vector_list = np.stack(normal_vector_list)
    direction_vector_list = np.stack(direction_vector_list)
    return direction_vector_list, normal_vector_list
def calc_start_point_direction_normal_vector(chromo):
    # Starting point variation translated from chromosome
    [[var_start_pt_x, var_start_pt_y, var_start_pt_z],
     variation_start_alpha, variation_start_beta, var_start_n_gamma] = translate_start_varriation_from_chomo(chromo)

    Start_point = np.concatenate(np.array(
        [[patch_start[0] + L_aim * 0.5 * var_start_pt_x], [patch_start[1] + L_aim * 0.5 * var_start_pt_y],
         [patch_start[2] + L_aim * 0.5 * var_start_pt_z]]))

    ##################################################
    # Rotate Direction Vector around alpha
    if variation_start_alpha < math.pi / 2:
        direction_rotation_alpha = Quaternion(axis=Start_normal_prep_fromstart,
                                              angle=(-variation_start_alpha - (math.pi) / 2)).rotate(
            Start_direction_prep_fromstart)

    else:
        direction_rotation_alpha = Quaternion(axis=Start_normal_prep_fromstart,
                                              angle=(-variation_start_alpha + (math.pi) / 2)).rotate(
            Start_direction_prep_fromstart)

    # Rotate new Start_direction_prep_fromstart around direction_rotation_alpha
    Start_direction = Quaternion(axis=direction_rotation_alpha, angle=(variation_start_beta)).rotate(
        Start_direction_prep_fromstart)
    Start_normal = Quaternion(axis=direction_rotation_alpha, angle=(variation_start_beta)).rotate(
        Start_normal_prep_fromstart)
    Start_normale_gamma = Quaternion(axis=Start_direction, angle=(var_start_n_gamma)).rotate(
        Start_normal)  # Comment_DB: start_n_strich rotated about start_r

    return Start_direction, Start_normale_gamma, Start_point
def translate_start_varriation_from_chomo(chromo):
    # From Gen Value(0-c_max) to Startvariation (0 +/- var_range)
    var_range = 0.8
    # Last vales in chromo ar start variation variables
    variation_start = [(0 - var_range + (var_range / (chromo_resolution / 2)) * gen_value) for gen_value in
                       chromo[-6:-3:1]]

    if chromo[-3] < chromo_resolution / 2:
        variation_start_alpha = (135 + (chromo[-3] * 45 / (chromo_resolution / 2))) * 2 * math.pi / 360
    else:
        variation_start_alpha = ((chromo[-3] - chromo_resolution / 2) * 45 / (
                chromo_resolution / 2)) * 2 * math.pi / 360  # Quadratische Vert. von 135°-180°

    variation_start_beta = (chromo[-2] * (180 / chromo_resolution) - 90) * 2 * math.pi / 360

    gamma_max = 10  # [degree] Maximum tilt angle for Start_n
    gamma_max_rad = gamma_max * (2 * math.pi / 360)

    var_start_n_gamma = -gamma_max_rad + gamma_max_rad / (chromo_resolution / 2) * chromo[-1]

    return variation_start, variation_start_alpha, variation_start_beta, var_start_n_gamma
def translate_alpha_beta_length_from_chromo(chromo):
    l_list = []
    alpha_list = []
    beta_list = []
    for i in range(0, len(chromo) - 5, 3):
        l_list.append(chromo[i] * l_factor)
    for i in range(1, len(chromo) - 4, 3):

        if chromo[i] < chromo_resolution / 2:
            alpha = (135 + (chromo[i] * 45 / (chromo_resolution / 2))) * 2 * math.pi / 360
        else:
            alpha = ((chromo[i] - chromo_resolution / 2) * 45 / (
                    chromo_resolution / 2)) * 2 * math.pi / 360  # Quadratische Vert. von 135°-180°
        alpha_list.append(alpha)
    for i in range(2, len(chromo) - 3, 3):
        beta = (chromo[i] * (180 / chromo_resolution) - 90) * 2 * math.pi / 360
        beta_list.append(beta)
    return alpha_list, beta_list, l_list  # beta in radians, length in mm


# Calculation of the fitness of a chromosome
def Fitness(chromo, gen_Num, test_fitness=False):  # Comment DKu_Wenzel L_aim=L_aim
    LoP = ListOfPoints(chromo)

    # Distance_fitness
    distance_fit, avg_dist, area_negativ_distances = calc_distance_fitness(L_aim, testpatch, LoP)

    # Lenght_fittnes
    length_fit = calc_length_fitness(L_aim, chromo, l_factor)

    # Border_fittnes start and end
    border_fit_end, border_fit_start = calc_border_fitness(patch_start, patch_end, LoP)

    # Adaptiv gamma
    if test_fitness:
        gamma_d_hat, gamma_l_hat, gamma_pe_hat, gamma_ps_hat = 1, 1, 1, 1
    else:
        gamma_d_hat, gamma_l_hat, gamma_pe_hat, gamma_ps_hat = evalute_adaptiv_gamma_(gen_Num)

    fitness = distance_fit * gamma_d_hat + length_fit * gamma_l_hat + border_fit_end * gamma_pe_hat + border_fit_start * gamma_ps_hat

    # Penalty functions
    penalty_l_min, l_min_list = calc_length_penalty(LoP)

    penalty_alpha_max = calc_alpha_penalty(LoP)

    penalty_beta_max = calc_beta_penalty(LoP)

    penalty_negativ_position = calc_negativ_position_penalty(area_negativ_distances)

    penalty_beta_minimize = calc_beta_minimize_penalty(LoP, l_min_list)
    # Calculate Fitness
    # todo: penalty factor just for positiv fitness. Is negativ fitness wanted?

    if fitness > 0:
        if with_production_restriction:
            fitness = penalty_negativ_position * penalty_l_min * penalty_alpha_max * penalty_beta_max * fitness
        elif with_minimize_beta:
            # fitness = fitness - (1 - penalty_beta_minimize) * 500
            fitness = fitness * penalty_beta_minimize

    return fitness, distance_fit, length_fit, border_fit_start, border_fit_end, avg_dist, penalty_l_min, penalty_alpha_max, penalty_beta_max, penalty_negativ_position, penalty_beta_minimize
# Functions in Fitness
def patch_length_in_mm(chromo, l_factor_chromo_mm):
    # Calculates the length of a patch. Can read out chomosomes as class chromosomes or as a simple list.

    if inspect.isclass(chromo):
        lengt = 0
        for i in range(0, len(chromo) - 5, 3):
            lengt = lengt + chromo.genes[i]
    else:
        lengt = 0
        for i in range(0, len(chromo) - 5, 3):
            lengt = lengt + chromo[i]
    return lengt * l_factor_chromo_mm
def calc_border_fitness(patch_start, patch_end, LoP):
    ###PARABOLIC###
    # k_p = (100 - 90) / (5 ** 2)  # Comment_DB: k_p = 0.4
    # border_fit_start = 100 - (stlprep3_6.calc_distance_between_two_points(LoP[1], patch_start) ** 2) * k_p
    # border_fit_end = 100 - (stlprep3_6.calc_distance_between_two_points(LoP[2], patch_end) ** 2) * k_p  # Comment_DB: trial and error for k_p
    ###LINEAR###
    k_p_lin = (100 - 90) / 5
    border_fit_start = 100 - abs((stlprep3_6.calc_distance_between_two_points(LoP[1], patch_start)) * k_p_lin)
    border_fit_end = 100 - abs((stlprep3_6.calc_distance_between_two_points(LoP[2], patch_end)) * k_p_lin)
    ###GAUSSIAN###
    # k_p_gauss = -math.log(9 / 10) / (5 ** 2) #Comment_DB: deviation of 5 mm--> 90
    # border_fit_start = 100 * math.exp(-k_p_gauss * (stlprep3_6.calc_distance_between_two_points(LoP[1], patch_start)) ** 2)
    # border_fit_end = 100 * math.exp(-k_p_gauss * (stlprep3_6.calc_distance_between_two_points(LoP[2], patch_end)) ** 2)
    return border_fit_end, border_fit_start
def calc_length_fitness(L_aim, chromo, l_factor_chromo_mm):
    L = patch_length_in_mm(chromo, l_factor_chromo_mm)
    ###PARABOLIC###
    # k_l = (100 - 50) / ((L_aim * 0.2) ** 2)  # Comment_DB: = 1/128 for L_aim = 400. Higher L_aim yields lower k_l
    # length_fit = 100 - ((L - L_aim) ** 2) * k_l
    ###LINEAR###
    k_l_lin = (100 - 50) / (L_aim * 0.2)
    length_fit = 100 - abs((L - L_aim) * k_l_lin)
    ###GAUSSIAN###
    # k_l_gauss = -math.log(5/10)/((0.2*L_aim) ** 2) #Comment_DB: deviation of 0.2*L_aim --> 50
    # length_fit = 100 * math.exp(-k_l_gauss * (L - L_aim) ** 2)
    return length_fit
def calc_avg_dist(testpatch, LoP):
    points = LoP[0]

    distances_testpatch_currentpatch = trimesh.proximity.closest_point(testpatch, points)[
        1]  # closest_point(testpatch, points)[1] #DKu_Wenzel: signed_distance

    distances_testpatch_currentpatch, area_negativ_distances = calc_direction_sign_distance(
        distances_testpatch_currentpatch, points)

    """# todo: distance per tape section
        distances_testpatch_sections = []
        points_sections = LoP[-2]

        for i in range(len(points_sections)):
            #distances_testpatch_currentpatch = trimesh.proximity.signed_distance(testpatch, points_sections[i])  # closest_point(testpatch, points)[1] #DKu_Wenzel: signed_distance
            try: distances_testpatch_sections.append(trimesh.proximity.signed_distance(testpatch,
                                                                             points_sections[i]))
            except: break
    """
    # Comment_DKu_Wenzel trimesh.proximity.closest_point(..)[1] gives back distances
    avg_dist = sum(abs(distances_testpatch_currentpatch)) / len(distances_testpatch_currentpatch)
    max_abs_dist = max(abs(distances_testpatch_currentpatch))

    return avg_dist, max_abs_dist, area_negativ_distances
def calc_direction_sign_distance(distances_testpatch_currentpatch, points):
    points_trendline_KOS = stlprep3_6.translate_and_rotate_points_from_OLD_to_trendline_KOS(points,
                                                                                            stlprep3_6.trendline_global_KOS,
                                                                                            stlprep3_6.center_point_of_cloud_weighted)

    # Step size
    y_0_grid_point_index = np.asarray(np.round(max_y / (max_y - min_y) * grid_resolution), dtype=np.int32)
    x_0_grid_point_index = np.asarray(np.round(max_x / (max_x - min_x) * grid_resolution), dtype=np.int32)
    dy = (max_y - min_y) / grid_resolution
    dx = (max_x - min_x) / grid_resolution

    area_negativ_distances = 0

    for i in range(len(points_trendline_KOS)):
        y_index = np.asarray(
            np.round(np.add(np.divide(points_trendline_KOS[i][1], dy), (grid_resolution - y_0_grid_point_index))),
            dtype=np.int32)
        x_index = np.asarray(
            np.round(np.add(np.divide(points_trendline_KOS[i][0], dx), (grid_resolution - x_0_grid_point_index))),
            dtype=np.int32)

        if turn_around_normal:  # point > z_grid_values
            try:
                with np.errstate(invalid='ignore'):  # ignore invalid errer (occurs if z_grid_value is NaN)
                    true_if_points_below_surface = (
                                points_trendline_KOS[i][2] > z_grid_values_trendline_KOS[x_index, y_index])
            except:
                true_if_points_below_surface = True  # Exception for index out of bounds → Points outside of geometry. Distance set negativ → higher penalty
        else:
            try:
                with np.errstate(invalid='ignore'):  # ignore invalid errer (occurs if z_grid_value is NaN)
                    true_if_points_below_surface = (
                                points_trendline_KOS[i][2] < z_grid_values_trendline_KOS[x_index, y_index])
            except:
                true_if_points_below_surface = True  # Exception for index out of bounds → Points outside of geometry. Distance set negativ → higher penalty

        if true_if_points_below_surface:
            if i > 0:
                dx_i = (points_trendline_KOS[i][0] - points_trendline_KOS[i - 1][0])
                dy_i = (points_trendline_KOS[i][1] - points_trendline_KOS[i - 1][1])
                width_dx_dy_section = math.pow(dx_i * dx_i + dy_i * dy_i, 1 / 2)

                distances_testpatch_currentpatch[i] = - distances_testpatch_currentpatch[i]
                area_negativ_distances += distances_testpatch_currentpatch[i] * width_dx_dy_section

    return distances_testpatch_currentpatch, area_negativ_distances
def calc_distance_fitness(L_aim, testpatch, LoP):
    #  Average distance calculation
    avg_dist, max_dist, area_negativ_distances = calc_avg_dist(testpatch, LoP)

    if max_dist < 16:  # 16 = min_length of a section. For max_dist = 16, max_factor =~ 0.36
        max_factor = math.exp(-max_dist / 16)  # Dku_Wenzel: If max_dist close to 0, max_factor =~ 1.
    else:
        max_factor = 0.33

    ###PARABOLIC###
    k_d = (100 - 90) / (0.005 ** 2)  # DKu_Wenzel: Changed from 0.05 to 0.005, higher sensitivity
    distance_fit = (100 - k_d * (
                (avg_dist) / L_aim) ** 2) * max_factor  # Comment_DB: max distance fitness is 100 (avg_dist = 0)
    ###LINEAR###
    # k_d_lin = (100 - 90) / 0.005
    # distance_fit = (100 - abs(k_d_lin * (avg_dist/L_aim))) * max_factor
    ###GAUSSIAN###
    # k_d_gauss = -math.log(9/10)/(0.005**2) #Comment_DB: deviation of 0.005L_aim --> 90
    # distance_fit = (100 * math.exp(-k_d_gauss*(avg_dist / L_aim) ** 2)) * max_factor
    return distance_fit, avg_dist, area_negativ_distances
def evalute_adaptiv_gamma_(gen_Num):
    # if not 'p' in globals():
    #    pass
    # else:
    if num_gen_set2 <= gen_Num + 1 < num_gen_set3:
        gamma_d_hat = gamma_d2
        gamma_l_hat = gamma_l2
        gamma_ps_hat = gamma_ps2
        gamma_pe_hat = gamma_pe2
    elif num_gen_set3 <= gen_Num + 1 < num_gen_set4:
        gamma_d_hat = gamma_d3
        gamma_l_hat = gamma_l3
        gamma_ps_hat = gamma_ps3
        gamma_pe_hat = gamma_pe3
    elif num_gen_set4 <= gen_Num + 1 <= num_gen:
        gamma_d_hat = gamma_d4
        gamma_l_hat = gamma_l4
        gamma_ps_hat = gamma_ps4
        gamma_pe_hat = gamma_pe4
    else:  # Comment_DB: If inputs don't make sense, revert back to original gammas for all sets
        gamma_d_hat = gamma_d
        gamma_l_hat = gamma_l
        gamma_ps_hat = gamma_ps
        gamma_pe_hat = gamma_pe
    return gamma_d_hat, gamma_l_hat, gamma_pe_hat, gamma_ps_hat


# Penalty
def calc_beta_penalty(LoP):
    beta_list = LoP[6]
    beta_max = max(beta_list)
    beta_min = min(beta_list)
    if abs(beta_max) > abs(beta_min):  # Roboter restriction for beta: +90° to -55°.
        beta_max_value = 90 / 180 * math.pi
        beta_min_value = -55 / 180 * math.pi
    else:
        beta_max_value = 55 / 180 * math.pi
        beta_min_value = -90 / 180 * math.pi
    a = math.log(1 / 2, math.e) * (-1) / (
                math.pi / 4)  # for delta_beta = 1/4 Pi → penalty_factor = 1/2 #todo: Check behaviour of penalty factor
    penalty_beta_max = 1
    for i in range(len(beta_list)):

        if beta_list[i] > beta_max_value:
            penalty_factor = math.exp((-a * (beta_list[i] - beta_max_value)))

            if penalty_factor < penalty_beta_max:
                penalty_beta_max = penalty_factor

        elif beta_list[i] < beta_min_value:
            penalty_factor = math.exp((a * (beta_list[i] - beta_min_value)))

            if penalty_factor < penalty_beta_max:
                penalty_beta_max = penalty_factor

    return penalty_beta_max
def calc_alpha_penalty(LoP):
    alpha_list = LoP[5]
    penalty_alpha_max = 1
    for i in range(len(alpha_list)):

        if alpha_list[i] > math.pi / 2:
            if alpha_list[i] < 3 * math.pi / 4:
                a = math.log(1 / 3, math.e) * (-1) / (math.pi / 2)  # for alpha = 1/2 Pi → penalty_factor = 1/3
                penalty_factor = math.exp((-a * (3 * math.pi / 4 - alpha_list[i])))
                if penalty_factor < penalty_alpha_max:
                    penalty_alpha_max = penalty_factor

        else:  # alpha < math.pi/2
            if alpha_list[i] > math.pi / 4:
                a = math.log(1 / 3, math.e) * (-1) / (math.pi / 2)  # for alpha 1/2 Pi → penalty_factor = 1/3
                penalty_factor = math.exp((-a * (alpha_list[i] - math.pi / 4)))
                if penalty_factor < penalty_alpha_max:
                    penalty_alpha_max = penalty_factor
    return penalty_alpha_max
def calc_length_penalty(LoP):
    TCP_offset = 8  # mm
    arc_length = 3  # mm
    width_tool = 5  # mm
    delta_L = LoP[-1]
    length_list = LoP[4]
    # l_min = TCP_offset + arc_length + width_tool + delta_L_start + delta_L_end
    penalty_l_max = 1
    l_min_list = []
    for i in range(1, len(length_list) - 1):  # length or first an last tape section can be shorter
        l_min_i = TCP_offset + arc_length + width_tool + abs(delta_L[i]) + abs(delta_L[i + 1])
        l_min_list.append(l_min_i)
        if length_list[i] < l_min_i:
            a = 2.2  # for l_segment = 1/2 l_min → penalty_factor = 1/3 , (-2*ln(1/3))
            penalty_factor = math.exp((-a * (l_min_i - length_list[i]) / l_min_i))
            if penalty_factor < penalty_l_max:
                penalty_l_max = penalty_factor
    return penalty_l_max, l_min_list
def calc_negativ_position_penalty(area_negativ_distances):
    penalty_negativ_position = 1

    if area_negativ_distances < 0:  # if negativ
        # if area_negativ_distances < -L_aim * 2: #  (L_aim * 1mm) = Area, if Tape is around 1mm below surface
        #    penalty_negativ_position = 0.1  # Penalty limitation
        # else:
        penalty_negativ_position = math.exp((area_negativ_distances / (2 * L_aim * 1)))  # With a negativ distance
    return penalty_negativ_position
def calc_beta_minimize_penalty(LoP, l_min_list):
    beta_list = LoP[6]
    length_list = LoP[4]
    length_ratios = []
    length_ratio_setting = 0.9

    for i in range(len(length_list) - 2):  # start and end section doesnt have a l_min
        length_ratio_i = length_list[i + 1] / l_min_list[i]
        length_ratios.append(length_ratio_i)

    sum_beta = 0
    nr_penalty_sections = 0
    a = math.log(1 / 4, math.e) / (math.pi * 10 / 180)  # for beta = 10° → penalty_factor = 1/3

    for i in range(len(length_ratios)):
        beta_1, beta_2 = 1, 1


        if length_ratios[i] < length_ratio_setting and length_ratios[
            i] > 0:  # todo: try with different length_ratios. 1, just smaller sections, 5, also bigger sections strait!?
            # max is length_ratio.max?

            beta_1 = abs(beta_list[i])
            if i != len(length_ratios) - 1:
                beta_2 = abs(beta_list[i + 1])

            smaller_beta = beta_1 if beta_1 < beta_2 else beta_2

            if smaller_beta > math.pi / 180 * 5:  # If beta is small enough to get deleted, dont add it to the sum
                sum_beta += smaller_beta
            nr_penalty_sections += 1

    sum_beta_avg = sum_beta / nr_penalty_sections if nr_penalty_sections > 0 else 0
    penalty_beta_minimize = math.exp((a * (sum_beta_avg)))
    return penalty_beta_minimize


# Creation of start chromosom(Comment_DB: independent of chromominmaxvalue limitation from galileo)
def create_start_chromo(start_2D_or_3D):
    # Takes the start parameters from the geometry preprocessing and converts them into a chromosome with the
    # corresponding resolution. The return value is the chromosome of the starting solution.
    start_chromo = []

    # Fill length1, alpha1, beta1, length2...
    for i in range(len(start_lengths[start_2D_or_3D])):
        start_chromo.append(int(start_lengths[start_2D_or_3D][i] / l_factor))
        if i < len(start_betas[
                       start_2D_or_3D]):  # Comment_DB: range of beta_list compared to range of l_list is smaller by 1

            # start_chromo.append(int(chromo_resolution / 2))  # Comment_DB: Alphas -> zero on default
            if start_alphas[start_2D_or_3D][i] > math.pi / 2:
                start_chromo.append(int(round(
                    (start_alphas[start_2D_or_3D][i] - ((3 / 4) * math.pi)) * (4 / (math.pi)) * (
                                chromo_resolution / 2))))
            else:
                start_chromo.append(int(round(
                    (chromo_resolution / 2) + (start_alphas[start_2D_or_3D][i] / (math.pi / 4)) * (
                                chromo_resolution / 2))))

            beta_chromo = (start_betas[start_2D_or_3D][i] / (math.pi) + 1 / 2) * chromo_resolution
            start_chromo.append(int(round(beta_chromo)))

    # Variation parameters are set to chromo_resolution/2
    for i in range(6):
        start_chromo.append(int(chromo_resolution / 2))
    return start_chromo
def show_chromo(chromo, Titel="Result"):
    LoP = ListOfPoints(chromo)
    points_all_filled_up = LoP[0]
    patch_visualisation_points = LoP[3]

    #############PLOTTING###########
    figure = pyplot.figure()  # Comment_DB: create a new figure
    figure.suptitle(Titel, fontsize=12)

    axes = mplot3d.Axes3D(figure)
    patch_visual = mplot3d.art3d.Poly3DCollection(stlprep3_6.triangle_vectors_of_stl, linewidths=1, alpha=0.5)

    # Comment_DB: stl mesh. Added to show point cloud
    axes.scatter(points_all_filled_up[:, 0], points_all_filled_up[:, 1], points_all_filled_up[:, 2], c='y')
    # Plotting the patch The bend edge points are plotted with triangles.

    for i in range(len(patch_visualisation_points) - 2):
        verts = [
            list(zip(
                [patch_visualisation_points[i][0], patch_visualisation_points[i + 1][0],
                 patch_visualisation_points[i + 2][0]], \
                [patch_visualisation_points[i][1], patch_visualisation_points[i + 1][1],
                 patch_visualisation_points[i + 2][1]], \
                [patch_visualisation_points[i][2], patch_visualisation_points[i + 1][2],
                 patch_visualisation_points[i + 2][2]]))]
        axes.add_collection3d(Poly3DCollection(verts), zs='z')

    # Paint the bends red
    axes.scatter(patch_visualisation_points[:, 0], patch_visualisation_points[:, 1], patch_visualisation_points[:, 2],
                 c='r')

    axes.scatter(patch_start[0], patch_start[1], patch_start[2], c='black')
    axes.scatter(patch_end[0], patch_end[1], patch_end[2], c='black')
    face_color = (0.5, 0.5, 1, 1)  # alternative: matplotlib.colors.rgb2hex([0.5, 0.5, 1])
    patch_visual.set_facecolor(face_color)
    patch_visual.set_edgecolor((0.5, 0.5, 1, 0.5))
    axes.add_collection3d(patch_visual)

    # Show the plot to the screen
    axes.autoscale(enable=False, axis='both')  # you will need this line to change the Z-axis
    axes.set_xbound(-100, 100)
    axes.set_ybound(-50, 150)
    axes.set_zbound(-100, 100)
    pyplot.axis('off')
    pyplot.show(figure)
def initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, startchromo2D, startchromo2D_edge, startchromo3D):
    # Create an object of the Population class
    p = Population(pop_size)  # Comment_DB: pop_size is user input in dialog box
    # Start values from preprocessing are given to the population
    p.startchromo3D = startchromo3D
    p.startchromo2D = startchromo2D
    p.startchromo2D_edge = startchromo2D_edge
    # Should the initialization be done with the start values of the pre-process?
    p.preprocessedInit = (
                use_2D_with_edge_detection or use_2D_Solution or use_3D_Solution)  # Comment_DB: init_preprocess also user input, preprocessedInit is from chromosome class in Galileo module
    p.initrange = 3  # Comment_DB: in chromosome class in Galileo module. Interval for variation of preprocessed gene
    p.p_randInit = 8
    p.p_prepInit = 70
    # Set the fitness function
    p.evalFunc = Fitness  # Comment_DB: The FUNCTION Fitness is assigned, not the lowercase equation fitness! Stores the function. p.evalFunc = Fitness() stores the return value
    # Set the minimum & maximum values and length of a chromosome depending on the number of kinks
    p.chromoMinValues = [0] * (3 * amount_of_bends + 7)
    p.chromoMaxValues = [chromo_resolution] * (3 * amount_of_bends + 7)
    # Integer alleles -> useInteger = 1, Float values -> useInteger = 0
    useInteger = 1
    p.useInteger = useInteger
    # p.useInteger = 0
    # Select the selection method -> Roulette wheel, Rankedselect, Eliteranked
    # p.selectFunc = p.select_Roulette
    # p.selectFunc = p.select_Ranked
    p.selectFunc = p.select_EliteRanked  # Comment_DB: function returns elites[k-1]
    # p.selectFunc = p.select_Roulette #Comment_DKu_Wenzel: Convergence problem with elites?
    # How many chromosomes may survive? #Comment_DKu_Wenzel: That is the number of children!!!
    # p.replacementSize = p.numChromosomes
    p.replacementSize = p.numChromosomes * 3  # 3x more children then parents todo: explain:
    # Crossover probability
    p.crossoverRate = p_crossover
    # Selection of the crossover function -> Flat
    p.crossoverFunc = p.crossover_Flat
    # p.crossoverFunc = p.crossover_OnePoint
    # p.crossoverFunc = p.crossover_Uniform
    # Mutations probability
    p.mutationRate = p_mutation
    # Selection of the mutation function
    if (calc_2D_with_edge_detection or calc_2D_Solution or calc_3D_Solution):
        # p.mutateFunc = p.mutate_Default
        p.mutateFunc = p.mutate_Uniform
        p.mutationRange = p_mutate_range
    else:
        p.mutateFunc = p.mutate_Uniform
        # p.mutateFunc = p.mutate_Gauss
        p.mutationRange = p_mutate_range
    # Replacementfunktion: SteadyState, SteadyState without double chromosomes & Generational (only children survive)
    # p.replaceFunc = p.replace_SteadyState
    p.replaceFunc = p.replace_SteadyStateNoDuplicates
    # p.replaceFunc = p.replace_SteadyStateNoDuplicates_just_offsprings
    # p.replaceFunc = p.replace_Generational
    p.maxGenerations = num_gen

    p.testpatch = testpatch
    p.num_gen_set2 = num_gen_set2
    p.num_gen_set3 = num_gen_set3
    p.num_gen_set4 = num_gen_set4

    return p
def save_fitness_over_gen():
    file2 = open("fitness_over_generation.txt", "w")
    list_gen_index = [fitness_list_gen_index, distance_fit_list_gen_index, length_fit_list_gen_index, border_fit_start_list_gen_index, border_fit_end_list_gen_index, mutation_rate_list_gen_index, penalty_alpha_list_gen_index, penalty_beta_list_gen_index, penalty_length_list_gen_index, penalty_negativ_list_gen_index, penalty_minimise_beta_list_gen_index]
    file2.write(str(fitness_list_gen_index[0]) + "\n")
    for list_1 in list_gen_index:
        file2.write(str(list_1[1])+ "\n")

    file2.close()

# Repair and delete bendpoints
def repair_start_chromo(startchromo3D, startchromo2D):
    while len(ListOfPoints(startchromo3D)[4]) != len(ListOfPoints(startchromo2D)[4]):
        # length lists
        listlength3D = ListOfPoints(startchromo3D)[4]
        listlength2D = ListOfPoints(startchromo2D)[4]
        # total length
        totallength3D = sum(listlength3D)
        totallength2D = sum(listlength2D)

        # length normed to total length
        listlength2D_norm = length_normed(listlength2D, totallength2D)
        listlength3D_norm = length_normed(listlength3D, totallength3D)

        if len(listlength3D) < len(listlength2D):
            # 3D is shorter
            startchromo3D = calc_point_of_diff_and_insert_neutral_bendpoint(listlength2D_norm, listlength3D,
                                                                            listlength3D_norm, startchromo3D)
        else:  # 2D shorter
            startchromo2D = calc_point_of_diff_and_insert_neutral_bendpoint(listlength3D_norm, listlength2D,
                                                                            listlength2D_norm, startchromo2D)

    return startchromo3D, startchromo2D

# Functions in repair start chromo
def length_normed(list_lengths, total_length):
    listlength2D_norm = list(np.zeros(len(list_lengths)))
    for i in range(len(list_lengths)):
        listlength2D_norm[i] = list_lengths[i] / total_length
    return listlength2D_norm
def calc_point_of_diff_and_insert_neutral_bendpoint(listlength_longer_norm, listlength_shorter, listlength_shorter_norm,chromo_shorter):
    point_of_differenz = calc_points_of_differenz(listlength_longer_norm, listlength_shorter, listlength_shorter_norm)
    # Two possibilities
    option_1_added_lengths = listlength_longer_norm[point_of_differenz] + listlength_longer_norm[point_of_differenz - 1]
    option_2_added_lengths = listlength_longer_norm[point_of_differenz] + listlength_longer_norm[point_of_differenz + 1]
    if option_1_added_lengths / listlength_shorter_norm[point_of_differenz - 1] > 0.85 and \
            option_1_added_lengths / listlength_shorter_norm[point_of_differenz - 1] < 1.15:

        chromo_shorter_updated = insert_new_bendpoint_and_update_old(listlength_longer_norm, listlength_shorter,
                                                                     option_1_added_lengths,
                                                                     point_of_differenz - 1,
                                                                     chromo_shorter)

    else:  # option_2/listlength3D_norm[point_of_differenz]>0.85:

        chromo_shorter_updated = insert_new_bendpoint_and_update_old(listlength_longer_norm, listlength_shorter,
                                                                     option_2_added_lengths,
                                                                     point_of_differenz,
                                                                     chromo_shorter)
    return chromo_shorter_updated
def calc_points_of_differenz(listlength_longer_norm, listlength_shorter, listlength_shorter_norm):
    compare_2D_3D = list(np.zeros(len(listlength_shorter)))
    point_of_differenz = 0
    for i in range(len(listlength_shorter)):
        compare_2D_3D[i] = listlength_longer_norm[i] / listlength_shorter_norm[i]
        # compare_end_2D_3D[-i-1] = listlength3D_norm[-i-1]/listlength2D_norm[-i-1]
        if compare_2D_3D[i] < 0.80:  # mehr als 25% Abweichung
            point_of_differenz = i
            break
    return point_of_differenz
def insert_new_bendpoint_and_update_old(listlength_longer_norm, listlength_shorter, option_added_lengths,point_of_insertion, chromo_shorter):
    division = listlength_longer_norm[point_of_insertion] / option_added_lengths
    # Insert the new bendingpoint:
    chromo_shorter.insert((point_of_insertion) * 3 + 3,
                          chromo_shorter[(point_of_insertion) * 3 + 2])  # old alpha
    chromo_shorter.insert((point_of_insertion) * 3 + 3,
                          chromo_shorter[(point_of_insertion) * 3 + 1])  # old beta
    chromo_shorter.insert((point_of_insertion) * 3 + 3, int(
        listlength_shorter[point_of_insertion] * (1 - division) / l_factor))  # (1-division) * length
    # Update old/previous bendingpoint with 0° angels:
    chromo_shorter[(point_of_insertion) * 3] = int(
        listlength_shorter[point_of_insertion] * division / l_factor)  # (division) * length
    chromo_shorter[(point_of_insertion) * 3 + 1] = int(chromo_resolution / 2)  # 0° angel
    chromo_shorter[(point_of_insertion) * 3 + 2] = int(chromo_resolution / 2)  # 0° angel

    return chromo_shorter


# Remove
def remove_bendingpoint(chromo, bend_nr):
    LoP = ListOfPoints(chromo)
    array_bendingpoints = LoP[-3]  # bendpoints/midpoints, from start-endpoint
    direction_vectors = LoP[-5]  # directions at each bendpoint, start and endpoint dont have a direction!
    normal_vectors = LoP[-4]  # normal at each bendpoint, start and endpoint dont have a direction!
    length_list = LoP[4]
    new_direction_vector = array_bendingpoints[bend_nr + 1] - array_bendingpoints[bend_nr - 1]
    new_length = np.linalg.norm(new_direction_vector)
    new_direction_vector = stlprep3_6.norm_vector(new_direction_vector)
    new_normal_vector = stlprep3_6.norm_vector(
        (normal_vectors[bend_nr - 1] * (length_list[bend_nr - 1] / new_length) + normal_vectors[bend_nr]) * (
                    length_list[bend_nr] / new_length))
    new_side_vector = np.cross(new_normal_vector, new_direction_vector)
    new_side_vector = stlprep3_6.norm_vector(new_side_vector)

    trendline_patch_new = np.stack([new_direction_vector, new_side_vector, new_normal_vector])

    if bend_nr > 1:
        # Old Vectors, they dont change
        direction_before = direction_vectors[bend_nr - 2]
        normal_before = normal_vectors[bend_nr - 2]
        side_direction_before = np.cross(normal_before, direction_before)
        trendline_patch_before = np.stack([direction_before, side_direction_before, normal_before])

        if bend_nr < len(direction_vectors) - 1:  # When it is not the last bendpoint
            direction_after = direction_vectors[bend_nr + 1]
            normal_after = normal_vectors[bend_nr + 1]
            side_direction_after = np.cross(normal_after, direction_after)
            trendline_patch_after = np.stack([direction_after, side_direction_after, normal_after])

        new_alpha_bend_before, new_beta_bend_before = calc_alpha_beta_before(
            array_bendingpoints, bend_nr, trendline_patch_before, trendline_patch_new)

        # trendline_patch_new = np.stack([new_direction_vector, new_side_direction, new_normal_vector_alpha_beta])

        if bend_nr < len(direction_vectors) - 1:
            new_alpha_bend_after, new_beta_bend_after = calc_alpha_beta_after(array_bendingpoints, bend_nr,
                                                                              trendline_patch_after,
                                                                              trendline_patch_new)

        # Updating Chromosom
        chromo = updating_length_angel_before(bend_nr, chromo, new_alpha_bend_before, new_beta_bend_before, new_length)

        if bend_nr < len(direction_vectors) - 1:
            chromo = updating_angel_after(bend_nr, chromo, new_alpha_bend_after, new_beta_bend_after)

        # Delete the bendpoint in between

        chromo.pop(3 * (bend_nr) - 2)
        chromo.pop(3 * (bend_nr) - 2)
        chromo.pop(3 * (bend_nr) - 2)

    else:  # bend_nr == 1:
        # Old Vectors, they dont change
        direction_before = Start_direction_prep_fromstart
        normal_before = Start_normal_prep_fromstart
        side_direction_before = np.cross(normal_before, direction_before)
        trendline_patch_before = np.stack([direction_before, side_direction_before, normal_before])

        if len(direction_vectors) > 2:
            direction_after = direction_vectors[bend_nr + 1]
            normal_after = normal_vectors[bend_nr + 1]
            side_direction_after = np.cross(direction_after, normal_after)
            trendline_patch_after = np.stack([direction_after, side_direction_after, normal_after])

        new_alpha_bend_before, new_beta_bend_before = calc_alpha_beta_before(
            array_bendingpoints, bend_nr, trendline_patch_before, trendline_patch_new)

        if len(direction_vectors) > 2:
            new_alpha_bend_after, new_beta_bend_after = calc_alpha_beta_after(array_bendingpoints, bend_nr,
                                                                              trendline_patch_after,
                                                                              trendline_patch_new)

        # Updating Chromosom
        # new length before
        chromo[3 * (bend_nr - 1)] = int(round(new_length / l_factor))
        # new alpha before, changing the start_variation parameters
        if new_alpha_bend_before > math.pi / 2:
            chromo[-3] = int(
                round((new_alpha_bend_before - ((3 / 4) * math.pi)) * (4 / (math.pi)) * (chromo_resolution / 2)))
        else:
            chromo[-3] = int(
                round((chromo_resolution / 2) + new_alpha_bend_before / (math.pi / 4) * (chromo_resolution / 2)))
        # new beta before, changing the start_variation parameters
        chromo[-2] = int(round((-new_beta_bend_before / (math.pi) + 1 / 2) * chromo_resolution))

        if len(direction_vectors) > 2:
            chromo = updating_angel_after(bend_nr, chromo, new_alpha_bend_after, new_beta_bend_after)

        # Delete the bendpoint in between

        chromo.pop(3 * (bend_nr) - 2)
        chromo.pop(3 * (bend_nr) - 2)
        chromo.pop(3 * (bend_nr) - 2)

    return chromo
# Functions in remove bendingpoint
def calc_alpha_beta_before(array_bendingpoints, bend_nr, trendline_patch_before, trendline_patch_new):
    new_beta_bend_before = math.acos(np.dot(trendline_patch_new[2], trendline_patch_before[2]))

    try:
        point_on_bending_edge = calc_plane_line_intersection(trendline_patch_before[2], array_bendingpoints[bend_nr - 1],
                                                         trendline_patch_new[1], array_bendingpoints[bend_nr + 1],
                                                         epsilon=1e-6)
        vector_along_bending_edge = stlprep3_6.norm_vector(point_on_bending_edge - array_bendingpoints[bend_nr - 1])

        try:
            new_alpha_bend_before = math.acos(np.dot(trendline_patch_before[1], vector_along_bending_edge))

        except:
            new_alpha_bend_before = 0

    except RuntimeError as error:
        new_alpha_bend_before = 0


    if new_alpha_bend_before > math.pi / 2:
        new_alpha_bend_before = (new_alpha_bend_before - math.pi)

    # new direction vector in patch_KOS
    direction_vector_before_patch_KOS = stlprep3_6.translate_and_rotate_points_from_OLD_to_trendline_KOS(
        np.asarray([trendline_patch_new[0], trendline_patch_before[0]]),
        trendline_patch_before,
        np.asarray([0, 0, 0]))
    # sign beta
    if direction_vector_before_patch_KOS[0, 2] > -0.001:
        new_beta_bend_before = -new_beta_bend_before

    # For chromosom translation
    if new_alpha_bend_before < 0:
        new_alpha_bend_before = math.pi + new_alpha_bend_before
    return new_alpha_bend_before, new_beta_bend_before
def calc_alpha_beta_after(array_bendingpoints, bend_nr, trendline_patch_after, trendline_patch_new):
    new_beta_bend_after = -math.acos(np.dot(trendline_patch_after[2], trendline_patch_new[2]))

    point_on_bending_edge = calc_plane_line_intersection(trendline_patch_new[2], array_bendingpoints[bend_nr + 1],
                                                         trendline_patch_after[1], array_bendingpoints[bend_nr + 2],
                                                         epsilon=1e-6)
    vector_along_bending_edge = stlprep3_6.norm_vector(point_on_bending_edge - array_bendingpoints[bend_nr + 1])

    try:
        new_alpha_bend_after = math.acos(np.dot(trendline_patch_new[1], vector_along_bending_edge))
    except:
        new_alpha_bend_after = 0
    if new_alpha_bend_after > math.pi / 2:
        new_alpha_bend_after = (new_alpha_bend_after - math.pi)

    # new direction vector in patch_KOS
    new_direction_vector_patch_KOS = stlprep3_6.translate_and_rotate_points_from_OLD_to_trendline_KOS(
        np.asarray(
            [trendline_patch_after[0], trendline_patch_after[1], trendline_patch_after[2], trendline_patch_new[0],
             trendline_patch_new[1], trendline_patch_new[2]]),
        trendline_patch_new,
        np.asarray([0, 0, 0]))
    # sign beta
    if new_direction_vector_patch_KOS[0, 2] < -0.001:
        new_beta_bend_after = -new_beta_bend_after
    # sign alpha
    # if new_direction_vector_patch_KOS[0, 1] < -0.001:
    #   new_alpha_bend_after = - new_alpha_bend_after

    if new_alpha_bend_after < 0:
        new_alpha_bend_after = math.pi + new_alpha_bend_after

    return new_alpha_bend_after, new_beta_bend_after
def updating_length_angel_before(bend_nr, chromo, new_alpha_bend_before, new_beta_bend_before, new_length):
    # new length before
    chromo[3 * (bend_nr - 1)] = int(round(new_length / l_factor))
    # new alpha before
    if new_alpha_bend_before > math.pi / 2:
        chromo[3 * (bend_nr - 2) + 1] = int(
            round((new_alpha_bend_before - ((3 / 4) * math.pi)) * (4 / (math.pi)) * (chromo_resolution / 2)))
    else:
        chromo[3 * (bend_nr - 2) + 1] = int(
            round((chromo_resolution / 2) + new_alpha_bend_before / (math.pi / 4) * (chromo_resolution / 2)))
    # new beta before    (chromo_resolution/2)+ (start_alphas[start_2D_or_3D][i]/(math.pi/4)) * (chromo_resolution / 2)
    chromo[3 * (bend_nr - 2) + 2] = int(round((new_beta_bend_before / (math.pi) + 1 / 2) * chromo_resolution))
    return chromo
def updating_angel_after(bend_nr, chromo, new_alpha_bend_after, new_beta_bend_after):
    # new_alpha_bend_after
    if new_alpha_bend_after > math.pi / 2:
        chromo[3 * (bend_nr) + 1] = int(
            round((new_alpha_bend_after - ((3 / 4) * math.pi)) * (4 / (math.pi)) * (chromo_resolution / 2)))
    else:
        chromo[3 * (bend_nr) + 1] = int(
            round((chromo_resolution / 2) + new_alpha_bend_after / (math.pi / 4) * (chromo_resolution / 2)))
    # new beta after
    chromo[3 * (bend_nr) + 2] = int(round((new_beta_bend_after / (math.pi) + 1 / 2) * chromo_resolution))

    return chromo
def calc_plane_line_intersection(plane_normal, plane_point, line_direction, line_point, epsilon=1e-6):
    ndotu = plane_normal.dot(line_direction)
    if abs(ndotu) < epsilon:
        raise RuntimeError("no intersection or line is within plane")

    w = line_point - plane_point
    si = -plane_normal.dot(w) / ndotu
    Psi = w + si * line_direction + plane_point
    return Psi

#Main function
def main():
    global  grid_resolution, max_x, max_y, min_x, min_y, z_grid_values_trendline_KOS, turn_around_normal, use_last_setting, load_preproc_results, step_size, testpatch, tape_type, width, pointspersection, equidistant_pts_between_bendpts, adap_mutation, chromo_resolution, gamma_d, gamma_d2, gamma_d3, gamma_d4, gamma_l, gamma_l2, gamma_l3, gamma_l4, gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, gamma_ps, gamma_ps2, gamma_ps3, gamma_ps4, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_crossover, p_mutate_range, p_mutation, pop_size, useInteger, start_lengths, L_aim, start_betas, start_alphas, patch_start, patch_end, Start_direction_prep_fromstart, Start_normal_prep_fromstart, amount_of_bends, startchromo3D, startchromo2D, startchromo2D_edge, p, time_start, time_end, fitness_list_gen_index, distance_fit_list_gen_index, length_fit_list_gen_index, border_fit_start_list_gen_index, border_fit_end_list_gen_index, mutation_rate_list_gen_index, end , calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution,l_factor, use_2D_with_edge_detection, use_2D_Solution, use_3D_Solution, individual_optimization

    #################### Preprocessor ####################

    # GUI Preprocessor Settings
    [use_last_setting, load_preproc_results, tape_type, width, input_file, max_distance, grid_resolution, width_for_edge_detection, \
    calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution]=get_Preprocessor_Vars_from_GUI()

    # Default values to plot preprocessor results
    global turn_around_normal, equidistant_pts_between_bendpts, step_size, pointspersection, with_production_restriction, with_minimize_beta
    turn_around_normal = False
    with_production_restriction = False
    with_minimize_beta = False
    equidistant_pts_between_bendpts = True
    step_size = 1
    pointspersection = 10

    if use_last_setting:
        try:
            load_settings("./")
        except:
            print("Settings are missing - Check folder")
            exit()

    elif load_preproc_results:
        sub_dir = GUI_select_folder_directory()
        delete_old_population_and_start_chromo()  # Just delete population and start chromosoms, depending on setting they change. Settings get overwritten.
        copy_loaded_settings_to_main_folder(sub_dir) # Copy the loaded settingssheet and startparameter into the current folder to work on it, without changing the initial settings
        try:
            load_settings("./") # Since we copy the settings into our main folder, we can load the settings directly from here.
        except:
            print("Settings are missing - Check folder")
            exit()

    else: # Calculation of start-parameter wih preprocessor :
        delete_old_population_and_start_chromo()
        testpatch = trimesh.load(input_file)

        continue_check = True
        while continue_check:
            [start_lengths,
             L_aim,  # Comment_DB: already in [mm]
             start_betas,
             start_alphas,

             patch_start,
             patch_end,

             Start_direction_prep_fromstart,
             Start_normal_prep_fromstart,

             amount_of_bends,
             z_grid_values_trendline_KOS,
             max_x, max_y, min_x, min_y] = stlprep3_6.startparam(input_file, max_distance, width_for_edge_detection, grid_resolution,  # todo: Start_vector/_normal
                                                      calc_2D_with_edge_detection, calc_3D_Solution)
            l_factor = 0.75 * float(L_aim) / chromo_resolution
            save_startparameter(L_aim,l_factor, Start_direction_prep_fromstart,
                                Start_normal_prep_fromstart,patch_end,patch_start,amount_of_bends)
            # Default value for displaying the results of the preprocessor.
            show_and_save_preprocessor_result(calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution, patch_end,
                                              patch_start)
            continue_check, max_distance, width_for_edge_detection = change_amount_of_bendpoints_GUI(max_distance, width_for_edge_detection)



    # Preproc results in chromosom format
    # Initialize
    startchromo3D = []
    startchromo2D = []
    startchromo2D_edge = []

    if use_last_setting or load_preproc_results:
        try:
            load_preprocessor_start_chromosoms(calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution, "./")
        except:
            print("Start chromosoms are missing - Check folder and settingssheet")
            exit()

    else:  # If not use last settings
        if calc_3D_Solution:
            startchromo3D = create_start_chromo(0)  # 0, for 3D.
        if calc_2D_Solution:
            startchromo2D = create_start_chromo(1)  # 1, for 2D start solution
        if calc_2D_with_edge_detection:
            startchromo2D_edge = create_start_chromo(2)  # 2, for 2D start with edge detection
    #################### Evolutionary Algorithm ####################

    # GUI EA Settings
    individual_optimization = True
    while individual_optimization:
        [use_2D_with_edge_detection, use_2D_Solution, use_3D_Solution, individual_optimization, turn_around_normal, adap_mutation,
     gamma_d, gamma_d2,
     gamma_d3, gamma_d4, gamma_l, gamma_l2, gamma_l3, gamma_l4, gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4, gamma_ps,
     gamma_ps2, gamma_ps3, gamma_ps4, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_crossover, p_mutate_range,
     p_mutation, pop_size,pointspersection, equidistant_pts_between_bendpts, step_size, input_pop1, input_pop2, input_pop3] = get_EA_Vars_from_GUI(calc_3D_Solution, calc_2D_Solution, calc_2D_with_edge_detection)
        # Individual optimization
        if_individual_optimization(individual_optimization, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution)

    # EA initialized with loaded Population
    if input_pop1 != [] or input_pop2 != [] or input_pop3 != []:
        main_EA_with_loaded_pop(adap_mutation, amount_of_bends, input_pop1, input_pop2, input_pop3, num_gen, num_gen_set2,
                            num_gen_set3, num_gen_set4, p_mutation, pop_size)

    # EA initialized with start chromosom
    else:

        # Repair Startchromosoms to get same amount of bendpoints
        repair_start_chromosoms_to_same_amount_of_bendpoints(use_2D_Solution, use_2D_with_edge_detection,
                                                             use_3D_Solution)
        # Main EA
        main_EA(adap_mutation, amount_of_bends, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_mutation, pop_size,
                startchromo2D, startchromo2D_edge, startchromo3D, use_2D_Solution, use_2D_with_edge_detection,
                use_3D_Solution)

    GUI_End_Save_Patch()

# Functions in main
def repair_start_chromosoms_to_same_amount_of_bendpoints(use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution):
    global startchromo3D, startchromo2D, startchromo2D_edge
    if use_3D_Solution & use_2D_Solution:
        if len(ListOfPoints(startchromo3D)[4]) != len(ListOfPoints(startchromo2D)[4]):
            startchromo3D, startchromo2D = repair_start_chromo(startchromo3D,
                                                                                                startchromo2D)
    if use_2D_with_edge_detection & use_2D_Solution:
        if len(ListOfPoints(startchromo2D_edge)[4]) != len(ListOfPoints(startchromo2D)[4]):
            startchromo2D_edge, startchromo2D = repair_start_chromo(startchromo2D_edge, startchromo2D)
    if use_3D_Solution & use_2D_with_edge_detection:
        if len(ListOfPoints(startchromo3D)[4]) != len(ListOfPoints(startchromo2D_edge)[4]):
            startchromo3D, startchromo2D_edge = repair_start_chromo(startchromo3D, startchromo2D_edge)
def find_and_remove_unnecessary_bendingpoints(chromo):
    betas = ListOfPoints(chromo)[6][0:-1]
    remove_bendpoint_nr = []
    for i in range(len(betas)):
        degree = 5 # remove bendpoints for 5° > beta > -5°
        if betas[i] > (-degree*math.pi/180) and betas[i] < (degree*math.pi/180):
            remove_bendpoint_nr.append(i)
    if remove_bendpoint_nr != []:
        bend_pt_nr_variable = 1 # if a point is deleted, the next Bendpoint is one "position" closer
        for i in remove_bendpoint_nr:
            chromo = remove_bendingpoint(chromo, i+bend_pt_nr_variable)
            bend_pt_nr_variable -= 1
    return chromo, len(remove_bendpoint_nr)
def show_and_save_preprocessor_result(calc_2D_Solution, calc_2D_with_edge_detection, calc_3D_Solution, patch_end, patch_start):

    if calc_3D_Solution:
        stlprep3_6.show_startstrip(ListOfPoints(create_start_chromo(0))[3], patch_start, patch_end, "3D")
        save_start_chromo(str(create_start_chromo(0)), "_3D")
    if calc_2D_Solution:
        stlprep3_6.show_startstrip(ListOfPoints(create_start_chromo(1))[3], patch_start, patch_end, "2D")
        save_start_chromo(str(create_start_chromo(1)), "_2D")
    if calc_2D_with_edge_detection:
        stlprep3_6.show_startstrip(ListOfPoints(create_start_chromo(2))[3], patch_start, patch_end,
                                   "2D with Edge detection")
        save_start_chromo(str(create_start_chromo(2)), "_2DE")
def show_startchromos_results():
    if calc_3D_Solution:
        print_penalties(startchromo3D, "Start solution 3D ")
        show_chromo(startchromo3D, "Start solution 3D ")
    if calc_2D_Solution:
        print_penalties(startchromo2D, "Start solution 2D")
        show_chromo(startchromo2D, "Start solution 2D")
    if calc_2D_with_edge_detection:
        print_penalties(startchromo2D_edge, "Start solution 2DE")
        show_chromo(startchromo2D_edge, "Start solution 2DE")
def print_penalties(start_chromo, start_solution ):
    len_pen, alpha_pen, beta_pen, neg_pen, beta_min_pen = Fitness(start_chromo, 1, True)[6:]
    print("\nProduction - Penalties: length, alpha and beta for " + start_solution +
    "\n\t length penalty:" + str(len_pen) + (
        ". Section are too short. Suggestion: increase distance to geometry in preprocessor " if len_pen < 0.3 else "") +
    "\n\t alpha penalty:" + str(alpha_pen) + (
        ". Alphas are to big. Suggestion: select different start and end point " if alpha_pen < 0.4 else "") +
    "\n\t beta penalty:" + str(beta_pen) + (". Betas are to big, Suggestion: decrease distance to geometry in preprocessor " if beta_pen < 0.4 else "") +

    "\n\t negativ penalty:" + str(neg_pen) +
    "\n\t beta_min penalty:" + str(beta_min_pen)
        )
def if_individual_optimization(individual_optimization, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution,show_solution = True):
    global startchromo3D, startchromo2D_edge, startchromo2D, with_minimize_beta
    if individual_optimization:
        with_minimize_beta = True
        if use_3D_Solution:
            optimize_startchromo_and_remove_bendpoints(startchromo3D,show_solution, "_3D")
        if use_2D_with_edge_detection:
            optimize_startchromo_and_remove_bendpoints(startchromo2D_edge, show_solution, "_2DE")
        if use_2D_Solution:
            optimize_startchromo_and_remove_bendpoints(startchromo2D,show_solution, "_2D")
        with_minimize_beta = False
def optimize_startchromo_and_remove_bendpoints(startchromo3D,show_solution, string_dim):

    startchromo3D = individual_optimization_of_chromo(startchromo3D, string_dim)
    if show_solution: show_chromo(startchromo3D)
    startchromo3D, amount_removed_bends = find_and_remove_unnecessary_bendingpoints(startchromo3D)
    if amount_removed_bends > 0 and show_solution:
        show_chromo(startchromo3D, "Startchromo" + string_dim + " after individual Opti. with " + str(
            amount_removed_bends) + " bendpts removed.")
    save_start_chromo(str(startchromo3D), string_dim)
def individual_optimization_of_chromo(start_chromo, string_dimension):
    #show_chromo(start_chromo, "Start solution "+string_dimension)
    #Initialize and Prep
    p_individual = initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, start_chromo, [], [])
    p_individual.prepPopulation(True, False, False)
    # Save Start_Population
    save_current_population(p_individual, string_dimension+"_Start")
    # EA
    EA_loop(adap_mutation, num_gen, p_individual, p_mutation)
    # Save End_Population
    save_current_population(p_individual, string_dimension+"_after_EA")
    # Update startchromo
    start_chromo = p_individual.currentGeneration[0].genes

    return start_chromo


# Main EA loop
def main_EA(adap_mutation, amount_of_bends, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_mutation, pop_size, startchromo2D, startchromo2D_edge, startchromo3D, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution,with_show_results=True):
    p = initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, startchromo2D, startchromo2D_edge, startchromo3D)
    p.prepPopulation(use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution)
    main_EA_loop(adap_mutation, num_gen, p, p_mutation, with_show_results)
def main_EA_with_loaded_pop(adap_mutation, amount_of_bends, input_pop1, input_pop2, input_pop3, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_mutation,pop_size,with_show_results=True):
    print("Continue with loaded Population")
    use_pop_1, use_pop_2, use_pop_3 = False, False, False
    loaded_pop_1, loaded_pop_2, loaded_pop_3 = [[]], [[]], [[]]
    if input_pop1 != []:
        use_pop_1 = True
        loaded_pop_1 = load_population(pop_size, input_pop1)
    if input_pop2 != []:
        use_pop_2 = True
        loaded_pop_2 = load_population(pop_size, input_pop2)
    if input_pop3 != []:
        use_pop_3 = True
        loaded_pop_3 = load_population(pop_size, input_pop3)

    # Initialize
    p_loaded = initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, loaded_pop_1[0],
                                                          loaded_pop_2[0], loaded_pop_3[0])  # pop_2D[0] is the best solution of that generation
    # Prep population with loaded pop. If more then one selected, just the top 1/2 or 1/3 of every population is taken.
    p_loaded.prepPopulation_read_in_population(use_pop_1, use_pop_2, use_pop_3, loaded_pop_1, loaded_pop_2,
                                               loaded_pop_3)
    main_EA_loop(adap_mutation, num_gen, p_loaded, p_mutation, with_show_results)
def main_EA_loop(adap_mutation, num_gen, p, p_mutation, with_show_results):
    global with_production_restriction
    with_production_restriction = True
    save_current_population(p, "_main_Start")  # Save Start_Population
    EA_loop(adap_mutation, num_gen, p, p_mutation)
    save_current_population(p, "_main_after_EA")  # Save End_Population

    # Save fitness over generation lists
    save_fitness_over_gen()

    # Print and save results
    print_consol_output_end(p)
    if with_show_results:
        show_fitness_and_subfitness_over_generations_end()
        show_chromo(p.bestFitIndividual.genes,"Result of main EA")

    # Remove unnecessary bendingpoints (remove bendpoints for 5° > beta > -5° )  # todo: insert production accuracy
    global best_fit_with_optimized_amount_of_bendingpoints
    best_fit_with_optimized_amount_of_bendingpoints = list(p.bestFitIndividual.genes)
    best_fit_with_optimized_amount_of_bendingpoints, amount_removed_bends = find_and_remove_unnecessary_bendingpoints(
        best_fit_with_optimized_amount_of_bendingpoints)
    if best_fit_with_optimized_amount_of_bendingpoints != p.bestFitIndividual.genes:
        print(str(amount_removed_bends) + " Bendpoints have been removed from the end result.")
        if with_show_results:
            show_chromo(best_fit_with_optimized_amount_of_bendingpoints, "Result after removed bendpoints")
        # save_start_chromo(best_fit_with_removed_bendingpoints,"_removed_bendpoint_result")

# EA Loop
def EA_loop(adap_mutation, num_gen, p, p_mutation):
    global time_start, time_end, fitness_list_gen_index, distance_fit_list_gen_index, length_fit_list_gen_index, border_fit_start_list_gen_index, border_fit_end_list_gen_index, mutation_rate_list_gen_index,penalty_alpha_list_gen_index,penalty_beta_list_gen_index ,penalty_length_list_gen_index, penalty_negativ_list_gen_index,penalty_minimise_beta_list_gen_index
    global border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, mutation_rate_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list

    ##Comment_DB: initialize arrays of the fitness values (Saving values in the arrays)
    border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, mutation_rate_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list = initialize_lists_fitness_over_generation()

    ## Initializazion
    time_start = timer()  # Comment_DB: start timer
    p.currentGeneration.sort()  # Comment_DB: sort and reverse randomly initialized pop from best fit to worst
    p.currentGeneration.reverse()

    # Print settings
    print_statement_beginning(adap_mutation, p)

    # Comment_DB: EA Loop (Only this for loop determines the number of generations to iterate! This is not determined in galileo!)
    for i in range(num_gen):
        print_statement_within_EA_loop_begin(adap_mutation, i, p)

                # Comment_DB: Begin Evolutionary Algorithm
        if i != num_gen - 1:
            # Calculate Fitness for each individum

            p.evaluate() # todo: is that needed? We calculated fitness while replacement

            p_mutation_rate_adding_adaptive = 0.2
            p_difference_avg_fit_to_best_fit = 50

            if adap_mutation == 1 and p.generationNumber > num_gen / 2:  # Comment_DB: Increase/decrease mutation rate by 0.5 for adaptive mutation for mutationrate less than 0.5. 0.999 for > 0.5
                if p.avgFitness >= p.bestFitIndividual.fitness - p_difference_avg_fit_to_best_fit and p.mutationRate == p_mutation:
                    if p_mutation >= 0.8:
                        p.mutationRate = 0.999
                    else:
                        p.mutationRate = p.mutationRate + p_mutation_rate_adding_adaptive
                elif p.avgFitness >= p.bestFitIndividual.fitness - p_difference_avg_fit_to_best_fit and (
                        p.mutationRate == p_mutation + p_mutation_rate_adding_adaptive or p.mutationRate == 0.999):
                    pass
                else:
                    p.mutationRate = p_mutation

            # Selection of parents
            p.select()
            # Recombination
            p.crossover()
            # Mutation
            p.mutate()
            # Replace old with new generation
            p.replace() # DKu_Wenzel: Within replace_func we sort the current generation. (Takes a long time since we calculate here the fitness for sorting)


        # DKu_Wenzel: New bestFitIndividual after crossover and mutation
        p.bestFitIndividual = p.currentGeneration[0]



        # print the best fit individual, and its fitness
        fitness_best_gen = Fitness(p.bestFitIndividual.genes, p.generationNumber)
        print_statement_within_EA_loop_end(fitness_best_gen, i, p)

        p.generationNumber = p.generationNumber + 1  # Comment_DB: This is one plus the output (Gen 0 has generationNumber 1)

        # Comment_DB: append determined values into arrays after each iteration
        border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list,mutation_rate_list = append_to_list_fitness_over_generation(
            border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_best_gen, fitness_list, i,
            length_fit_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list,
            penalty_minimise_beta_list, penalty_negativ_list,mutation_rate_list,p)

    time_end = timer()  # Comment_DB: end timer

    stack_lists_fitness_over_generation(adap_mutation, border_fit_end_list, border_fit_start_list, distance_fit_list,
                                        fitness_list, length_fit_list, mutation_rate_list, num_gen_list,
                                        penalty_alpha_list, penalty_beta_list, penalty_length_list,
                                        penalty_minimise_beta_list, penalty_negativ_list)
def print_statement_within_EA_loop_end(fitness_best_gen, i, p):
    print("\nBest Fit Member of Generation ", i, " :", p.bestFitIndividual, "\n\tFitness:",
          fitness_best_gen[0],
          "\n\t\tdistance Fit:", fitness_best_gen[1],
          "\n\t\tLength Fit:", fitness_best_gen[2], "\n\t\tBorder Fit Start:",
          fitness_best_gen[3],
          "\n\t\tBorder Fit End:", fitness_best_gen[4])
    print("\t\tAverage distance", fitness_best_gen[5])
    print_penalties(p.bestFitIndividual.genes, "Best Fit")
    print("\n")
def print_statement_within_EA_loop_begin(adap_mutation, i, p):
    if adap_mutation == 0:
        print("\n#####Elite Population Members of Generation", i, "\b#####")
    else:
        if not i == 0:
            print("\n#####Elite Population Members of Generation", i, "\b.", "Mutation Rate:", p.mutationRate,
                  "\b#####")
        else:
            print("\n#####Elite Population Members of Generation", i, "\b#####")
    # Comment_DB: print elite population members and their overall fitness, distance fit, and average distance
    for j in range(p.selectionSize):
        curren_gen_fit = Fitness(p.currentGeneration[j].genes, p.generationNumber)
        print("\n\tPopulation Member ", j, " :",
              p.currentGeneration[j].genes,
              "\n\t\tMember Fitness:",
              curren_gen_fit[0],
              "\tMember distance Fit:",
              curren_gen_fit[1],
              "\tMember Average distance:",
              curren_gen_fit[5]
              )
def print_statement_beginning(adap_mutation, p):
    if p.preprocessedInit == False:  # Comment_DB: print if random or preprocessor init
        print("##########Random Gene Initialization for Generation 0##########")
    else:
        print("##########Preprocessor Gene Initialization for Generation 0##########")
    print("##########Generations Sorted from Highest to Lowest Fitness##########")
    if adap_mutation == 1:
        print("##########Adaptive Mutation Selected##########")
# EA Fitness over Generation
def initialize_lists_fitness_over_generation():
    num_gen_list = np.array([])
    fitness_list = np.array([])
    distance_fit_list = np.array([])
    length_fit_list = np.array([])
    border_fit_start_list = np.array([])
    border_fit_end_list = np.array([])
    mutation_rate_list = np.array([])
    penalty_alpha_list = np.array([])
    penalty_beta_list = np.array([])
    penalty_length_list = np.array([])
    penalty_negativ_list = np.array([])
    penalty_minimise_beta_list = np.array([])
    return border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, mutation_rate_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list
def append_to_list_fitness_over_generation(border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_best_gen, fitness_list, i, length_fit_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list, mutation_rate_list, p):
    num_gen_list = np.append(num_gen_list, [i])
    fitness_list = np.append(fitness_list, [fitness_best_gen[0]])
    distance_fit_list = np.append(distance_fit_list, [fitness_best_gen[1]])
    length_fit_list = np.append(length_fit_list, [fitness_best_gen[2]])
    if adap_mutation == 1:
        if i == 0:
            mutation_rate_list = np.append(mutation_rate_list, [0])
        else:
            mutation_rate_list = np.append(mutation_rate_list, [p.mutationRate])
    border_fit_start_list = np.append(border_fit_start_list,
                                      [fitness_best_gen[3]])
    border_fit_end_list = np.append(border_fit_end_list,
                                    [fitness_best_gen[4]])
    penalty_alpha_list = np.append(penalty_alpha_list,
                                   [fitness_best_gen[7]])
    penalty_beta_list = np.append(penalty_beta_list,
                                  [fitness_best_gen[8]])
    penalty_length_list = np.append(penalty_length_list,
                                    [fitness_best_gen[6]])
    penalty_negativ_list = np.append(penalty_negativ_list,
                                     [fitness_best_gen[9]])
    penalty_minimise_beta_list = np.append(penalty_minimise_beta_list,
                                           [fitness_best_gen[10]])
    return border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, num_gen_list, penalty_alpha_list, penalty_beta_list, penalty_length_list, penalty_minimise_beta_list, penalty_negativ_list,mutation_rate_list
def stack_lists_fitness_over_generation(adap_mutation, border_fit_end_list, border_fit_start_list, distance_fit_list, fitness_list, length_fit_list, mutation_rate_list, num_gen_list,penalty_alpha_list, penalty_beta_list, penalty_length_list,penalty_minimise_beta_list, penalty_negativ_list):
    global fitness_list_gen_index, distance_fit_list_gen_index, length_fit_list_gen_index, border_fit_start_list_gen_index, border_fit_end_list_gen_index, mutation_rate_list_gen_index, penalty_alpha_list_gen_index, penalty_beta_list_gen_index, penalty_length_list_gen_index, penalty_negativ_list_gen_index, penalty_minimise_beta_list_gen_index
    fitness_list_gen_index = np.stack((num_gen_list, fitness_list))  # Comment_DB: stack gen_list with list of fitnesses
    distance_fit_list_gen_index = np.stack((num_gen_list, distance_fit_list))
    length_fit_list_gen_index = np.stack((num_gen_list, length_fit_list))
    border_fit_start_list_gen_index = np.stack((num_gen_list, border_fit_start_list))
    border_fit_end_list_gen_index = np.stack((num_gen_list, border_fit_end_list))
    if adap_mutation == 1:
        mutation_rate_list_gen_index = np.stack((num_gen_list, mutation_rate_list))
    penalty_alpha_list_gen_index = np.stack((num_gen_list, penalty_alpha_list))
    penalty_beta_list_gen_index = np.stack((num_gen_list, penalty_beta_list))
    penalty_length_list_gen_index = np.stack((num_gen_list, penalty_length_list))
    penalty_negativ_list_gen_index = np.stack((num_gen_list, penalty_negativ_list))
    penalty_minimise_beta_list_gen_index = np.stack((num_gen_list, penalty_minimise_beta_list))

# Save and print parameters
def print_consol_output_end(p):
    print("\n\nEnd Patch length: ", patch_length_in_mm(p.bestFitIndividual.genes, l_factor),
          "L_Aim (From Preprocessor):", L_aim)
    print("End Fitness: ", p.bestFitIndividual.getFitness(),
          "\n\tEnd distance Fit:", Fitness(p.bestFitIndividual.genes, p.generationNumber)[1],
          "\n\tEnd Length Fit:", Fitness(p.bestFitIndividual.genes, p.generationNumber)[2],
          "\n\tEnd Border Fit Start:", Fitness(p.bestFitIndividual.genes, p.generationNumber)[3],
          "\n\tEnd Border Fit End:", Fitness(p.bestFitIndividual.genes, p.generationNumber)[4])
    if adap_mutation == 1:
        print("\tEnd Mutation Rate: ", p.mutationRate)
    print("\nSettings Used: ")
    print("Set 1 gamma_d: ", gamma_d, "\tSet 2 gamma_d: ", gamma_d2, "\tSet 3 gamma_d: ", gamma_d3, "\tSet 4 gamma_d: ",
          gamma_d4,
          "\nSet 1 gamma_l: ", gamma_l, "\tSet 2 gamma_l: ", gamma_l2, "\tSet 3 gamma_l: ", gamma_l3,
          "\tSet 4 gamma_l: ",
          gamma_l4,
          "\nSet 1 gamma_ps: ", gamma_ps, "\tSet 2 gamma_ps: ", gamma_ps2, "\tSet 3 gamma_ps: ", gamma_ps3,
          "\tSet 4 gamma_ps: ", gamma_ps4,
          "\nSet 1 gamma_pe: ", gamma_pe, "\tSet 2 gamma_pe: ", gamma_pe2, "\tSet 3 gamma_pe: ", gamma_pe3,
          "\tSet 4 gamma_pe: ", gamma_pe4,
          "\n\nSet 2 Gen Start Point: ", num_gen_set2, "\tSet 3 Gen Start Point: ", num_gen_set3,
          "\tSet 4 Gen Start Point: ", num_gen_set4)
    print("\nPop Size: ", pop_size,
          "\n# of Gens: ", num_gen,
          "\nMutation Rate: ", p_mutation, "\tMutation Rate (Final):",
          p.mutationRate if not adap_mutation == 0 else "N/A",
          "\tMutation Range: ", p_mutate_range,
          "\nCrossover Rate: ", p_crossover if not p.crossoverFunc == p.crossover_Flat else "N/A")
    print("\n\nElapsed Time [s]: ", time_end - time_start)
def show_fitness_and_subfitness_over_generations_end():
    ##Comment_DB: Create plots for fitness and subfitnesses
    plt.figure(1)
    plt.plot(fitness_list_gen_index[0], fitness_list_gen_index[1], linewidth=2)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.title('Overall Fitness 1')
    plt.figure(2)
    plt.plot(distance_fit_list_gen_index[0], distance_fit_list_gen_index[1])
    plt.xlabel('Generation')
    plt.ylabel('distance Fitness')
    plt.title('distance Fitness 1')
    plt.figure(3)
    plt.plot(length_fit_list_gen_index[0], length_fit_list_gen_index[1])
    plt.xlabel('Generation')
    plt.ylabel('Length Fitness')
    plt.title('Length Fitness 1')
    plt.figure(4)
    plt.plot(border_fit_start_list_gen_index[0], border_fit_start_list_gen_index[1])
    plt.xlabel('Generation')
    plt.ylabel('Border Fitness at Start')
    plt.title('Border Fitness at Start 1')
    plt.figure(5)
    plt.plot(border_fit_end_list_gen_index[0], border_fit_end_list_gen_index[1])
    plt.xlabel('Generation')
    plt.ylabel('Border Fitness at End')
    plt.title('Border Fitness at End 1')
    if adap_mutation == 1:
        plt.figure(6)
        plt.stem(mutation_rate_list_gen_index[0], mutation_rate_list_gen_index[1])
        plt.xlabel('Generation')
        plt.ylabel('Mutation Rate')
        plt.title('Mutation Rate')
    plt.show()
def GUI_End_Save_Patch():
    global end
    end = Tk()
    Label(end, text="Do you want to save the result?").grid(row=10, column=1, )
    Label(end, justify=LEFT, text=" ").grid(row=11, sticky=W)
    Button(end, text="End", command=sys.exit).grid(row=30, column=0, )
    Button(end, text="Save used settings", command=save_used_settings).grid(row=30, column=1, )
    Button(end, text="Save patch parameters", command=save_patch_file).grid(row=30, column=2, )

    Label(end, justify=LEFT, text=" ").grid(row=11, sticky=W)



    mainloop()
def save_patch_file():
    name = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    bestPatch_parameter_l = ListOfPoints(best_fit_with_optimized_amount_of_bendingpoints)[4]
    bestPatch_parameter_alpha = ListOfPoints(best_fit_with_optimized_amount_of_bendingpoints)[5]
    bestPatch_parameter_beta = ListOfPoints(best_fit_with_optimized_amount_of_bendingpoints)[6]
    name.write("width=" + str(width) + "\n")
    name.write("length=" + str(sum(bestPatch_parameter_l)) + "\n")
    name.write("type=" + str(tape_type) + "\n")
    l = 0
    for i in range(amount_of_bends):
        l = int(l) + int(bestPatch_parameter_l[i])
        l = str(l)
        if bestPatch_parameter_alpha[i] > 90 * 2 * math.pi / (360):
            alpha = str(int((bestPatch_parameter_alpha[i] - math.pi) * 360 / (2 * math.pi)))
        else:
            alpha = str(int(bestPatch_parameter_alpha[i] * 360 / (2 * math.pi)))
        beta = str(int(bestPatch_parameter_beta[i] * 360 / (2 * math.pi)))
        name.write(l + ";" + alpha + ";" + beta + "\n")
    name.close
    end.destroy()
def save_used_settings():
    """
    Save settingssheet, startpara and settinssheet_EA to separate folder

    Also save used populations into that folder.

    Populations just make sense with their settingssheets!!
    """

    # get the current script path and go to subdir where the results are getting saved
    here = os.path.dirname(os.path.realpath(__file__))
    subdir = "settings_population"
    subdir_path = os.path.join(here, subdir)

    if not os.path.isdir("./settings_population"):
        os.mkdir(os.path.join(here, subdir))

    create_subdir_with_results(subdir_path)
def create_subdir_with_results(subdir_path):
    # Name of the new folder. Date at the begin and time at the end for uniqunes.
    use2D, use2DE, use3D, with_individual_optimization = "", "", "", ""
    if use_2D_with_edge_detection: use2DE = "_2DE"
    if use_2D_Solution: use2D = "_2D"
    if use_3D_Solution: use3D = "_3D"
    if individual_optimization: with_individual_optimization = "_indOpt"
    current_pop_settings = str(datetime.date.today()) + time.strftime("_%H_%M_%S", time.localtime()) + "_gen" + str(
        num_gen) + "_pop" + str(
        pop_size) + with_individual_optimization + use2D + use2DE + use3D
    # Creat the new folder
    os.mkdir(os.path.join(subdir_path, current_pop_settings))
    # Patch to new folder
    path_to_settings = os.path.join(subdir_path, current_pop_settings)
    # copy settings in new folder
    shutil.copy('settingssheet_EA.txt', path_to_settings)
    shutil.copy('settingssheet.txt', path_to_settings)
    shutil.copy('startparameter.txt', path_to_settings)
    if use_2D_Solution:
        shutil.copy('start_chromo_2D.txt', path_to_settings)
        if individual_optimization:
            shutil.copy('population_2D_Start.txt', path_to_settings)
            shutil.copy('population_2D_after_EA.txt', path_to_settings)
    if use_2D_with_edge_detection:
        shutil.copy('start_chromo_2DE.txt', path_to_settings)
        if individual_optimization:
            shutil.copy('population_2DE_Start.txt', path_to_settings)
            shutil.copy('population_2DE_after_EA.txt', path_to_settings)
    if use_3D_Solution:
        shutil.copy('start_chromo_3D.txt', path_to_settings)
        if individual_optimization:
            shutil.copy('population_3D_Start.txt', path_to_settings)
            shutil.copy('population_3D_after_EA.txt', path_to_settings)
    shutil.copy('population_main_Start.txt', path_to_settings)
    shutil.copy('population_main_after_EA.txt', path_to_settings)
    shutil.copy('fitness_over_generation.txt', path_to_settings)

chromo_resolution = 2000  # Default value
if __name__ == '__main__':
    main()
else:
    global turn_around_normal, equidistant_pts_between_bendpts, step_size, pointspersection, with_production_restriction, with_minimize_beta
    turn_around_normal = False
    with_production_restriction = False
    with_minimize_beta = False
    equidistant_pts_between_bendpts = True
    step_size = 1
    pointspersection = 10



