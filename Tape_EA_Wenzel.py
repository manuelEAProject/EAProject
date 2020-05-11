import numpy as np
import inspect
import sys
import os
import pyquaternion
from pyquaternion import Quaternion
import math
import random
import math
import trimesh
from stl import mesh
import timeit
from trimesh import proximity
from galileo_EA_Edit6 import Population
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import matplotlib.ticker as ticker
import stl_preprocessing_Edit3_Tape_EA_Edit6 as stlprep3_6
from tkinter import *
from tkinter import filedialog
import plotly.graph_objects as go
import xlsxwriter
from timeit import default_timer as timer

############## GUI  &  SETTINGS #######################
# GUI-Settings

def callback(input_file):
    if os.path.isfile('./settingssheet.txt'):
        settingssheet = open('./settingssheet.txt')
        name = filedialog.askopenfilename(initialdir=settingssheet.readline(), title="Select file",
                                          filetypes=(("stl files", "*.stl"), ("all files", "*.*")))
    else:
        name = filedialog.askopenfilename(initialdir="/", title="Select file",
                                          filetypes=(("stl files", "*.stl"), ("all files", "*.*")))
    input_file.delete(0, 'end')
    input_file.insert(0, name)
    if len(input_file.get()) < 1:
        input_file.insert(0, "Datei auswählen...")

    settingssheet = open('./settingssheet.txt', 'w+')
    settingssheet.write(input_file.get())
    settingssheet.close()


def save_settings(settings_list):
    settings_sheet = open('./settingssheet.txt', 'r+')

    for listitem in settings_list:
        settings_sheet.write('%s\n' % listitem)
    settings_sheet.close()


master = Tk()
master.protocol("WM_DELETE_WINDOW", sys.exit)

Label(master, text="Settings für Tape - Algorithmus").grid(row=10, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=11, sticky=W)

input_file = Entry(master)

if os.path.isfile('./settingssheet.txt'):
    settingssheet = open('./settingssheet.txt')
    input_file.insert(0, settingssheet.readline())
    settingssheet.close()
else:
    input_file.insert(0, "Datei auswählen...")

input_file.grid(row=20, columnspan=3, sticky=W + E + N + S)
Button(text='Browse..', command=lambda: callback(input_file)).grid(row=20, column=3, sticky=W)

Label(master, text="Verwendung des Preprozess oder Randominitialisierung").grid(row=30, sticky=W)
init_preprocess = BooleanVar()
init_preprocess.set(True)
init_preprocess_cb = Checkbutton(master, text="Preprozess", variable=init_preprocess,
                                 command=lambda: init_random.set(False))
init_preprocess_cb.grid(row=31, column=0, sticky=W)
init_random = BooleanVar()
init_random_cb = Checkbutton(master, text="Random", variable=init_random, command=lambda: init_preprocess.set(False))
init_random_cb.grid(row=31, column=1, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=31, sticky=W)

prep_start_end = BooleanVar()
prep_start_end.set(False)  # Comment_DB: unchecked as default
prep_start_end_cb = Checkbutton(master, text="Start/Endpunkt aus Präprozess", variable=prep_start_end,
                                command=lambda: manual_start_end.set(False))
prep_start_end_cb.grid(row=33, column=0, sticky=W)
manual_start_end = BooleanVar()
manual_start_end.set(False)
manual_start_end_cb = Checkbutton(master, text="Start/Endpunkt manuell festlegen", variable=manual_start_end,
                                  command=lambda: prep_start_end.set(False))
manual_start_end_cb.grid(row=33, column=1, sticky=W)

Label(master, text="Startpunkt [mm]").grid(row=34, column=2, sticky=W)
Label(master, text="Endpunkt [mm]").grid(row=34, column=3, sticky=W)
Label(master, text="x: ").grid(row=35, column=1, sticky=E)
Label(master, text="y: ").grid(row=36, column=1, sticky=E)
Label(master, text="z: ").grid(row=37, column=1, sticky=E)
x_start = Entry(master)
x_start.grid(row=35, column=2, sticky=W)
y_start = Entry(master)
y_start.grid(row=36, column=2, sticky=W)
z_start = Entry(master)
z_start.grid(row=37, column=2, sticky=W)
x_end = Entry(master)
x_end.grid(row=35, column=3, sticky=W)
y_end = Entry(master)
y_end.grid(row=36, column=3, sticky=W)
z_end = Entry(master)
z_end.grid(row=37, column=3, sticky=W)

Label(master, text="Savitkzy-Golay").grid(row=42, sticky=W)
Label(master, text="Ordnungsgrad").grid(row=43, sticky=W)
poly_order = Entry(master)
poly_order.insert(0, 3)
poly_order.grid(row=43, column=1, sticky=W)
Label(master, text="Fensterquotient").grid(row=43, column=2, sticky=W)
window_quotient = Entry(master)
window_quotient.insert(0, 5)
window_quotient.grid(row=43, column=3, sticky=W)
Label(master, text='Maximaler Abstand zur SavGol-Kurve (Zunaechst)').grid(row=44, sticky=W)
max_distance = Entry(master)
max_distance.insert(0, 10)  # Comment_DB: insert 10 at front of list
max_distance.grid(row=44, column=1, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=45, sticky=W)

Label(master, text="Tape-Typ :").grid(row=46, sticky=W)
tape_type = Entry(master)
tape_type.insert(0, "BASF_CFK_Tape")
tape_type.grid(row=46, column=1, sticky=W)
Label(master, text="Tapebreite [mm]:").grid(row=47, sticky=W)
width = Entry(master)
width.insert(0, 20)
width.grid(row=47, column=1, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=48, sticky=W)

Label(master, text="Punktauflösung Patch").grid(row=55, sticky=W)

fix_number_of_pts = BooleanVar()
fix_number_of_pts_cb = Checkbutton(master, text="Fixe Punktanzahl zwischen Biegestellen", variable=fix_number_of_pts,
                                   command=lambda: equidistant_pts_between_bendpts.set(False))
fix_number_of_pts_cb.grid(row=60, column=0, sticky=W)
Label(master, text="Anzahl Punkte:").grid(row=60, column=1, sticky=W)
fix_number_of_pts.set(True)
pointspersection = Entry(master)
pointspersection.insert(0, 7)
pointspersection.grid(row=60, column=2, sticky=W)

equidistant_pts_between_bendpts = BooleanVar()
equidistant_pts_between_bendpts_cb = Checkbutton(master, text="Äquidistante Punkte",
                                                 variable=equidistant_pts_between_bendpts,
                                                 command=lambda: fix_number_of_pts.set(False),
                                                 disabledforeground="gray")
equidistant_pts_between_bendpts_cb.grid(row=70, column=0, sticky=W)
Label(master, text="Abstand Punkte [mm]:").grid(row=70, column=1, sticky=W)
step_size = Entry(master)
step_size.insert(0, 6)
step_size.grid(row=70, column=2, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=71, sticky=W)

#####Comment_DB: Different weightings for 4 equally-divided Population Sets!#####
Label(master, text="Feintuning Fitnessfunktion (4 Population Sets)").grid(row=73, sticky=W)
Label(master, text="Gewichtung Abstandsfitness [gamma_d]: ").grid(row=74, column=0, sticky=W)
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

Label(master, text="Gewichtung Längenfitness [gamma_l] : ").grid(row=75, column=0, sticky=W)
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

Label(master, text="Gewichtung Startfitness [gamma_ps] : ").grid(row=76, column=0, sticky=W)
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

Label(master, text="Gewichtung Endfitness [gamma_pe] : ").grid(row=77, column=0, sticky=W)
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
num_gen_set2.grid(row=78, column=2)

num_gen_set3 = Entry(master)
num_gen_set3.grid(row=78, column=3)

num_gen_set4 = Entry(master)
num_gen_set4.grid(row=78, column=4)

Label(master, justify=LEFT, text=" ").grid(row=79, sticky=W)
Label(master, text="Einstellungen für den Evolutionären Algorithmus").grid(row=80, sticky=W)

Label(master, justify=LEFT, text="Populationsgröße:").grid(row=90, column=0, sticky=W)
pop_size = Entry(master)
pop_size.insert(0, 12)
pop_size.grid(row=90, column=1, sticky=W)

Label(master, text="Anzahl Generationen:").grid(row=100, column=0, sticky=W)
num_gen = Entry(master)
num_gen.insert(0, 80)
num_gen.grid(row=100, column=1, sticky=W)

Label(master, text="Auflösung der Allelparameter:").grid(row=110, column=0, sticky=W)
chromo_resolution = Entry(master)
chromo_resolution.insert(0, 100)
chromo_resolution.grid(row=110, column=1, sticky=W)

Label(master, text="Mutation Rate:").grid(row=120, sticky=W)
p_mutation = Entry(master)
p_mutation.insert(0, 0.1)
p_mutation.grid(row=120, column=1, sticky=W)
Label(master, text="Mutationsrange:").grid(row=120, column=2, sticky=W)
p_mutate_range = Entry(master)
p_mutate_range.insert(0, 0.05)
p_mutate_range.grid(row=120, column=3, sticky=W)

adap_mutation = BooleanVar()
adap_mutation.set(False)
adap_mutation_cb = Checkbutton(master, text="Adaptive Mutation", variable=adap_mutation)
adap_mutation_cb.grid(row=110, column=3, sticky=W)

Label(master, text="Crossoverrate:").grid(row=130, sticky=W)
p_crossover = Entry(master)
p_crossover.insert(0, 0.7)
p_crossover.grid(row=130, column=1, sticky=W)

if os.path.isfile('./settingssheet.txt'):

    t = "True" + '\n'
    try:
        settingssheet = open('./settingssheet.txt')
        if not settingssheet.readline() == "Datei auswählen...":
            settingssheet.seek(0)
            input_file.delete(0, 'end')
            input = settingssheet.readline()
            input_file.insert(0, input[0:-1])
            if settingssheet.readline() == t:
                init_preprocess.set(True)
            else:
                init_preprocess.set(False)

            if settingssheet.readline() == t:
                init_random.set(True)
            else:
                init_random.set(False)

            width.delete(0, 'end')
            width.insert(0, float(settingssheet.readline()))
            if settingssheet.readline() == t:
                fix_number_of_pts.set(True)
            else:
                fix_number_of_pts.set(False)

            pointspersection.delete(0, 'end')
            pointspersection.insert(0, int(settingssheet.readline()))
            if settingssheet.readline() == t:
                equidistant_pts_between_bendpts.set(True)
            else:
                equidistant_pts_between_bendpts.set(False)

            step_size.delete(0, 'end')
            step_size.insert(0, float(settingssheet.readline()))

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

            num_gen_set2.delete(0, 'end')
            num_gen_set2.insert(0, float(settingssheet.readline()))
            num_gen_set3.delete(0, 'end')
            num_gen_set3.insert(0, float(settingssheet.readline()))
            num_gen_set4.delete(0, 'end')
            num_gen_set4.insert(0, float(settingssheet.readline()))

            pop_size.delete(0, 'end')
            pop_size.insert(0, int(settingssheet.readline()))
            num_gen.delete(0, 'end')
            num_gen.insert(0, int(settingssheet.readline()))
            chromo_resolution.delete(0, 'end')
            chromo_resolution.insert(0, int(settingssheet.readline()))
            p_mutation.delete(0, 'end')
            p_mutation.insert(0, float(settingssheet.readline()))

            if settingssheet.readline() == t:
                adap_mutation.set(True)
            else:
                adap_mutation.set(False)

            p_mutate_range.delete(0, 'end')
            p_mutate_range.insert(0, float(settingssheet.readline()))
            p_crossover.delete(0, 'end')
            p_crossover.insert(0, float(settingssheet.readline()))
            poly_order.delete(0, 'end')
            poly_order.insert(0, int(settingssheet.readline()))
            window_quotient.delete(0, 'end')
            window_quotient.insert(0, int(settingssheet.readline()))
            max_distance.delete(0, 'end')
            max_distance.insert(0, int(settingssheet.readline()))
            settingssheet.close()
    except:
        print("Bitte settingssheet.txt löschen")
        settingssheet.close()

Button(master, text='Abbrechen', command=sys.exit).grid(row=1000, column=0, pady=4)
Button(master, text='Start', command=master.quit).grid(row=1000, column=1, pady=4)
mainloop()  # führt das GUI aus

input_file = input_file.get()
testpatch = trimesh.load(input_file)

tape_type = tape_type.get()
width = float(width.get())  # Tapebreite
# Sollen die Punkte zwischen den Biegestellen gleichmäßigverteilt werden (True,-> step_size) oder eine fixe Anzahl an Punkten
# Zwischen den Biegestellen gewählt werden (False,-> pointspersection)
equidistant_pts_between_bendpts = float(equidistant_pts_between_bendpts.get())
step_size = float(step_size.get())  # Abstand zwischen den Patchpunkten - NUR FALLS equidist_pts_between_bendpts=True

pointspersection = int(pointspersection.get())  # Anzahl der Punkte zwischen den Biegestellen

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



manual_start_end = bool(manual_start_end.get())  # Comment_DB: "Start/Endpunkt aus Praeprozess"
if manual_start_end:  # Comment_DB: If true
    x_start = float(x_start.get())
    y_start = float(y_start.get())
    z_start = float(z_start.get())
    x_end = float(x_end.get())
    y_end = float(y_end.get())
    z_end = float(z_end.get())


#### Settings für Evolutionären Algorithmus ####
num_gen = int(num_gen.get())  # Anzahl der Generationen       Comment_DB: User inputs
pop_size = int(pop_size.get())  # Populationsgröße (Anzahl Lösungen pro Generation)
chromo_resolution = int(chromo_resolution.get())  # Auflösung der Allelparameter
# Soll die Initialisierung mit den Startwerten des Preprozesses erfolgen? 1 = True
init_preprocess = bool(init_preprocess.get())
poly_order = int(poly_order.get())
window_quotient = int(window_quotient.get())
max_distance = int(max_distance.get())
# Ganzzahlige Allele -> useInteger = 1, Floatwerte -> useInteger = 0
useInteger = 1
# Crossover Wahrscheinlichkeit     Comment_DB: Probability nature not defined yet
p_crossover = float(p_crossover.get())
# Mutations Wahrscheinlichkeit
p_mutation = float(p_mutation.get())
###Comment_DB: Adaptive Mutation###
adap_mutation = float(adap_mutation.get())
# Mutationsrange (in welchem Prozentbereich darf das Allel mutieren?)
p_mutate_range = float(p_mutate_range.get())

save_settings([input_file, init_preprocess, init_random.get(), width, fix_number_of_pts.get(),
                 pointspersection, equidistant_pts_between_bendpts, step_size,
                 gamma_d, gamma_d2, gamma_d3, gamma_d4,
                 gamma_l, gamma_l2, gamma_l3, gamma_l4,
                 gamma_ps, gamma_ps2, gamma_ps3, gamma_ps4,
                 gamma_pe, gamma_pe2, gamma_pe3, gamma_pe4,
                 num_gen_set2, num_gen_set3, num_gen_set4,
                 pop_size, num_gen, chromo_resolution,
                 p_mutation, adap_mutation, p_mutate_range, p_crossover, poly_order,
                 window_quotient, max_distance])  # Speichert die gewählten Einstellungen

master.destroy()  # Schließt das Settings Fenster

############Vorverarbeitung der Geometriedaten###################

# Ruft das stl_preprocessing modul auf und übergibt die stl-Datei an die Funktion startparam
# Startparam gibt eine Liste mit den berechneten Startparametern zurück:
# STARTPARAMETER: [Start_p_mid, Start_r_mid, Start_n_mid, l_list, totallength, beta_list[°], Start_p_ID]
start_parameter = stlprep3_6.startparam(input_file, poly_order, window_quotient, max_distance)

Start_p_prep = start_parameter[0] # Comment_DKu_Wenzel: No usage anymore!!
Start_r_prep = start_parameter[1] # Comment_DKu_Wenzel: No usage anymore!!
Start_n = start_parameter[2] # Comment_DKu_Wenzel: No usage anymore!!
Start_p_id = start_parameter[6]  # Comment_DKu_Wenzel: No usage anymore!! # Comment_DB: The I.D. of the next bending point from start p mid in trendline direction
L_aim = start_parameter[4]  # Comment_DB: already in [mm]
patch_start = start_parameter[7]
patch_end = start_parameter[8]
AnzahlKnicke = len(start_parameter[3]) - 1

# COMMENT_DB: New parameters
Start_r_prep_fromstart = start_parameter[10]
Start_n_atstart = start_parameter[11]
Start_p_id_fromstart = start_parameter[9]
AnzahlKnicke_too = start_parameter[12]  # Comment_DKu_Wenzel: Verschieden von AnzahlKnicke(+2) und keine Benutzung

if manual_start_end:
    patch_start = np.asarray([x_start, y_start, z_start])
    patch_end = np.asarray([x_end, y_end, z_end])

# Faktor für das Längenallel in den Chromosomen -> eine Länge kann maximal L_aim lang werden
l_factor = 0.5 * L_aim / chromo_resolution  # Comment_DB: already in [mm]


# Chromosom der Startlösung (Comment_DB: independent of chromominmaxvalue limitation from galileo)
def preprocessed_chromo(startparameter, chromoresolution):
    # Nimmt die Startparameter aus der Geometriedatenvorverarbeitung und wandelt diese in ein Chromosom mit der
    # entsprechenden Auflösung um. Rückgabewert ist das Chromosom der Startlösung.
    preprocessed_chromo = []

    # Lengths
    for i in range(len(startparameter[3])):
        preprocessed_chromo.append(int(startparameter[3][i] / l_factor))
        if i < len(startparameter[5]):  # Comment_DB: range of beta_list compared to range of l_list is smaller by 1
            preprocessed_chromo.append(int(chromoresolution / 2))  # Comment_DB: Alphas -> zero on default
            beta = start_parameter[5][i]  # Comment_DB: Betas
            beta_chromo = (beta + 90) * chromoresolution / 180
            preprocessed_chromo.append(int(beta_chromo))

    # Variable Startparameter werden standardmäßig auf chromo_resolution/2 gesetzt
    for i in range(7):
        preprocessed_chromo.append(int(chromo_resolution / 2))
    return preprocessed_chromo


# Kinematische Beschreibung des Patchs   COMMENT_DB: This is the translation of values suitable for the evolutionary algorithm!
def ListOfPoints(chromo):  # Comment_DB: chromo not defined elsewhere. chromo here is a parameter. Function definition.
    l_list = []  # Comment_DB: empty list
    alpha_list = []  # Comment_DB: empty list
    beta_list = []  # Comment_DB: empty list
    for i in range(0, len(startchromo) - 5, 3):  # Comment_DB: adjusted for reordered startchromo (lengths)
        if chromo[i] == 0:  # Comment_DB: chromo[i] is the "c" variable!
            l_list.append(1)
        else:
            l_list.append(chromo[i] * l_factor)

    for i in range(1, len(startchromo) - 4, 3):  # Comment_DB: adjusted for reordered startchromo (alphas)
        alpha = 0
        if chromo[i] < chromo_resolution / 2:
            alpha = (135 + (chromo[i] * 45 / (chromo_resolution / 2))) * 2 * math.pi / 360
        else:
            alpha = ((chromo[i] - chromo_resolution / 2) * 45 / (
                    chromo_resolution / 2)) * 2 * math.pi / 360  # Quadratische Vert. von 135°-180°
        alpha_list.append(alpha)

    for i in range(2, len(startchromo) - 3, 3):  # Comment_DB: adjusted for reordered startchromo (betas)
        beta = (chromo[i] * (180 / chromo_resolution) - 90) * 2 * math.pi / 360
        beta_list.append(beta)
    # print("beta_list in LoP", beta_list) #Comment_DB: compare with preprocessor beta_list. LoP one is in radians.
    # Variabler Startpunkt
    var_range = 0.8
    var_start_pt_x = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[
        3 * AnzahlKnicke + 1]  # Comment_DB: chromo[i] is the "c" variable! Do not have to change this w.r.t new order of starting chromosome
    var_start_pt_y = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * AnzahlKnicke + 2]
    var_start_pt_z = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * AnzahlKnicke + 3]

    var_start_r_x = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * AnzahlKnicke + 4]
    var_start_r_y = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * AnzahlKnicke + 5]
    var_start_r_z = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * AnzahlKnicke + 6]

    gamma_max = 10  # [Grad] Maximaler Kippwinkel für Start_n
    gamma_max_rad = gamma_max * (2 * math.pi / 360)
    var_start_n_gamma = -gamma_max_rad + gamma_max_rad / (chromo_resolution / 2) * chromo[3 * AnzahlKnicke + 7]

    Start_p = np.concatenate(np.array(
        [[patch_start[0] * var_start_pt_x], [patch_start[1] * var_start_pt_y],
         [patch_start[2] * var_start_pt_z]]))  # Comment_DB: AMENDED

    Start_r = np.concatenate(np.array(
        [[Start_r_prep_fromstart[0] * var_start_r_x], [Start_r_prep_fromstart[1] * var_start_r_y],
         [Start_r_prep_fromstart[2] * var_start_r_z]]))

    Start_r = 1 / np.linalg.norm(Start_r) * Start_r

    Start_q = np.cross(Start_r, Start_n_atstart)
    Start_q = 1 / np.linalg.norm(Start_q) * Start_q

    Start_n_strich = np.cross(Start_q, Start_r)

    Start_n_strich = Quaternion(axis=Start_r, angle=(var_start_n_gamma)).rotate(
        Start_n_strich)  # Comment_DB: start_n_strich rotated about start_r
    Start_n_strich = 1 / np.linalg.norm(Start_n_strich) * Start_n_strich

    Start_v = np.array([0, 0, 0])

    r_list = [Start_r]
    v_list = [Start_v]
    n_list = [Start_n_strich]

    #### Vektorenberechnung für Tapeseite nach dem Startpunkt
    if not len(alpha_list) == 0:
        for i in range(1, len(l_list)):
            if alpha_list[i - 1] < math.pi / 2:
                v_new = Quaternion(axis=n_list[i - 1], angle=(alpha_list[i - 1] - (math.pi) / 2)).rotate(
                    r_list[i - 1])
            else:
                v_new = Quaternion(axis=n_list[i - 1], angle=(alpha_list[i - 1] - 3 * (math.pi) / 2)).rotate(
                    r_list[i - 1])
            v_list.append(v_new)
            r_rotated = Quaternion(axis=v_list[i], angle=beta_list[i - 1]).rotate(r_list[i - 1])
            r_list.append(r_rotated)
            n_new = Quaternion(axis=v_list[i], angle=beta_list[i - 1]).rotate(n_list[i - 1])
            n_list.append(n_new)
    else:
        print("lenalpha0")

    # print("unstacked r_list", r_list) #Comment_DB: TEST
    r_list = np.stack(r_list)
    # print("stacked r_list", r_list) #Comment_DB: TEST
    ## Mittellinie after startpoint
    p_list = [Start_p]
    old_p = Start_p

    for i, l in enumerate(l_list):
        dist = l * r_list[i, :]
        new_p = old_p + dist
        p_list.append(new_p)
        old_p = new_p

    # Auffüllen von Punkten zwischen den Biegestellen auf Mittellinie. Entweder fixe Anzahl an Punkten oder äquidistant
    if equidistant_pts_between_bendpts:
        middle_pts = []
        for i, l in enumerate(l_list):
            a_new = p_list[i][np.newaxis, :] + np.outer(
                np.linspace(0, l, math.floor(l / step_size), endpoint=True),
                r_list[i])
            middle_pts.append(a_new)
    else:
        middle_pts = []
        for i, l in enumerate(l_list):
            a_new = p_list[i][np.newaxis, :] + np.outer(np.linspace(0, l, pointspersection, endpoint=True),
                                                        r_list[i])
            middle_pts.append(a_new)

    middle_pts = np.concatenate(middle_pts)

    ########## Linker Rand nach Startpunkt ######################
    ## Anpassung Startpunkt Seitenrand links
    alpha_start = alpha_list[Start_p_id_fromstart - 1]  # Comment_DB: Might not be necessary
    if alpha_start > math.pi / 2:
        delta_l_l_start = (width / 2) * math.tan(math.pi - alpha_start)
    if alpha_start < math.pi / 2:
        delta_l_l_start = - (width / 2) * math.tan(alpha_start)

    # Linker Rand
    delta_l_list = [delta_l_l_start]
    l_left_list = []
    if len(alpha_list) == 0:
        l_left_list = [l_list[0]]
    for i in range(1, len(l_list)):
        if alpha_list[i - 1] > math.pi / 2:
            delta_l_l = (width / 2) * math.tan(math.pi - alpha_list[i - 1])
            l_left_new = l_list[i - 1] + delta_l_l - delta_l_list[i - 1]
            l_left_list.append(l_left_new)
            delta_l_list.append(delta_l_l)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_l
            if i == len(l_list) - 1:
                l_left_list.append((l_list[-1] + delta_l_end))

        if alpha_list[i - 1] < math.pi / 2:
            delta_l_l = - (width / 2) * math.tan(alpha_list[i - 1])
            l_left_new = l_list[i - 1] + delta_l_l - delta_l_list[i - 1]
            l_left_list.append(l_left_new)
            delta_l_list.append(delta_l_l)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_l
            if i == len(l_list) - 1:
                l_left_list.append((l_list[-1] + delta_l_end))

    # Ausnahme für keinen Knick nach Start:
    if len(l_list) == 1:
        l_left_list.append((l_list[-1]) + delta_l_l_start)

    # Eckpunkte Left (COMMENT_DB: Misleading comment. This is the modeling of left edge bend points, just like the middle_pts portion above)
    Start_p_left = Start_p - np.cross(Start_r,
                                      Start_n_strich) * width / 2 + delta_l_l_start * Start_r  # Comment_DB: np.cross(Start_r, Start_n_strich) == -Start_q rotated
    p_left_list = [Start_p_left]
    old_p_left = Start_p_left
    for i, l in enumerate(l_left_list):
        dist = l * r_list[i, :]
        new_p_left = old_p_left + dist
        p_left_list.append(new_p_left)
        old_p_left = new_p_left

    # Auffüllen mit Punkten ( Fixe Punktanzahl oder fixer Punktabstand):
    left_pts = []
    if equidistant_pts_between_bendpts:

        for i, l in enumerate(l_left_list):
            if l < 0:
                continue
            else:
                a_new_l = p_left_list[i][np.newaxis, :] + np.outer(np.linspace(0, l, math.floor(l / step_size),
                                                                               endpoint=False), r_list[i])
                left_pts.append(a_new_l)
        if len(left_pts) == 0:
            left_pts.append(Start_p)

    else:

        for i, l in enumerate(l_left_list):
            if l < 0:
                continue
            else:
                a_new_l = p_left_list[i][np.newaxis, :] + np.outer(np.linspace(0, l, pointspersection,
                                                                               endpoint=False), r_list[i])
                left_pts.append(a_new_l)
        if len(left_pts) == 0:
            left_pts.append(Start_p)

    left_pts = np.concatenate(left_pts)
    ######################## Rechter Rand nach Startpunkt######################
    ## Anpassung Startpunkt Seitenrand rechts
    alpha_start = alpha_list[Start_p_id_fromstart - 1]  # Comment_DB: Might not be necessary
    if alpha_start > math.pi / 2:
        delta_l_r_start = -(width / 2) * math.tan(math.pi - alpha_start)
    if alpha_start < math.pi / 2:
        delta_l_r_start = (width / 2) * math.tan(alpha_start)

    delta_r_list = [delta_l_r_start]
    l_right_list = []
    if len(alpha_list) == 0:
        l_right_list = [l_list[0]]
    for i in range(1, len(l_list)):
        if alpha_list[i - 1] > math.pi / 2:
            delta_l_r = - (width / 2) * math.tan(math.pi - alpha_list[i - 1])
            l_right_new = l_list[i - 1] + delta_l_r - delta_r_list[i - 1]
            l_right_list.append(l_right_new)
            delta_r_list.append(delta_l_r)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_r
            if i == len(l_list) - 1:
                l_right_list.append((l_list[-1] + delta_l_end))

        if alpha_list[i - 1] < math.pi / 2:
            delta_l_r = (width / 2) * math.tan(alpha_list[i - 1])
            l_right_new = l_list[i - 1] + delta_l_r - delta_r_list[i - 1]
            l_right_list.append(l_right_new)
            delta_r_list.append(delta_l_r)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_r
            if i == len(l_list) - 1:
                l_right_list.append((l_list[-1] + delta_l_end))

    # Ausnahme für keinen Knick nach Start:
    if len(l_list) == 1:
        l_right_list.append((l_list[-1]) + delta_l_r_start)
    # print("längen", l_left_list_a,'\n',l_list_a,'\n', l_right_list_a)

    # Eckpunkte Rechts
    Start_p_right = Start_p + np.cross(Start_r, Start_n_strich) * width / 2 + delta_l_r_start * Start_r
    p_right_list = [Start_p_right]
    old_p_right = Start_p_right
    for i, l in enumerate(l_right_list):
        dist = l * r_list[i, :]
        new_p_right = old_p_right + dist
        p_right_list.append(new_p_right)
        old_p_right = new_p_right

    right_pts = []
    if equidistant_pts_between_bendpts:

        for i, l in enumerate(l_right_list):
            if l < 0:
                continue
            else:
                a_new_r = p_right_list[i][np.newaxis, :] + np.outer(
                    np.linspace(0, l, math.floor(l / step_size), endpoint=False), r_list[i])
                right_pts.append(a_new_r)
        if len(right_pts) == 0:
            right_pts.append(Start_p)
    else:

        for i, l in enumerate(l_right_list):
            if l < 0:
                continue
            else:
                a_new_r = p_right_list[i][np.newaxis, :] + np.outer(
                    np.linspace(0, l, pointspersection, endpoint=False), r_list[i])
                right_pts.append(a_new_r)
        if len(right_pts) == 0:
            right_pts.append(Start_p)

    right_pts = np.concatenate(right_pts)

    result = np.concatenate((middle_pts, left_pts, right_pts), axis=0)
    start = middle_pts[0]  # Comment_DB: Change to beginning of list
    end = middle_pts[-1]

    # Nur Eck- und Biegepunkte des Tapes für einfache Visualisierung
    patch_visualisation_points = []

    for i in range(len(p_left_list)):
        patch_visualisation_points.append(p_left_list[i])
        patch_visualisation_points.append(p_right_list[i])

    patch_visualisation_points = np.stack(patch_visualisation_points, axis=0)

    return result, start, end, patch_visualisation_points, l_list, alpha_list, beta_list, Start_p, Start_r  # Comment_DB: Not dependent on preprocessed_chromo


# Berechnung der Patchlänge (für Fitness)
def PatchLength(chromo, AnzahlKnicke, l_factor):
    # Berechnet die Länge eines Patches. Kann Chomosome als class chromosome oder auch als einfache Liste auslesen.

    if inspect.isclass(chromo):
        lengt = 0
        for i in range(0, len(startchromo) - 5,3):  # Comment_DB: fixed to correspond to new ordered starting chromosome, however, this if statement is not used
            lengt = lengt + chromo.genes[i] * l_factor
    else:
        lengt = 0
        for i in range(0, len(startchromo) - 5,3):  # Comment_DB: fixed to correspond to new ordered starting chromosome
            lengt = lengt + chromo[i] * l_factor
    return lengt


# Berechnung der Fitness eines Chromosoms
def Fitness(chromo):

    LoP = ListOfPoints(chromo)
    A = trimesh.proximity.closest_point(testpatch, LoP[0])  # Comment_DB: LoP[0] refers to "result" returned from LoP
    # print("for reordered chromo",A[1]) #Comment_DB: test
    # print("for reordered chromo",LoP[0]) #Comment_DB: test
    L = 0
    for i in range(0, len(startchromo) - 5, 3):  # Comment_DB: fixed to correspond to new orderd starting chromosome
        L = L + chromo[i] * l_factor  # Comment_DB: patch length

    A_norm = []

    for i in range(len(A[1])):  # Comment_DB: A[1] is returned variable "result_distance" from proximity.py module
        A_norm.append(math.sqrt((A[1][i]) ** 2))

    L_aim = start_parameter[4]  # Comment_DB: determined from Sav-Gol curve
    L_aim = L_aim + 45  # todo  Versuch mit L_aim korrigiert

    avg_dist = sum(A_norm) / len(A[1])  # Comment_DB: always positive

    ###PARABOLIC###
    k_d = (100 - 90) / (0.005 ** 2)  # Comment_DB: k_d = 400000
    distance_fit = 100 - k_d * (avg_dist / L_aim) ** 2  # Comment_DB: max distance fitness is 100 (avg_dist = 0)

    ###LINEAR###
    # k_d_lin = (100 - 90) / 0.005
    # distance_fit = 100 - abs(k_d_lin * (avg_dist/L_aim))

    ###GAUSSIAN###
    # k_d_gauss = -math.log(9/10)/(0.005**2) #Comment_DB: deviation of 0.005L_aim --> 90
    # distance_fit = 100 * math.exp(-k_d_gauss*(avg_dist / L_aim) ** 2)

    ###PARABOLIC###
    k_l = (100 - 50) / ((L_aim * 0.2) ** 2)  # Comment_DB: = 1/128 for L_aim = 400. Higher L_aim yields lower k_l
    length_fit = 100 - ((L - L_aim) ** 2) * k_l

    ###LINEAR###
    # k_l_lin = (100-50)/(L_aim*0.2)
    # length_fit = 100 - abs((L - L_aim) * k_l_lin)

    ###GAUSSIAN###
    # k_l_gauss = -math.log(5/10)/((0.2*L_aim) ** 2) #Comment_DB: deviation of 0.2*L_aim --> 50
    # length_fit = 100 * math.exp(-k_l_gauss * (L - L_aim) ** 2)

    ###PARABOLIC###
    k_p = (100 - 90) / (5 ** 2)  # Comment_DB: k_p = 0.4
    border_fit_start = 100 - (stlprep3_6.distance(LoP[1], patch_start) ** 2) * k_p
    # border_fit_start = 100 - 10 * abs(stlprep3_6.distance(LoP[1], patch_start))

    # Comment_DB: Removed k_p. LoP[1] is real patch, patch_start from preproc
    # border_fit_start = 100 * math.exp(-(stlprep3_6.distance(LoP[1], patch_start)/10) ** 2)
    # border_fit_start = abs(600/stlprep3_6.distance(LoP[1],patch_start))

    border_fit_end = 100 - (stlprep3_6.distance(LoP[2], patch_end) ** 2) * k_p  # Comment_DB: trial and error for k_p
    # border_fit_end = abs(600 / stlprep3_6.distance(LoP[2], patch_end))
    # border_fit_end = 100 - 10 * abs(stlprep3_6.distance(LoP[2], patch_end))
    # border_fit_end = 100 * math.exp(-(stlprep3_6.distance(LoP[2], patch_end)/10) ** 2)  #Comment_DB: k_p = 0.4

    ###LINEAR###
    # k_p_lin = (100-90)/5
    # border_fit_start = 100 - abs((stlprep3_6.distance(LoP[1], patch_start)) * k_p_lin)
    # border_fit_end = 100 - abs((stlprep3_6.distance(LoP[2], patch_end)) * k_p_lin)

    ###GAUSSIAN###
    # k_p_gauss = -math.log(9 / 10) / (5 ** 2) #Comment_DB: deviation of 5 mm--> 90
    # border_fit_start = 100 * math.exp(-k_p_gauss * (stlprep3_6.distance(LoP[1], patch_start)) ** 2)
    # border_fit_end = 100 * math.exp(-k_p_gauss * (stlprep3_6.distance(LoP[2], patch_end)) ** 2)

    #####Comment_DB: different gammas!!#####
    global gamma_d
    global gamma_l
    global gamma_ps
    global gamma_pe
    global gamma_d2
    global gamma_l2
    global gamma_ps2
    global gamma_pe2
    global gamma_d3
    global gamma_l3
    global gamma_ps3
    global gamma_pe3
    global gamma_d4
    global gamma_l4
    global gamma_ps4
    global gamma_pe4
    global gamma_d_hat
    global gamma_l_hat
    global gamma_ps_hat
    global gamma_pe_hat

    if not 'p' in globals():
        pass
    else:
        if num_gen_set2 <= p.generationNumber + 1 < num_gen_set3:
            gamma_d_hat = gamma_d2
            gamma_l_hat = gamma_l2
            gamma_ps_hat = gamma_ps2
            gamma_pe_hat = gamma_pe2
        elif num_gen_set3 <= p.generationNumber + 1 < num_gen_set4:
            gamma_d_hat = gamma_d3
            gamma_l_hat = gamma_l3
            gamma_ps_hat = gamma_ps3
            gamma_pe_hat = gamma_pe3
        elif num_gen_set4 <= p.generationNumber + 1 <= num_gen:
            gamma_d_hat = gamma_d4
            gamma_l_hat = gamma_l4
            gamma_ps_hat = gamma_ps4
            gamma_pe_hat = gamma_pe4
        else:  # Comment_DB: If inputs don't make sense, revert back to original gammas for all sets
            gamma_d_hat = gamma_d
            gamma_l_hat = gamma_l
            gamma_ps_hat = gamma_ps
            gamma_pe_hat = gamma_pe

    fitness = distance_fit * gamma_d_hat + length_fit * gamma_l_hat + border_fit_end * gamma_pe_hat + border_fit_start * gamma_ps_hat

    # print("gamma_d", gamma_d)
    # print("fitness",fitness, "border_fit_end", border_fit_end)
    return fitness, distance_fit, length_fit, border_fit_start, border_fit_end, avg_dist


# Erstellung Chromosom der Startlösung
startchromo = preprocessed_chromo(start_parameter, chromo_resolution)

"""
#Ausgabe der Ergebnisse aus Vorverarbeitung der Geometriedaten

print("Anzahl Biegestellen:", AnzahlKnicke)
print("L_aim",L_aim)
print("Start Chromo",startchromo)
print("\tStart Length",PatchLength(startchromo,AnzahlKnicke,l_factor))
print("\tStart Fitness", Fitness(startchromo)[0])
print("\t\tStart Distance Fit",Fitness(startchromo)[1])
print("\t\tStart Length Fit",Fitness(startchromo)[2])
print("\t\tStart Border Fit Start",Fitness(startchromo)[3])
print("\t\tStart Border Fit End",Fitness(startchromo)[4])
print("\t\tStart Average Distance", Fitness(startchromo)[5])
"""

# Ergebnisse aus Vorverarbeitung visualisieren:
stlprep3_6.show_startstrip(input_file, ListOfPoints(startchromo)[0], poly_order, window_quotient, max_distance,
                           ListOfPoints(startchromo)[3], patch_start, patch_end)

####################Evolutionärer Algorithmus####################
## Erzeuge ein Objekt der Klasse Population:
p = Population(pop_size)  # Comment_DB: pop_size is user input in dialog box
## Startwerte aus Preprocessing werden an die Population gegeben
p.startchromo = startchromo
# Soll die Initialisierung mit den Startwerten des Preprozesses erfolgen?
p.preprocessedInit = init_preprocess  # Comment_DB: init_preprocess also user input, preprocessedInit is from chromosome class in Galileo module
p.initrange = 10  # Comment_DB: in chromosome class in Galileo module. Interval for variation of preprocessed gene
p.p_randInit = 8
p.p_prepInit = 70
# Festlegen der Fitnessfunktion
p.evalFunc = Fitness  # Comment_DB: The FUNCTION Fitness is assigned, not the lowercase equation fitness! Stores the function. p.evalFunc = Fitness() stores the return value
## Festlegen der Minimal - & Maximalwerte und der Länge eines Chromosoms in Abh. der Knickanzahl
p.chromoMinValues = [0] * (3 * AnzahlKnicke + 8)
p.chromoMaxValues = [chromo_resolution] * (3 * AnzahlKnicke + 8)

## Ganzzahlige Allele -> useInteger = 1, Floatwerte -> useInteger = 0
p.useInteger = useInteger
# p.useInteger = 0

## Wahl der Selektionsmethode -> Rouletterad, Rankedselect, Eliteranked
# p.selectFunc = p.select_Roulette
# p.selectFunc = p.select_Ranked
# p.selectFunc = p.select_EliteRanked #Comment_DB: function returns elites[k-1]
p.selectFunc = p.select_Roulette  # MW: Convergenz problem mit Elite?

## Wie viele Chromosome dürfen überleben? #Comment_DKu_Wenzel: Das ist die Anzahl an Kinder!!!
# p.replacementSize = p.numChromosomes
p.replacementSize = p.numChromosomes * 5
## Crossover Wahrscheinlichkeit
p.crossoverRate = p_crossover
## Wahl der Crossoverfunktion -> Flat
p.crossoverFunc = p.crossover_Flat
# p.crossoverFunc = p.crossover_Uniform
# p.crossoverFunc = p.crossover_Uniform
## Mutations Wahrscheinlichkeit
p.mutationRate = p_mutation

## Wahl der Mutationsfunktion
if init_preprocess == 0:
    # p.mutateFunc = p.mutate_Default
    p.mutateFunc = p.mutate_Uniform
    p.mutationRange = p_mutate_range
else:
    p.mutateFunc = p.mutate_Uniform
    # p.mutateFunc = p.mutate_Gauss
    p.mutationRange = p_mutate_range

## Replacementfunktion: SteadyState, SteadyState ohne doppelten Chromosome & Generationell(nur Kinder überleben)
# p.replaceFunc = p.replace_SteadyState
p.replaceFunc = p.replace_SteadyStateNoDuplicates
# p.replaceFunc = p.replace_Generational
p.maxGenerations = num_gen


def show_current_solution(ListOfPoints):  # Comment_DB: PATCH OF STARTING CHROMOSOME DETERMINED FROM DATA FROM PREPROCESSOR. SHOWING THE "ListOfPoints"!
    bestPatch = ListOfPoints[0]
    bestPatch_patternpoints = ListOfPoints[3]
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)
    your_mesh = mesh.Mesh.from_file(input_file)  # Comment_DB: This is the target shape!
    patch_visual = mplot3d.art3d.Poly3DCollection(your_mesh.vectors, linewidths=3, alpha=0.5)
    axes.scatter(bestPatch[:, 0], bestPatch[:, 1], bestPatch[:, 2], c='y')

    # Plotten des Patches. Die Knickkantenpunkte werden mit Dreiecken geplottet.
    patch_meshpoints = []
    verts = [list(zip(bestPatch_patternpoints[:, 0], bestPatch_patternpoints[:, 1], bestPatch_patternpoints[:, 2]))]
    for i in range(len(bestPatch_patternpoints) - 2):
        verts = [list(
            zip([bestPatch_patternpoints[i][0], bestPatch_patternpoints[i + 1][0], bestPatch_patternpoints[i + 2][0]],
                [bestPatch_patternpoints[i][1], bestPatch_patternpoints[i + 1][1], bestPatch_patternpoints[i + 2][1]],
                [bestPatch_patternpoints[i][2], bestPatch_patternpoints[i + 1][2], bestPatch_patternpoints[i + 2][2]]))]
        axes.add_collection3d(Poly3DCollection(verts), zs='z')
        patch_meshpoints.append(verts)
    # Biegestellen rot färben:
    axes.scatter(bestPatch_patternpoints[:, 0], bestPatch_patternpoints[:, 1], bestPatch_patternpoints[:, 2], c='r')
    axes.scatter(patch_start[0], patch_start[1], patch_start[2], c='black')
    axes.scatter(patch_end[0], patch_end[1], patch_end[2], c='black')
    face_color = [0.5, 0.5, 1]  # alternative: matplotlib.colors.rgb2hex([0.5, 0.5, 1])
    patch_visual.set_facecolor(face_color)
    axes.add_collection3d(patch_visual)
    axes.autoscale(enable=False, axis='both')  # you will need this line to change the Z-axis
    axes.set_xbound(-100, 100)
    axes.set_ybound(-50, 150)
    axes.set_zbound(-100, 100)
    pyplot.axis('off')
    pyplot.show(figure)
    return


##Comment_DB: initialize arrays of the fitness values (Saving values in the arrays)
num_gen_list = np.array([])
fitness_list = np.array([])
distance_fit_list = np.array([])
length_fit_list = np.array([])
border_fit_start_list = np.array([])
border_fit_end_list = np.array([])
mutation_rate_list = np.array([])
## Initialisierung
time_start = timer()  # Comment_DB: start timer
p.prepPopulation()
p.currentGeneration.sort()  # Comment_DB: sort and reverse randomly initialized pop from best fit to worst
p.currentGeneration.reverse()

if p.preprocessedInit == False:  # Comment_DB: print if random or preprocessor init
    print("##########Random Gene Initialization for Generation 0##########")
else:
    print("##########Preprocessor Gene Initialization for Generation 0##########")

print("##########Generations Sorted from Highest to Lowest Fitness##########")

if adap_mutation == 1:
    print("##########Adaptive Mutation Selected##########")
# Comment_DB: EA Loop (Only this for loop determines the number of generations to iterate! This is not determined in galileo!)
for i in range(num_gen):
    if adap_mutation == 0:
        print("\n#####Elite Population Members of Generation", i, "\b#####")
    else:
        if not i == 0:
            print("\n#####Elite Population Members of Generation", i, "\b.", "Mutation Rate:", p.mutationRate,
                  "\b#####")
        else:
            print("\n#####Elite Population Members of Generation", i, "\b#####")

    ###Comment_DB: (variable weighting) sort the chromosomes again at set gens, since the new gammas change the fitnesses###
    if p.generationNumber == num_gen_set2 - 1 or p.generationNumber == num_gen_set3 - 1 or p.generationNumber == num_gen_set4 - 1:
        p.currentGeneration.sort(key=p.evaluate())
        p.currentGeneration.reverse()

    # Comment_DB: print elite population members and their overall fitness, distance fit, and average distance
    for j in range(p.selectionSize):
        print("\n\tPopulation Member ", j, " :",
              p.currentGeneration[j].genes,
              "\n\t\tMember Fitness:",
              Fitness(p.currentGeneration[j].genes)[0],
              "\tMember Distance Fit:",
              Fitness(p.currentGeneration[j].genes)[1],
              "\tMember Average Distance:",
              Fitness(p.currentGeneration[j].genes)[5]
              )

    if adap_mutation == 1:
        if i == 0:
            mutation_rate_list = np.append(mutation_rate_list, [0])
        else:
            mutation_rate_list = np.append(mutation_rate_list, [p.mutationRate])

    # Comment_DB: Begin Evolutionary Algorithm
    if i != num_gen - 1:
        # Bewertung
        p.evaluate()

        if adap_mutation == 1 and p.generationNumber > num_gen / 2:  # Comment_DB: Increase/decrease mutation rate by 0.5 for adaptive mutation for mutationrate less than 0.5. 0.999 for > 0.5
            if p.avgFitness >= p.bestFitIndividual.fitness - 50 and p.mutationRate == p_mutation:
                if p_mutation >= 0.5:
                    p.mutationRate = 0.999
                else:
                    p.mutationRate = p.mutationRate + 0.5
            elif p.avgFitness >= p.bestFitIndividual.fitness - 50 and (
                    p.mutationRate == p_mutation + 0.5 or p.mutationRate == 0.999):
                pass
            else:
                p.mutationRate = p_mutation

        # Paarungsselektion
        p.select()
        # Rekombination
        p.crossover()
        # Mutation
        p.mutate()
        # Umweltselektion
        p.replace()

    p.generationNumber = p.generationNumber + 1  # Comment_DB: This is one plus the output (Gen 0 has generationNumber 1)

    # print the best fit individual, and its fitness
    print("\nBest Fit Member of Generation ", i, " :", p.bestFitIndividual, "\n\tFitness:", p.bestFitIndividual.fitness,
          "\n\t\tDistance Fit:", Fitness(p.bestFitIndividual.genes)[1],
          "\n\t\tLength Fit:", Fitness(p.bestFitIndividual.genes)[2], "\n\t\tBorder Fit Start:",
          Fitness(p.bestFitIndividual.genes)[3],
          "\n\t\tBorder Fit End:", Fitness(p.bestFitIndividual.genes)[4])
    print("\t\tAverage Distance", Fitness(p.bestFitIndividual.genes)[5])

    print("\n")

    ###Comment_DB: Show iterations for Total Gen 100 or Total Gen 50###
    if num_gen == 100:
        if p.generationNumber == 1 or p.generationNumber == 5 or p.generationNumber == 50 or p.generationNumber == 90:
            show_current_solution(ListOfPoints(p.bestFitIndividual.genes))
    if num_gen == 50:
        if p.generationNumber == 1 or p.generationNumber == 5 or p.generationNumber == 25 or p.generationNumber == 40:
            show_current_solution(ListOfPoints(p.bestFitIndividual.genes))

    # Comment_DB: append determined values into arrays after each iteration
    num_gen_list = np.append(num_gen_list, [i])

    fitness_list = np.append(fitness_list, [p.bestFitIndividual.fitness])
    distance_fit_list = np.append(distance_fit_list, [Fitness(p.bestFitIndividual.genes)[1]])
    length_fit_list = np.append(length_fit_list, [Fitness(p.bestFitIndividual.genes)[2]])
    border_fit_start_list = np.append(border_fit_start_list, [Fitness(p.bestFitIndividual.genes)[3]])
    border_fit_end_list = np.append(border_fit_end_list, [Fitness(p.bestFitIndividual.genes)[4]])

time_end = timer()  # Comment_DB: end timer
fitness_list_gen_index = np.stack((num_gen_list, fitness_list))  # Comment_DB: stack gen_list with list of fitnesses
distance_fit_list_gen_index = np.stack((num_gen_list, distance_fit_list))
length_fit_list_gen_index = np.stack((num_gen_list, length_fit_list))
border_fit_start_list_gen_index = np.stack((num_gen_list, border_fit_start_list))
border_fit_end_list_gen_index = np.stack((num_gen_list, border_fit_end_list))
if adap_mutation == 1:
    mutation_rate_list_gen_index = np.stack((num_gen_list, mutation_rate_list))

print("\n\nEnd Patch length: ", PatchLength(p.bestFitIndividual.genes, AnzahlKnicke, l_factor),
      "L_Aim (From Preprocessor):", L_aim)
print("End Fitness: ", p.bestFitIndividual.getFitness(),
      "\n\tEnd Distance Fit:", Fitness(p.bestFitIndividual.genes)[1],
      "\n\tEnd Length Fit:", Fitness(p.bestFitIndividual.genes)[2],
      "\n\tEnd Border Fit Start:", Fitness(p.bestFitIndividual.genes)[3],
      "\n\tEnd Border Fit End:", Fitness(p.bestFitIndividual.genes)[4])

if adap_mutation == 1:
    print("\tEnd Mutation Rate: ", p.mutationRate)

print("\nSettings Used: ")
print("Set 1 gamma_d: ", gamma_d, "\tSet 2 gamma_d: ", gamma_d2, "\tSet 3 gamma_d: ", gamma_d3, "\tSet 4 gamma_d: ",
      gamma_d4,
      "\nSet 1 gamma_l: ", gamma_l, "\tSet 2 gamma_l: ", gamma_l2, "\tSet 3 gamma_l: ", gamma_l3, "\tSet 4 gamma_l: ",
      gamma_l4,
      "\nSet 1 gamma_ps: ", gamma_ps, "\tSet 2 gamma_ps: ", gamma_ps2, "\tSet 3 gamma_ps: ", gamma_ps3,
      "\tSet 4 gamma_ps: ", gamma_ps4,
      "\nSet 1 gamma_pe: ", gamma_pe, "\tSet 2 gamma_pe: ", gamma_pe2, "\tSet 3 gamma_pe: ", gamma_pe3,
      "\tSet 4 gamma_pe: ", gamma_pe4,
      "\n\nSet 2 Gen Start Point: ", num_gen_set2, "\tSet 3 Gen Start Point: ", num_gen_set3,
      "\tSet 4 Gen Start Point: ", num_gen_set4)

print("\nPop Size: ", pop_size,
      "\n# of Gens: ", num_gen,
      "\nMutation Rate: ", p_mutation, "\tMutation Rate (Final):", p.mutationRate if not adap_mutation == 0 else "N/A",
      "\tMutation Range: ", p_mutate_range,
      "\nCrossover Rate: ", p_crossover if not p.crossoverFunc == p.crossover_Flat else "N/A")

print("\n\nElapsed Time [s]: ", time_end - time_start)

##Comment_DB: Create plots for fitness and subfitnesses
plt.figure(1)
plt.plot(fitness_list_gen_index[0], fitness_list_gen_index[1], linewidth=2)
plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.title('Overall Fitness 1')

plt.figure(2)
plt.plot(distance_fit_list_gen_index[0], distance_fit_list_gen_index[1])
plt.xlabel('Generation')
plt.ylabel('Distance Fitness')
plt.title('Distance Fitness 1')

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

bestPatch = ListOfPoints(p.bestFitIndividual.genes)[0]
bestPatch_patternpoints = ListOfPoints(p.bestFitIndividual.genes)[3]

#############PLOTTING########### (#Comment_DB: Final patch)

figure = pyplot.figure()  # Comment_DB: create a new figure
axes = mplot3d.Axes3D(figure)
your_mesh = mesh.Mesh.from_file(input_file)
patch_visual = mplot3d.art3d.Poly3DCollection(your_mesh.vectors, linewidths=1,
                                              alpha=0.5)  # edgecolor = [1, 1, 1] #Comment_DB: added edgecolor to make the edges visible

testpatch_vector = mesh.Mesh.from_file(input_file)  # Comment_DB: stl mesh. Added to show point cloud
triangles = testpatch_vector.vectors
patch_pc = stlprep3_6.patch_pointcloud(triangles)
axes.scatter(bestPatch[:, 0], bestPatch[:, 1], bestPatch[:, 2], c='y')

# Plotten des Patches. Die Knickkantenpunkte werden mit Dreiecken geplottet.
patch_meshpoints = []
verts = [list(zip(bestPatch_patternpoints[:, 0], bestPatch_patternpoints[:, 1], bestPatch_patternpoints[:, 2]))]
for i in range(len(bestPatch_patternpoints) - 2):
    verts = [
        list(zip([bestPatch_patternpoints[i][0], bestPatch_patternpoints[i + 1][0], bestPatch_patternpoints[i + 2][0]], \
                 [bestPatch_patternpoints[i][1], bestPatch_patternpoints[i + 1][1], bestPatch_patternpoints[i + 2][1]], \
                 [bestPatch_patternpoints[i][2], bestPatch_patternpoints[i + 1][2],
                  bestPatch_patternpoints[i + 2][2]]))]
    axes.add_collection3d(Poly3DCollection(verts), zs='z')
    patch_meshpoints.append(verts)
# Biegestellen rot färben:
axes.scatter(bestPatch_patternpoints[:, 0], bestPatch_patternpoints[:, 1], bestPatch_patternpoints[:, 2], c='r')

patch_meshpoints = np.concatenate(np.asarray(patch_meshpoints), axis=0)

axes.scatter(patch_start[0], patch_start[1], patch_start[2], c='black')
axes.scatter(patch_end[0], patch_end[1], patch_end[2], c='black')

face_color = [0.5, 0.5, 1]  # alternative: matplotlib.colors.rgb2hex([0.5, 0.5, 1])
patch_visual.set_facecolor(face_color)
axes.add_collection3d(patch_visual)

# Show the plot to the screen

axes.autoscale(enable=False, axis='both')  # you will need this line to change the Z-axis
axes.set_xbound(-100, 100)
axes.set_ybound(-50, 150)
axes.set_zbound(-100, 100)
pyplot.axis('off')

pyplot.show(figure)

######## Abspeichern der Biegeparameter #########
def save_patch_file():
    name = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    bestPatch_parameter_l = ListOfPoints(p.bestFitIndividual.genes)[4]
    bestPatch_parameter_alpha = ListOfPoints(p.bestFitIndividual.genes)[5]
    bestPatch_parameter_beta = ListOfPoints(p.bestFitIndividual.genes)[6]
    name.write("width=" + str(width) + "\n")
    name.write("length=" + str(sum(bestPatch_parameter_l)) + "\n")
    name.write("type=" + str(tape_type) + "\n")
    l = 0
    for i in range(AnzahlKnicke):
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

# TODO ####Comment_DB: Save End Fitness Values#####

end = Tk()

Label(end, text="Sind Sie mit dem Patch zufrieden?").grid(row=10, column=1, )
Label(end, justify=LEFT, text=" ").grid(row=11, sticky=W)
Button(end, text="Abbrechen", command=sys.exit).grid(row=30, column=0, )
Button(end, text="Patchparameter speichern", command=save_patch_file).grid(row=30, column=2, )
Label(end, justify=LEFT, text=" ").grid(row=11, sticky=W)
mainloop()