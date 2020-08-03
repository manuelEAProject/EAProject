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
from Galileo.galileo import Population
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits import mplot3d
from matplotlib import pyplot
import stl_preprocessing
from tkinter import *
from tkinter import filedialog
import pandas as pd
from pandas import ExcelWriter


#TODO: STL Creater aus Chromo

############## GUI  &  SETTINGS #######################
#GUI-Settings


def callback(input_file):
    if os.path.isfile('./settingssheet.txt'):
        settingssheet = open('./settingssheet.txt')
        name = filedialog.askopenfilename(initialdir=settingssheet.readline(),title="Select file",
                                      filetypes=(("stl files", "*.stl"), ("all files", "*.*")))
    else:
        name = filedialog.askopenfilename(initialdir="/", title="Select file",
                                          filetypes=(("stl files", "*.stl"), ("all files", "*.*")))
    input_file.delete(0, 'end')
    input_file.insert(0, name)
    if len(input_file.get()) < 1:
        input_file.insert(0, "Datei auswählen...")

    settingssheet = open('./settingssheet.txt','w+')
    settingssheet.write(input_file.get())
    settingssheet.close

    #print(name)


def save_settings(settings_list):
    settings_sheet = open('./settingssheet.txt', 'r+')
    #for i in range(len(settings_list)):
        #settings_sheet.write(settings_sheet[i]+'/n')
    for listitem in settings_list:
        settings_sheet.write('%s\n' % listitem)
    settings_sheet.close


master = Tk()
master.protocol("WM_DELETE_WINDOW", sys.exit)

Label(master, text="Settings für Tape - Algorithmus").grid(row=10, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=11, sticky=W)
Label(master, text="Verwendung des Preprozess oder Randominitialisierung").grid(row=30, sticky=W)

input_file = Entry(master)
#if os.path.isfile('C:\\Users\\Lukas\\PycharmProjects\\PatchGeometry\\settingssheet.txt'):
if os.path.isfile('./settingssheet.txt'):
    settingssheet = open('./settingssheet.txt')
    input_file.insert(0, settingssheet.readline())
    settingssheet.close()
else:
    input_file.insert(0,"Datei auswählen...")

input_file.grid(row=20,columnspan=3,sticky=W+E+N+S)
Button(text='Browse..', command=lambda : callback(input_file)).grid(row=20,column=3, sticky=W)

init_preprocess = BooleanVar()
init_preprocess.set(True)
init_preprocess_cb = Checkbutton(master, text="Preprozess", variable=init_preprocess,command = lambda :init_random.set(False))
init_preprocess_cb.grid(row=40, column=0, sticky=W)
init_random = BooleanVar()
init_random_cb = Checkbutton(master, text="Random", variable=init_random,command = lambda :init_preprocess.set(False))
init_random_cb.grid(row=40, column=1, sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=41, sticky=W)

Label(master,text="Savitkzy-Golay").grid(row=42,sticky=W)
Label(master,text="Ordnungsgrad").grid(row=43,sticky=W)
poly_order=Entry(master)
poly_order.insert(0,3)
poly_order.grid(row=43,column=1,sticky=W)
Label(master,text="Fensterquotient").grid(row=43,column = 2,sticky = W)
window_quotient=Entry(master)
window_quotient.insert(0,5)
window_quotient.grid(row=43,column=3,sticky=W)
Label(master,text='Maximaler Abstand zur SavGol-Kurve').grid(row=44,sticky=W)
max_distance=Entry(master)
max_distance.insert(0,10)
max_distance.grid(row=44,column=1,sticky=W)



Label(master, text="Tapebreite [mm]:").grid(row=46, sticky=W)
#width = Scale(master, from_=15, to=40,orient=HORIZONTAL)
width=Entry(master)
width.insert(0,20)
width.grid(row=46,column=1,sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=47, sticky=W)

Label(master, text="Punktauflösung Patch").grid(row=50, sticky=W)

fix_number_of_pts = BooleanVar()
fix_number_of_pts_cb = Checkbutton(master, text="Fixe Punktanzahl zwischen Biegestellen", variable=fix_number_of_pts,command=lambda :equidistant_pts_between_bendpts.set(False))
fix_number_of_pts_cb.grid(row=60, column=0, sticky=W)
Label(master,text="Anzahl Punkte:").grid(row=60,column=1,sticky=W)
fix_number_of_pts.set(True)
pointspersection = Entry(master)
pointspersection.insert(0,7)
pointspersection.grid(row=60,column=2,sticky=W)

equidistant_pts_between_bendpts = BooleanVar()
equidistant_pts_between_bendpts_cb = Checkbutton(master, text="Äquidistante Punte", variable=equidistant_pts_between_bendpts,command=lambda :fix_number_of_pts.set(False), disabledforeground="gray")
equidistant_pts_between_bendpts_cb.grid(row=70, column=0, sticky=W)
Label(master,text="Abstand Punkte [mm]:").grid(row=70,column=1,sticky=W)
step_size = Entry(master)
step_size.insert(0,6)
step_size.grid(row=70,column=2,sticky=W)
Label(master, justify=LEFT, text=" ").grid(row=71, sticky=W)

Label(master, text="Einstellungen für den Evolutionären Algorithmus").grid(row=80, sticky=W)

Label(master, justify=LEFT, text="Populationsgröße:").grid(row=90,column=0,sticky=W)
pop_size = Entry(master)
pop_size.insert(0,12)
pop_size.grid(row=90, column=1)

Label(master, text="Anzahl Generationen:").grid(row=100,column=0,sticky=W)
num_gen = Entry(master)
num_gen.insert(0, 80)
num_gen.grid(row=100, column=1)

Label(master,text="Auflösung der Allelparameter:").grid(row=110,column=0,sticky=W)
chromo_resolution = Entry(master)
chromo_resolution.insert(0, 100)
chromo_resolution.grid(row=110, column=1)

Label(master, text="Mutationsrate:").grid(row=120, sticky=W)
p_mutation=Entry(master)
p_mutation.insert(0,0.1)
p_mutation.grid(row=120,column=1,sticky=W)
Label(master, text="Mutationsrange:").grid(row=120,column=2, sticky=W)
p_mutate_range=Entry(master)
p_mutate_range.insert(0,0.05)
p_mutate_range.grid(row=120,column=3,sticky=W)


Label(master, text="Crossoverrate:").grid(row=130, sticky=W)
p_crossover=Entry(master)
p_crossover.insert(0,0.7)
p_crossover.grid(row=130,column=1,sticky=W)

if os.path.isfile('./settingssheet.txt'):

    t = "True"+'\n'

    settingssheet = open('./settingssheet.txt')
    if not settingssheet.readline()=="Datei auswählen...":
        settingssheet.seek(0)
        input_file.delete(0, 'end')
        input = settingssheet.readline()
        input_file.insert(0, input[0:-1])
        if settingssheet.readline()== t:
            init_preprocess.set(True)
        else:
            init_preprocess.set(False)

        if settingssheet.readline()== t:
            init_random.set(True)
        else:
            init_random.set(False)
        #init_random.set(bool(settingssheet.readline()))
        width.delete(0, 'end')
        width.insert(0,float(settingssheet.readline()))
        if settingssheet.readline()== t:
            fix_number_of_pts.set(True)
        else:
            fix_number_of_pts.set(False)

        pointspersection.delete(0, 'end')
        pointspersection.insert(0,int(settingssheet.readline()))
        if settingssheet.readline()== t:
            equidistant_pts_between_bendpts.set(True)
        else:
            equidistant_pts_between_bendpts.set(False)

        step_size.delete(0, 'end')
        step_size.insert(0,float(settingssheet.readline()))
        pop_size.delete(0, 'end')
        pop_size.insert(0,int(settingssheet.readline()))
        num_gen.delete(0, 'end')
        num_gen.insert(0,int(settingssheet.readline()))
        chromo_resolution.delete(0, 'end')
        chromo_resolution.insert(0,int(settingssheet.readline()))
        p_mutation.delete(0, 'end')
        p_mutation.insert(0,float(settingssheet.readline()))
        p_mutate_range.delete(0, 'end')
        p_mutate_range.insert(0,float(settingssheet.readline()))
        p_crossover.delete(0, 'end')
        p_crossover.insert(0,float(settingssheet.readline()))
        poly_order.delete(0,'end')
        poly_order.insert(0,int(settingssheet.readline()))
        window_quotient.delete(0,'end')
        window_quotient.insert(0,int(settingssheet.readline()))
        max_distance.delete(0,'end')
        max_distance.insert(0,int(settingssheet.readline()))
    settingssheet.close()

Button(master, text='Abbrechen', command=sys.exit).grid(row=1000, column=0,pady=4)
Button(master, text='Start', command=master.quit).grid(row=1000, column=1,  pady=4)
mainloop()#führt das GUI aus

settings_list = [input_file.get(),init_preprocess.get(),init_random.get(), width.get(),fix_number_of_pts.get(),\
                 pointspersection.get(),equidistant_pts_between_bendpts.get(),step_size.get(),pop_size.get(),num_gen.get(),chromo_resolution.get(),\
                 p_mutation.get(),p_mutate_range.get(),p_crossover.get(),poly_order.get(),window_quotient.get(),max_distance.get()]




input_file=input_file.get()
print("Inputfile: ",input_file)
#Import stl-file to trimesh and numpy-stl:
testpatch = trimesh.load(input_file)

width = float(width.get())             # Tapebreite
# Sollen die Punkte zwischen den Biegestellen gleichmäßigverteilt werden (True,-> step_size) oder eine fixe Anzahl an Punkten
# zwischen den Biegestellen gewählt werden (False,-> pointspersection)
equidistant_pts_between_bendpts = float(equidistant_pts_between_bendpts.get())
step_size = float(step_size.get())           # Abstand zwischen den Patchpunkten - NUR FALLS equidist_pts_between_bendpts=True

pointspersection = int(pointspersection.get())   # Anzahl der Punkte zwischen den Biegestellen



#### Settings für Evolutionären Algorithmus ####
num_gen=int(num_gen.get())             # Anzahl der Generationen
pop_size=int(pop_size.get())            # Populationsgröße (Anzahl Lösungen pro Generation)
chromo_resolution=int(chromo_resolution.get())    # Auflösung der Allelparameter
# Soll die Initialisierung mit den Startwerten des Preprozesses erfolgen? 1 = True
init_preprocess = bool(init_preprocess.get())
poly_order=int(poly_order.get())
window_quotient=int(window_quotient.get())
max_distance=int(max_distance.get())
# Ganzzahlige Allele -> useInteger = 1, Floatwerte -> useInteger = 0
useInteger = 1
# Crossover Wahrscheinlichkeit
p_crossover = float(p_crossover.get())
# Mutations Wahrscheinlichkeit
p_mutation = float(p_mutation.get())
# Mutationsrange (in welchem Prozentbereich darf das Allel mutieren?)
p_mutate_range = float(p_mutate_range.get())

save_settings(settings_list) # Speichert die gewählten Einstellungen
master.destroy() # Schließt das Settings Fenster

############Vorverarbeitung der Geometriedaten###################

#Ruft das stl_preprocessing modul auf und übergibt die stl-Datei an die Funktion startparam
#Startparam gibt eine Liste mit den berechneten Startparametern zurück:
#STARTPARAMETER: [Start_p_mid, Start_r_mid, Start_n_mid, l_list, totallength, beta_list[°], Start_p_ID]
start_parameter=stl_preprocessing.startparam(input_file,poly_order,window_quotient,max_distance)

Start_p_prep = start_parameter[0]
Start_r_prep = start_parameter[1]
Start_n = start_parameter[2]
Start_p_id = start_parameter[6]
L_aim = start_parameter[4]
patch_start = start_parameter[7]
patch_end = start_parameter[8]
amount_of_bends = len(start_parameter[3]) - 1
# Faktor für das Längenallel in den Chromosomen -> eine Länge kann maximal L_aim lang werden
l_factor = 0.5*L_aim/chromo_resolution

# Chromosom der Startlösung
def preprocessed_chromo(startparameter,chromoresolution):
    # Nimmt die Startparameter aus der Geometriedatenvorverarbeitung und wandelt diese in ein Chromosom mit der
    # entsprechenden Auflösung um. Rückgabewert ist das Chromosom der Startlösung.
    preprocessed_chromo = []
    # Lengths
    for i in range(len(startparameter[3])):
        preprocessed_chromo.append(int(startparameter[3][i]/l_factor))
    # Alphas -> zero on default
    for i in range(len(startparameter[5])):
        preprocessed_chromo.append(int(chromoresolution/2))
    # Betas
    for i in range(len(startparameter[5])):
        beta = start_parameter[5][i]
        beta_chromo = (beta+90)*chromoresolution/180
        preprocessed_chromo.append(int(beta_chromo))
    # Variable Startparameter werden standardmäßig auf chromo_resolution/2 gesetzt
    for i in range(7):
        preprocessed_chromo.append(int(chromo_resolution/2))
    return preprocessed_chromo

# Kinematische Beschreibung des Patchs
def ListOfPoints(chromo):
    l_list = []
    alpha_list = []
    beta_list = []
    for i in range(0, amount_of_bends + 1):
        if chromo[i] == 0:
            l_list.append(1)
        else:
            l_list.append(chromo[i] * l_factor)

    for i in range(amount_of_bends + 1, 2 * amount_of_bends + 1):
        alpha = 0
        if chromo[i] < chromo_resolution / 2:
            alpha = (135 + (chromo[i] * 45 / (chromo_resolution / 2))) * 2 * math.pi / 360
        else:
            # alpha =(90+chromo[i]*4.5)*2*math.pi/360
            # alpha =90+(((chromo[i]-9)**2)*(45/81))
            alpha = ((chromo[i] - chromo_resolution / 2) * 45 / (
                    chromo_resolution / 2)) * 2 * math.pi / 360  # Quadratische Vert. von 135°-180°
        alpha_list.append(alpha)

    for i in range(2 * amount_of_bends + 1, 3 * amount_of_bends + 1):
        beta = (chromo[i] * (180 / chromo_resolution) - 90) * 2 * math.pi / 360
        beta_list.append(beta)

    # Variabler Startpunkt
    var_range = 0.8
    var_start_pt_x = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * amount_of_bends + 1]
    var_start_pt_y = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * amount_of_bends + 2]
    var_start_pt_z = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * amount_of_bends + 3]

    var_start_r_x = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * amount_of_bends + 4]
    var_start_r_y = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * amount_of_bends + 5]
    var_start_r_z = 1 - var_range + (var_range / (chromo_resolution / 2)) * chromo[3 * amount_of_bends + 6]

    gamma_max = 10  # [Grad] Maximaler Kippwinkel für Start_n
    gamma_max_rad = gamma_max * (2 * math.pi / 360)
    var_start_n_gamma = -gamma_max_rad +  gamma_max_rad / (chromo_resolution / 2) * chromo[3 * amount_of_bends + 7]

    Start_p = np.concatenate(np.array(
        [[Start_p_prep[0] * var_start_pt_x], [Start_p_prep[1] * var_start_pt_y], [Start_p_prep[2] * var_start_pt_z]]))

    Start_r = np.concatenate(np.array(
        [[Start_r_prep[0] * var_start_r_x], [Start_r_prep[1] * var_start_r_y], [Start_r_prep[2] * var_start_r_z]]))

    Start_r = 1 / np.linalg.norm(Start_r) * Start_r

    Start_q = np.cross(Start_r, Start_n)
    Start_q = 1 / np.linalg.norm(Start_q) * Start_q

    Start_n_strich = np.cross(Start_q, Start_r)

    Start_n_strich = Quaternion(axis=Start_r, angle=(var_start_n_gamma)).rotate(Start_n_strich)
    Start_n_strich = 1 / np.linalg.norm(Start_n_strich) * Start_n_strich

    # After Startpoint:
    # Start_p_id tells us how many strip-sections (L1-L?) are before the startpoint
    Start_v = np.array([0, 0, 0])
    l_list_a = l_list[Start_p_id:]  # only the length sections after the startpoint
    alpha_list_a = alpha_list[Start_p_id:]

    beta_list_a = beta_list[Start_p_id:]
    r_list_a = [Start_r]
    v_list_a = [Start_v]
    n_list_a = [Start_n_strich]

    #### Vektorenberechnung für Tapeseite nach dem Startpunkt
    for i in range(1, len(l_list_a)):
        if alpha_list_a[i - 1] < math.pi / 2:
            v_new_a = Quaternion(axis=n_list_a[i - 1], angle=(alpha_list_a[i - 1] - (math.pi) / 2)).rotate(
                r_list_a[i - 1])
        else:
            v_new_a = Quaternion(axis=n_list_a[i - 1], angle=(alpha_list_a[i - 1] - 3 * (math.pi) / 2)).rotate(
                r_list_a[i - 1])
        v_list_a.append(v_new_a)
        r_rotated_a = Quaternion(axis=v_list_a[i], angle=beta_list_a[i - 1]).rotate(r_list_a[i - 1])
        r_list_a.append(r_rotated_a)
        n_new_a = Quaternion(axis=v_list_a[i], angle=beta_list_a[i - 1]).rotate(n_list_a[i - 1])
        n_list_a.append(n_new_a)

    r_list_a = np.stack(r_list_a)

    ## Mittellinie after startpoint
    p_list_a = [Start_p]
    old_p_a = Start_p
    # new_p = Start_p
    for i, l in enumerate(l_list_a):
        dist = l * r_list_a[i, :]
        new_p_a = old_p_a + dist
        p_list_a.append(new_p_a)
        old_p_a = new_p_a

    # Auffüllen von Punkten zwischen den Biegestellen auf Mittellinie. Entweder fixe Anzahl an Punkten oder äquidistant
    if equidistant_pts_between_bendpts:
        middle_pts_a = []
        for i, l in enumerate(l_list_a):
            a_new_a = p_list_a[i][np.newaxis, :] + np.outer(
                np.linspace(0, l, math.floor(l / step_size), endpoint=False),
                r_list_a[i])
            middle_pts_a.append(a_new_a)
    else:
        middle_pts_a = []
        for i, l in enumerate(l_list_a):
            a_new_a = p_list_a[i][np.newaxis, :] + np.outer(np.linspace(0, l, pointspersection, endpoint=False),
                                                            r_list_a[i])
            middle_pts_a.append(a_new_a)

    middle_pts_a = np.concatenate(middle_pts_a)

    ########## Linker Rand nach Startpunkt ######################
    ## Anpassung Startpunkt Seitenrand links
    alpha_start = alpha_list[Start_p_id - 1]
    if alpha_start > math.pi / 2:
        delta_l_l_start = (width / 2) * math.tan(math.pi - alpha_start)
    if alpha_start < math.pi / 2:
        delta_l_l_start = - (width / 2) * math.tan(alpha_start)

    # Linker Rand
    delta_l_list = [delta_l_l_start]
    l_left_list_a = []
    for i in range(1, len(l_list_a)):
        if alpha_list_a[i - 1] > math.pi / 2:
            delta_l_l = (width / 2) * math.tan(math.pi - alpha_list_a[i - 1])
            l_left_new = l_list_a[i - 1] + delta_l_l - delta_l_list[i - 1]
            l_left_list_a.append(l_left_new)
            delta_l_list.append(delta_l_l)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_l
            if i == len(l_list_a) - 1:
                l_left_list_a.append((l_list_a[-1] + delta_l_end))

        if alpha_list_a[i - 1] < math.pi / 2:
            delta_l_l = - (width / 2) * math.tan(alpha_list_a[i - 1])
            l_left_new = l_list_a[i - 1] + delta_l_l - delta_l_list[i - 1]
            l_left_list_a.append(l_left_new)
            delta_l_list.append(delta_l_l)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_l
            if i == len(l_list_a) - 1:
                l_left_list_a.append((l_list_a[-1] + delta_l_end))

    # Ausnahme für keinen Knick nach Start:
    if len(l_list_a) == 1:
        l_left_list_a.append((l_list_a[-1]) + delta_l_l_start)

    # Eckpunkte Left
    Start_p_left_a = Start_p - np.cross(Start_r, Start_n_strich) * width / 2 + delta_l_l_start * Start_r
    p_left_list_a = [Start_p_left_a]
    old_p_left_a = Start_p_left_a
    for i, l in enumerate(l_left_list_a):
        dist = l * r_list_a[i, :]
        new_p_left_a = old_p_left_a + dist
        p_left_list_a.append(new_p_left_a)
        old_p_left_a = new_p_left_a
    # Auffüllen mit Punkten ( Fixe Punktanzahl oder fixer Punktabstand):
    if equidistant_pts_between_bendpts:
        left_pts_a = []
        for i, l in enumerate(l_left_list_a):
            if l < 0:
                continue
            else:
                a_new_l = p_left_list_a[i][np.newaxis, :] + np.outer(np.linspace(0, l, math.floor(l / step_size),
                                                                                 endpoint=False), r_list_a[i])
                left_pts_a.append(a_new_l)
    else:
        left_pts_a = []
        for i, l in enumerate(l_left_list_a):
            if l < 0:
                continue
            else:
                a_new_l = p_left_list_a[i][np.newaxis, :] + np.outer(np.linspace(0, l, pointspersection,
                                                                                 endpoint=False), r_list_a[i])
                left_pts_a.append(a_new_l)
    left_pts_a = np.concatenate(left_pts_a)
    ######################## Rechter Rand nach Startpunkt######################
    ## Anpassung Startpunkt Seitenrand rechts
    alpha_start = alpha_list[Start_p_id - 1]
    if alpha_start > math.pi / 2:
        delta_l_r_start = -(width / 2) * math.tan(math.pi - alpha_start)
    if alpha_start < math.pi / 2:
        delta_l_r_start = (width / 2) * math.tan(alpha_start)

    delta_r_list = [delta_l_r_start]
    l_right_list_a = []
    for i in range(1, len(l_list_a)):
        if alpha_list_a[i - 1] > math.pi / 2:
            delta_l_r = - (width / 2) * math.tan(math.pi - alpha_list_a[i - 1])
            l_right_new_a = l_list_a[i - 1] + delta_l_r - delta_r_list[i - 1]
            l_right_list_a.append(l_right_new_a)
            delta_r_list.append(delta_l_r)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_r
            if i == len(l_list_a) - 1:
                l_right_list_a.append((l_list_a[-1] + delta_l_end))

        if alpha_list_a[i - 1] < math.pi / 2:
            delta_l_r = (width / 2) * math.tan(alpha_list_a[i - 1])
            l_right_new_a = l_list_a[i - 1] + delta_l_r - delta_r_list[i - 1]
            l_right_list_a.append(l_right_new_a)
            delta_r_list.append(delta_l_r)
            # falls keine weiteren Knicke:
            delta_l_end = - delta_l_r
            if i == len(l_list_a) - 1:
                l_right_list_a.append((l_list_a[-1] + delta_l_end))

    # Ausnahme für keinen Knick nach Start:
    if len(l_list_a) == 1:
        l_right_list_a.append((l_list_a[-1]) + delta_l_r_start)
    # print("längen", l_left_list_a,'\n',l_list_a,'\n', l_right_list_a)

    # Eckpunkte Rechts
    Start_p_right_a = Start_p + np.cross(Start_r, Start_n_strich) * width / 2 + delta_l_r_start * Start_r
    p_right_list_a = [Start_p_right_a]
    old_p_right_a = Start_p_right_a
    for i, l in enumerate(l_right_list_a):
        dist = l * r_list_a[i, :]
        new_p_right_a = old_p_right_a + dist
        p_right_list_a.append(new_p_right_a)
        old_p_right_a = new_p_right_a
    if equidistant_pts_between_bendpts:
        right_pts_a = []
        for i, l in enumerate(l_right_list_a):
            if l < 0:
                continue
            else:
                a_new_r = p_right_list_a[i][np.newaxis, :] + np.outer(
                    np.linspace(0, l, math.floor(l / step_size), endpoint=False), r_list_a[i])
                right_pts_a.append(a_new_r)
    else:
        right_pts_a = []
        for i, l in enumerate(l_right_list_a):
            if l < 0:
                continue
            else:
                a_new_r = p_right_list_a[i][np.newaxis, :] + np.outer(
                    np.linspace(0, l, pointspersection, endpoint=False), r_list_a[i])
                right_pts_a.append(a_new_r)
    right_pts_a = np.concatenate(right_pts_a)

    #################### Before Startpoint:
    # Start_p_id tells us how many strip-sections (L1-L?) are before the startpoint

    l_list_b = l_list[:Start_p_id]  # only the length sections before the startpoint
    l_list_b = l_list_b[::-1]  # reverse order -> from startpoint to end

    alpha_list_b = alpha_list[:Start_p_id]
    alpha_list_b = alpha_list_b[::-1]
    beta_list_b = beta_list[:Start_p_id]
    beta_list_b = beta_list_b[::-1]

    v_list_b = [Start_v]
    Start_r_b = -1 * Start_r
    r_list_b = [Start_r_b]
    n_list_b = [Start_n_strich]

    #### Vektorenberechnung für Tapeseite vor dem Startpunkt
    if not len(alpha_list_b) == 0:
        for i in range(1, len(l_list_b) + 1):
            # v-Vektoren müssen auf rechte Seite des Tapes zeigen
            if alpha_list_b[i - 1] > math.pi / 2:
                v_new_b = Quaternion(axis=n_list_b[i - 1], angle=(alpha_list_b[i - 1] - (math.pi) / 2)).rotate(
                    r_list_b[i - 1])
            else:
                v_new_b = Quaternion(axis=n_list_b[i - 1], angle=(alpha_list_b[i - 1] + (math.pi) / 2)).rotate(
                    r_list_b[i - 1])
            v_list_b.append(v_new_b)
            # Vorzeichen von Beta umkehren, da in umgedrehte Richtung gegangen wird
            r_rotated_b = Quaternion(axis=v_list_b[i], angle=- beta_list_b[i - 1]).rotate(r_list_b[i - 1])
            r_list_b.append(r_rotated_b)
            n_new_b = Quaternion(axis=v_list_b[i], angle=- beta_list_b[i - 1]).rotate(n_list_b[i - 1])
            n_list_b.append(n_new_b)
    else:
        print("lenalpha0")

    r_list_b = np.asarray(r_list_b)
    r_list_b = np.delete(r_list_b, 0, axis=0)

    v_list_b = np.asarray(v_list_b)
    v_list_b = np.delete(v_list_b, 0, axis=0)
    n_list_b = np.asarray(n_list_b)
    n_list_b = np.delete(n_list_b, 0, axis=0)

    p_list_b = [Start_p]
    old_p_b = Start_p

    for i, l in enumerate(l_list_b):
        dist = l * r_list_b[i, :]
        new_p_b = old_p_b + dist
        p_list_b.append(new_p_b)
        old_p_b = new_p_b

    if equidistant_pts_between_bendpts:
        middle_pts_b = []
        for i, l in enumerate(l_list_b):
            a_new_b = p_list_b[i][np.newaxis, :] + np.outer(
                np.linspace(0, l, math.floor(l / step_size), endpoint=False),
                r_list_b[i])
            middle_pts_b.append(a_new_b)
    else:
        middle_pts_b = []
        for i, l in enumerate(l_list_b):
            a_new_b = p_list_b[i][np.newaxis, :] + np.outer(np.linspace(0, l, pointspersection, endpoint=False),
                                                            r_list_b[i])
            middle_pts_b.append(a_new_b)

    middle_pts_b = np.concatenate(middle_pts_b)

    # alpha im Startpunkt eliminieren (wird durch delta_l_r / delta_l_l berücksichtigt)
    alpha_list_b = np.asarray(alpha_list_b)
    alpha_list_b = np.delete(alpha_list_b, 0, axis=0)

    ################ Linker Rand vor Startpunkt###########
    # Links im Sinne der Hauptrichtung (nach Start)
    delta_l_list = [-delta_l_l_start]
    l_left_list_b = []
    if len(alpha_list_b) == 0:
        l_left_list_b = [l_list_b[0]]
    for i in range(1, len(l_list_b)):
        if alpha_list_b[i - 1] > math.pi / 2:
            delta_l_l = -(width / 2) * math.tan(math.pi - alpha_list_b[i - 1])
            l_left_new_a = l_list_b[i - 1] + delta_l_l - delta_l_list[i - 1]
            l_left_list_b.append(l_left_new_a)
            delta_l_list.append(delta_l_l)
            delta_l_end = - delta_l_l
            if i == len(l_list_b) - 1:
                l_left_list_b.append((l_list_b[-1] + delta_l_end))

        if alpha_list_b[i - 1] < math.pi / 2:
            delta_l_l = (width / 2) * math.tan(alpha_list_b[i - 1])
            l_left_new_a = l_list_b[i - 1] + delta_l_l - delta_l_list[i - 1]
            l_left_list_b.append(l_left_new_a)
            delta_l_list.append(delta_l_l)
            delta_l_end = - delta_l_l
            if i == len(l_list_b) - 1:
                l_left_list_b.append((l_list_b[-1] + delta_l_end))




    # print("l_left_b",l_left_list_b)

    # Eckpunkte Left
    Start_p_left_b = Start_p_left_a
    p_left_list_b = [Start_p_left_b]
    old_p_left_b = Start_p_left_b
    for i, l in enumerate(l_left_list_b):
        dist = l * r_list_b[i, :]
        new_p_left_b = old_p_left_b + dist
        p_left_list_b.append(new_p_left_b)
        old_p_left_b = new_p_left_b

    if equidistant_pts_between_bendpts:
        left_pts_b = []
        for i, l in enumerate(l_left_list_b):
            if l < 0:
                continue
            else:
                a_new_l = p_left_list_b[i][np.newaxis, :] + np.outer(np.linspace(0, l, math.floor(l / step_size),
                                                                                 endpoint=False), r_list_b[i])
                left_pts_b.append(a_new_l)
        if len(left_pts_b) == 0:
            left_pts_b.append(Start_p)
            left_pts_b = np.asarray(left_pts_b)
        else:
            left_pts_b = np.concatenate(left_pts_b)
    else:
        left_pts_b = []
        for i, l in enumerate(l_left_list_b):
            if l < 0:
                continue
            else:
                a_new_l = p_left_list_b[i][np.newaxis, :] + np.outer(np.linspace(0, l, pointspersection,
                                                                                 endpoint=False), r_list_b[i])
                left_pts_b.append(a_new_l)
        if len(left_pts_b) == 0:
            left_pts_b.append(Start_p)
            left_pts_b = np.asarray(left_pts_b)
        else:
            left_pts_b = np.concatenate(left_pts_b)

    delta_r_list = [-delta_l_r_start]
    l_right_list_b = []
    if len(alpha_list_b) == 0:
        l_right_list_b = [l_list_b[0]]
    for i in range(1, len(l_list_b)):
        if alpha_list_b[i - 1] > math.pi / 2:
            delta_l_r = (width / 2) * math.tan(math.pi - alpha_list_b[i - 1])
            l_right_new_b = l_list_b[i - 1] + delta_l_r - delta_r_list[i - 1]
            l_right_list_b.append(l_right_new_b)
            delta_r_list.append(delta_l_r)
            delta_l_end = - delta_l_r
            if i == len(l_list_b) - 1:
                l_right_list_b.append((l_list_b[-1] + delta_l_end))

        if alpha_list_b[i - 1] < math.pi / 2:
            delta_l_r = - (width / 2) * math.tan(alpha_list_b[i - 1])
            l_right_new_b = l_list_b[i - 1] + delta_l_r - delta_r_list[i - 1]
            l_right_list_b.append(l_right_new_b)
            delta_r_list.append(delta_l_r)
            delta_l_end = - delta_l_r
            if i == len(l_list_b) - 1:
                l_right_list_b.append((l_list_b[-1] + delta_l_end))




    # Eckpunkte Rechts
    Start_p_right_b = Start_p_right_a
    p_right_list_b = [Start_p_right_b]
    old_p_right_b = Start_p_right_b
    for i, l in enumerate(l_right_list_b):
        dist = l * r_list_b[i, :]
        new_p_right_b = old_p_right_b + dist
        p_right_list_b.append(new_p_right_b)
        old_p_right_b = new_p_right_b

    if equidistant_pts_between_bendpts:
        right_pts_b = []
        for i, l in enumerate(l_right_list_b):
            if l < 0:
                continue
            else:
                a_new_r = p_right_list_b[i][np.newaxis, :] + np.outer(
                    np.linspace(0, l, math.floor(l / step_size), endpoint=False), r_list_b[i])
                right_pts_b.append(a_new_r)

        if len(right_pts_b) == 0:
            right_pts_b.append(Start_p)
            right_pts_b = np.asarray(right_pts_b)
        else:
            right_pts_b = np.concatenate(right_pts_b)
    else:
        right_pts_b = []
        for i, l in enumerate(l_right_list_b):
            if l < 0:
                continue
            else:
                a_new_r = p_right_list_b[i][np.newaxis, :] + np.outer(
                    np.linspace(0, l, pointspersection, endpoint=False), r_list_b[i])
                right_pts_b.append(a_new_r)

        if len(right_pts_b) == 0:
            right_pts_b.append(Start_p)
            right_pts_b = np.asarray(right_pts_b)
        else:
            right_pts_b = np.concatenate(right_pts_b)

    result = np.concatenate((middle_pts_a, left_pts_a, right_pts_a, middle_pts_b, left_pts_b, right_pts_b), axis=0)
    start = middle_pts_b[-1]
    end = middle_pts_a[-1]

    # Nur Eck- und Biegepunkte des Tapes für einfache Visualisierung
    patch_visualisation_points = []
    p_left_list_b = p_left_list_b[::-1]
    p_right_list_b = p_right_list_b[::-1]
    for i in range(len(p_left_list_b)):
        patch_visualisation_points.append(p_left_list_b[i])
        patch_visualisation_points.append(p_right_list_b[i])
    for i in range(len(p_left_list_a)):
        patch_visualisation_points.append(p_left_list_a[i])
        patch_visualisation_points.append(p_right_list_a[i])
    patch_visualisation_points = np.stack(patch_visualisation_points, axis=0)

    return result, start, end, patch_visualisation_points, l_list, alpha_list, beta_list, Start_p, Start_r, Start_r_b, \
           p_list_b, v_list_b, n_list_b, r_list_b

# Berechnung der Patchlänge (für Fitness)
def PatchLength(chromo,AnzahlKnicke,l_factor):
    # Berechnet die Länge eines Patches. Kann Chomosome als class chromosome oder auch als einfache Liste auslesen.

    if inspect.isclass(chromo):
        L = 0
        for i in range(int(AnzahlKnicke+1)):
            #print(chromo.genes[i])
            L = L + chromo.genes[i]*l_factor
    else:
        L = 0
        for i in range(int(AnzahlKnicke+1)):
            L = L + chromo[i]*l_factor
    return L

# Berechnung der Fitness eines Chromosoms
def Fitness(chromo):
    genes = chromo
    LoP = ListOfPoints(chromo)
    A = trimesh.proximity.closest_point(testpatch, LoP[0])

    L = 0
    for i in range(amount_of_bends + 1):
        L = L + genes[i] * l_factor

    A_norm = []
    A_end = []
    for i in range(len(A[1])):
        A_norm.append(math.sqrt((A[1][i]) ** 2))

    L_aim = start_parameter[4]

    max_dist = max(A_norm)
    avg_dist = sum(A_norm) / len(A[1])
    k_d = 10/2.25
    distance_fit = 100 - k_d * avg_dist ** 2
    if distance_fit < 0:
        distance_fit = 0.1

    k_l = (100-50)/((L_aim*0.2)**2)
    length_fit = 100 - ((L - L_aim) ** 2) *k_l
        #length_fit = 0
    if length_fit<0:
        length_fit=0.1
    k_p = (100-90)/(5**2)
    border_fit_start = 100 - (stl_preprocessing.distance(LoP[1], patch_start)**2)*k_p
    if border_fit_start < 0:
        border_fit_start = 0.1
    border_fit_end = 100 - (stl_preprocessing.distance(LoP[2], patch_end)**2)*k_p
    if border_fit_end < 0:
        border_fit_end = 0.1

    # print("dist",distance_fit)
    # print("maxdist",max_dist)
    # print("avgdist",avg_dist)
    # print("length",length_fit)
    k_l = 1.4
    k_p = 1.4
    k_d = 10 - 3 * k_l
    fitness = distance_fit * k_d + length_fit * k_l + border_fit_end * k_p + border_fit_start * k_p
    return fitness

# Erstellung Chromosom der Startlösung
startchromo = preprocessed_chromo(start_parameter,chromo_resolution)


#Ausgabe der Ergebnisse aus Vorverarbeitung der Geometriedaten

print("Anzahl Biegestellen:", amount_of_bends)
print("L_aim",L_aim)
print("startchromo",startchromo)
print("startlength", PatchLength(startchromo, amount_of_bends, l_factor))
print("startfitness", Fitness(startchromo))
# Ergebnisse aus Vorverarbeitung visualisieren:
#stl_preprocessing.show_startstrip(input_file,ListOfPoints(startchromo)[0],poly_order,window_quotient,max_distance,ListOfPoints(startchromo)[3])




####################Evolutionärer Algorithmus####################
## Erzeuge ein Objekt der Klasse Population:
p = Population(pop_size)
## Startwerte aus Preprocessing werden an die Population gegeben
p.startchromo = startchromo
# Soll die Initialisierung mit den Startwerten des Preprozesses erfolgen?
p.preprocessedInit = init_preprocess
p.initrange = 10
p.p_randInit = 8
p.p_prepInit = 70
# Festlegen der Fitnessfunktion
p.evalFunc = Fitness
## Festlegen der Minimal - & Maximalwerte und der Länge eines Chromosoms in Abh. der Knickanzahl
p.chromoMinValues = [0]*(3 * amount_of_bends + 8)
p.chromoMaxValues = [chromo_resolution]*(3 * amount_of_bends + 8)
## Ganzzahlige Allele -> useInteger = 1, Floatwerte -> useInteger = 0
p.useInteger = useInteger
#p.useInteger = 0

## Wahl der Selektionsmethode -> Rouletterad, Rankedselect, Eliteranked
#p.selectFunc = p.select_Roulette
#p.selectFunc = p.select_Ranked
p.selectFunc = p.select_EliteRanked

## Wie viele Chromosome dürfen überleben?
p.replacementSize = p.numChromosomes
## Crossover Wahrscheinlichkeit
p.crossoverRate = p_crossover
## Wahl der Crossoverfunktion -> Flat
p.crossoverFunc = p.crossover_Flat
#p.crossoverFunc = p.crossover_Uniform
## Mutations Wahrscheinlichkeit
p.mutationRate = p_mutation
## Wahl der Mutationsfunktion
if init_preprocess == 0 :
    p.mutateFunc = p.mutate_Default
    p.mutateFunc = p.mutate_Uniform
    p.mutationRange = p_mutate_range
else:
    p.mutateFunc = p.mutate_Uniform
    #p.mutateFunc = p.mutate_Gauss
    p.mutationRange = p_mutate_range

## Replacementfunktion: SteadyState, SteadyState ohne doppelten Chromosome & Generationell(nur Kinder überleben)
#p.replaceFunc = p.replace_SteadyState
p.replaceFunc = p.replace_SteadyStateNoDuplicates
#p.replaceFunc = p.replace_Generational
p.maxGenerations=num_gen


#df=pd.DataFrame({})

#with ExcelWriter('convergence.xlsx') as writer:
    #df.to_excel(writer)
Fitness_list = []
## Initialisierung
p.prepPopulation()
#EA Schleife
for i in range(num_gen):
    # Bewertung
    p.evaluate()
    # Paarungsselektion
    p.select()
    # Rekombination
    p.crossover()
    # Mutation
    p.mutate()
    # Umweltselektion
    p.replace()

    p.generationNumber=p.generationNumber+1
    # print the best fit individual, and its fitness
    print("Generation ",i," :",p.bestFitIndividual, p.bestFitIndividual.fitness)
    Fitness_list.append(p.bestFitIndividual.fitness)

fit1 = pd.DataFrame({'Fit1': Fitness_list})
with ExcelWriter('convergence.xlsx', mode='a') as writer:
    fit1.to_excel(writer, sheet_name='Convergence')





