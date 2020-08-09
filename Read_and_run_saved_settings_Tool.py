"""
For now we are reading in just Startchromosoms of saved population setings.
"""

from tkinter import *
from tkinter import filedialog
def GUI_Population():
    global folder_directory
    master = Tk()
    master.protocol("WM_DELETE_WINDOW", sys.exit)
    Label(master, text="Settings for Preprocessor").grid(row=10, sticky=W)
    Label(master, justify=LEFT, text=" ").grid(row=11, sticky=W)
    folder_directory = Entry(master)

    folder_directory.insert(0, "Select folder with settings and population...")
    folder_directory.grid(row=20, columnspan=3, sticky=W + E + N + S)
    Button(text='Browse..', command=lambda: select_folder(folder_directory)).grid(row=20, column=3, sticky=W)
    Button(master, text='Abort', command=sys.exit).grid(row=1000, column=0, pady=4)
    Button(master, text='Start', command=master.quit).grid(row=1000, column=2, pady=4)
    mainloop()  # Executes GUI

    folder_directory = str(folder_directory.get()+"/")

    return folder_directory
def select_folder(folder_directory):

    name_of_folder = filedialog.askdirectory(initialdir="./settings_population/", title="Select file")

    folder_directory.delete(0, 'end')
    folder_directory.insert(0, name_of_folder)
    if len(folder_directory.get()) < 1:
        folder_directory.insert(0, "Select stl-file...")
def create_chromo(lengths, alphas, betas):

    chromo = []

    # Fill length1, alpha1, beta1, length2...
    for i in range(len(lengths)):
        chromo.append(int(lengths[i] / l_factor))
        if i < len(betas):  # Comment_DB: range of beta_list compared to range of l_list is smaller by 1

            #start_chromo.append(int(chromo_resolution / 2))  # Comment_DB: Alphas -> zero on default
            if alphas[i] > math.pi/2:
                chromo.append(int(round( (alphas[i]-((3/4)*math.pi))*(4/(math.pi))* (chromo_resolution / 2))))
            else:
                chromo.append(int(round((chromo_resolution/2)+ (alphas[i]/(math.pi/4)) * (chromo_resolution / 2))))

            beta_chromo = (betas[i]/(math.pi) + 1/2)* chromo_resolution
            chromo.append(int(round(beta_chromo)))

    # Variable Startparameter werden standardmäßig auf chromo_resolution/2 gesetzt
    for i in range(6):
        chromo.append(int(chromo_resolution / 2))
    return chromo

############################## GUI and read in population ####################################################################
sub_dir = GUI_Population()
from Tape_EA_Wenzel import load_settings
load_settings(sub_dir)
from Tape_EA_Wenzel import read_in_last_start_chromosoms, pop_size, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution # import loaded settings
# Initialize Startchromos
startchromo3D = []
startchromo2D = []
startchromo2D_edge = []
# Read in Startchromosoms from sub_dir
read_in_last_start_chromosoms(
            pop_size, use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution,sub_dir)
##################################### Read in all Tape_EA Values ####################################################################################
# Import all variables from TapeEA (read_in_last_start_chromosoms creats startchromos in Tape_EA, have to be exportet afterwards)
from Tape_EA_Wenzel import *



import Tape_EA_Wenzel
############################################# Run EA on importet Population #############################################################################
p_loaded = initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, startchromo2D, startchromo2D_edge, startchromo3D)
                                                                    # pop_2D[0] sind die start_chromosome!!!
#p_loaded.prepPopulation_read_in_population(use_2D_Solution, use_2D_with_edge_detection, use_3D_Solution, pop_2D,pop_2DE,pop_3D)
p_loaded.prepPopulation(use_2D_Solution,use_2D_with_edge_detection, use_3D_Solution)
EA_loop(adap_mutation, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_loaded, p_mutation)


########################################################################################################################################
"""
# Loaded populations
try: pop_2D = read_in_population(pop_size,sub_dir+"population2_EA.txt")
except: pop_2D =[[]]
try: pop_2DE = read_in_population(pop_size, sub_dir+"population2_edge_EA.txt")
except: pop_2DE =[[]]
try: pop_3D = read_in_population(pop_size, sub_dir+"population3_EA.txt")
except: pop_3D =[[]]

# Correction of amount of bends if needed
if amount_of_bends != int((len(pop_2D[0])-7)/3): amount_of_bends = int((len(pop_2D[0])-7)/3)






p_loaded_2D = initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, pop_2D[0], pop_2DE[0], pop_3D[0])
                                                                    # pop_2D[0] sind die start_chromosome!!!
p_loaded_2D.prepPopulation_read_in_population(use_2D_Solution, False, False, pop_2D,pop_2DE,pop_3D)
#p_loaded.prepPopulation(calc_2D_Solution,calc_2D_with_edge_detection, calc_3D_Solution)
save_current_population(p_loaded_2D,"_loaded_2D")
EA_loop(adap_mutation, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_loaded_2D, p_mutation)

"""




"""
individuals = read_in_population(pop_size,"population3_EA.txt")


show_chromo(individuals[0])
show_chromo(individuals[1])
show_chromo(individuals[2])



print(Fitness(individuals[0],1))
print(Fitness(individuals[0], 15))

print(Fitness(individuals[1], 1))
print(Fitness(individuals[1], 15))


[lengths, alphas, betas] = ListOfPoints(individuals[0])[4:7]

lengths[0]=100

chromo = create_chromo(lengths, alphas, betas)

show_chromo(chromo)

print(Fitness(individuals[1], 1))

"""