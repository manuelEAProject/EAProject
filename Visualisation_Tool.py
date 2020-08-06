from Tape_EA_Wenzel import *

# Loaded populations
try: pop_2D = read_in_population(pop_size,"population2_EA.txt")
except: pop_2D =[[]]
try: pop_2DE = read_in_population(pop_size, "population2_edge_EA.txt")
except: pop_2DE =[[]]
try: pop_3D = read_in_population(pop_size, "population3_EA.txt")
except: pop_3D =[[]]

amount_of_bends  = int((len(pop_2D[0])-7)/3)

p_loaded = initialize_Population_with_global_Settings(pop_size, num_gen, amount_of_bends, pop_2D[0], # pop_2D[0] is the best solution of that generation
                                               pop_2DE[0], pop_3D[0])
#p_loaded.prepPopulation_read_in_population(calc_2D_Solution,calc_2D_with_edge_detection, calc_3D_Solution,pop_2D,pop_2DE,pop_3D)
p_loaded.prepPopulation(calc_2D_Solution,calc_2D_with_edge_detection, calc_3D_Solution)
EA_loop(adap_mutation, num_gen, num_gen_set2, num_gen_set3, num_gen_set4, p_loaded, p_mutation)

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