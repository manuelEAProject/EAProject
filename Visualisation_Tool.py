from Tape_EA_Wenzel import *

def read_in_population(file="population.txt"):
    generation = []
    file2 = open(file, "r")
    for j in range(pop_size):
        chromo_read = file2.readline()[1:-2]
        chromo_read = chromo_read.split(",")
        chromo_read = list(map(int, chromo_read))
        generation.append(chromo_read)
    file2.close()
    return generation

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

individuals = read_in_population()


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