from Tape_EA_Wenzel import *


generation = []

file2 = open("population.txt", "r")
for j in range(pop_size):
    chromo_read = file2.readline()[1:-2]
    chromo_read = chromo_read.split(",")
    chromo_read = list(map(int, chromo_read))
    generation.append(chromo_read)

file2.close()

show_chromo(generation[1])
show_chromo(generation[2])
show_chromo(generation[3])
show_chromo(generation[4])