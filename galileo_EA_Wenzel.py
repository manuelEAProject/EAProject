"""A construction kit for creating canonical genetic algorithms"""

# Galileo - the Python Genetic Algorithm module
# Version 1.0b
# Copyright 2003 Donald Goodman
# Based on ga2, the C++ canonical genetic algorithm lib, also by Donald Goodman
# In turn, based on Goldman's definitive (and indeed, canonical) book, /Genetic
# Algorithms/, an excellent reference.
# The author can be contacted at dgoodman@cs.msstate.edu
# the definitive version of this software can (for now) be found at
# http://www.cs.msstate.edu/~dgoodman
# and almost as certainly at
# http://www.sourceforge.net/projects/galileo
#
# for example usage, see the README file
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
import copy
from functools import total_ordering
from random import Random
from random import gauss

from itertools import repeat
import concurrent.futures

@total_ordering
class Chromosome:
    """The Chromosome class represents a single chromosome in a population.
    A Chromosome contains some number of genes (Python objects), and can
    be treated as a list, with indices and slices and all that good stuff
    """

    def __init__(self):
        """Constructs a new Chromosome instance"""
        self.genes = []
        self.geneMaxValues = []
        self.geneMinValues = []
        self.fitness = None
        self.evalFunc = None #Comment_DB: The Fitness(chromo) function
        self.parent = None
        self.initrange = None
        self.p_randInit = None
        self.p_prepInit = None
        self.current_generation_nr = 0 # For evaluating the fitness function, we have to know the current generation

    def __str__(self):
        return self.genes.__str__()

    def randomInit(self, generator, intValues=1):
        """Randomly initializes all genes within the ranges set with
        setMinValues and setMaxValues. generator should be an instance
        of Random from the random module. Doing things this way allows for
        thread safety. If intValues is set to 1, random
        integers will be created, else random floating point values will
        be generated.
        """
        # first, are the lists empty?
        minlen = len(self.geneMinValues)
        maxlen = len(self.geneMaxValues)
        if (minlen == 0) or (minlen != maxlen):
            return

        if intValues == 1:
            randFunc = generator.randint
        else:
            randFunc = generator.uniform
        self.genes = []
        for i in range(minlen):
            self.genes.append(randFunc(self.geneMinValues[i], self.geneMaxValues[i]))


        self.fitness = None

    def preprocessedInit(self, generator, startchromo, intValues=1):
        """Randomly initializes all genes within the ranges set with
        setMinValues and setMaxValues. generator should be an instance
        of Random from the random module. Doing things this way allows for
        thread safety. If intValues is set to 1, random
        integers will be created, else random floating point values will
        be generated. Appends startchromo as well (Comment_DB).
        """


        minlen = len(startchromo)

        if intValues == 1:
            randFunc = generator.randint
        else:
            randFunc = generator.uniform
        self.genes = []
        for i in range(minlen):
            randIntprob = generator.randint(0,100)
            #randVar = generator.random()*self.initrange-self.initrange
            randVar = generator.random() *2* self.initrange - self.initrange # DKu_Wenzel: Now randVar goes from +initrange to -initrange

            if randIntprob <= self.p_randInit:
                self.genes.append(randFunc(self.geneMinValues[i], self.geneMaxValues[i]))
                continue

            elif randIntprob > self.p_prepInit:
                self.genes.append(startchromo[i])

            else:
                if startchromo[i]+randVar/100*self.geneMaxValues[i] > self.geneMaxValues[i]:
                    self.genes.append(self.geneMaxValues[i])
                    continue
                if startchromo[i]+randVar/100*self.geneMaxValues[i]< self.geneMinValues[i]:
                    self.genes.append(self.geneMinValues[i])
                    continue
                else:
                    self.genes.append(int(startchromo[i]+randVar/100*self.geneMaxValues[i]))
                    continue


        self.fitness = None

    def evaluate(self):
        """Calls evalFunc for this chromosome, and caches the fitness value
        returned. Returns None if evalFunc is not yet defined.
        """

        if self.evalFunc != None:
            self.fitness = self.evalFunc(self.genes, self.current_generation_nr)[0] #Comment_DB: Calls Fitness function Fitness(chromo)[0]. Returns first value in Fitness(chromo)

            return self.fitness
        else:
            return None
    def calc_fitness_of_chromo(self):
        self.fitness = None  # Comment_DB: make sure getFitness() returns c.evaluate()
        self.getFitness()
    def getFitness(self):
        """Calls evaluate if there is no cached value, otherwise returns the cached
        fitness value.
        """
        if self.fitness != None: #or (self.fitness == self.evaluate()):
            return self.fitness
        else:
            return self.evaluate()

    def getGenes(self):
        return self.genes

    def copy(self):
        """Duplicates the chromosome.

        retval = Chromosome()
        for item in self.__dict__:
            retval.__dict__[item] = self.__dict__[item]
        return retval
        """
        return copy.deepcopy(self)

    def __len__(self):
        return len(self.genes)

    def __getitem__(self, key):
        retval = self.copy()
        retval.genes = self.genes[key]
        retval.geneMinValues = self.geneMinValues[key]
        retval.geneMaxValues = self.geneMaxValues[key]
        retval.fitness = None
        return retval

    def __setitem__(self, key, value):
        return self.genes.__setitem__(key, value)

    def __getslice__(self, i, j):
        retval = self.copy()
        retval.genes = self.genes[i:j]
        retval.geneMinValues = self.geneMinValues[i:j]
        retval.geneMaxValues = self.geneMaxValues[i:j]
        retval.fitness = None
        return retval
        return self.genes.__getslice__(i, j)

    def __contains__(self, item):
        return self.genes.__contains__(item)

    def __add__(self, other):
        retval = self.copy()
        retval.genes = self.genes + other.genes
        retval.geneMinValues = self.geneMinValues + other.geneMinValues
        retval.geneMaxValues = self.geneMaxValues + other.geneMaxValues
        retval.fitness = None
        return retval

    def __cmp__(self, other):
        s1 = self.getFitness()
        s2 = other.getFitness()
        return s1 - s2

    def __lt__(self, other):
        s1 = self.getFitness()
        s2 = other.getFitness()
        return s1 < s2

    def isIdentical(self, other):
        """If the genes in self and other are identical, returns 0
        """
        return (self.genes == other.genes)
    def isSimilar(self, other):
        """If the genes in self and other are identical, returns 0
        """
        deviation_in_percent = 0.20
        max_deviation = 1 + deviation_in_percent
        min_deviation = 1 - deviation_in_percent
        for i in [-6, -5, -4, -3, -2, -1]:
            if i == -3: #and self.genes[-2] == int(self.geneMaxValues[-2] / 2):
                pass
            else:
                if self.genes[i] < max_deviation * other.genes[i] and self.genes[i] > min_deviation * other.genes[i]:
                    pass
                else:
                    return False

        for i in range(len(self.genes)-6):
            if self.genes[i] < max_deviation* other.genes[i] and self.genes[i] > min_deviation* other.genes[i]:
                pass
            else:
                return False
        return True

"""The Population class represents an entire population of a single
    generation of Chromosomes. This population is replaced with each iteration
    of the algorithm. Functions are provided for storing generations for later
    analysis or retrieval, or for reloading the population from some point.
    All of the high level functionality is in this
    class: generally speaking, you will almost never call a function from any
    of the other classes.
    """
"""Constructs a population of chromosomes, with numChromosomes as the
size of the population. Note that prepPopulation must also be called
after all user defined variables have been set, to finish initialization.
"""

class Population:
    def __init__(self, numChromosomes):

        self.numChromosomes = numChromosomes
        self.generationNumber = 0
        self.maxGenerations = None
        self.currentGeneration = [] #Comment_DB: List of chromosomes in this list
        self.nextGeneration = [] #Comment_DB: for replacement function
        self.chromoMaxValues = []
        self.chromoMinValues = []
        self.startchromo3D = None
        self.startchromo2D = None
        self.startchromo2D_edge = None
        self.selectionSize = int(self.numChromosomes/3)

        self.num_gen_set2 = 0
        self.num_gen_set3 = 0
        self.num_gen_set4 = 0

        self.mutationRate = 0.0
        self.crossoverRate = 0.0
        self.replacementSize = 0
        self.replacementSizeES = 0
        self.useInteger = 0
        self.isSorted = 0

        self.crossoverCount = 0
        self.mutationCount = 0

        self.evalFunc = None
        self.mutateFunc = None
        self.mutationRange = None
        self.selectFunc = None
        self.crossoverFunc = None
        self.replaceFunc = None

        self.preprocessedInit = None
        self.initrange = None
        self.p_randInit = None
        self.p_prepInit = None
        self.generator = Random() #Comment_DB: self.generator is a Random() object

        self.minFitness = None
        self.maxFitness = None
        self.avgFitness = None
        self.sumFitness = None
        self.bestFitIndividual = None

        self.minFitnessNextGen = None
        self.maxFitnessNextGen = None
        self.avgFitnessNextGen = None
        self.sumFitnessNextGen = None
        self.bestFitIndividualNextGen = None
    def prepPopulation(self, use_start_2D, use_start_2D_edge, use_start_3D,use_best_Solution = False): #Comment_DB: For both random and preprocessed initialization
        """Radnomly initializes each chromosome according to the values in
        chromosMinValues and chromosMaxValues.
        """
        if (len(self.chromoMinValues) != len(self.chromoMaxValues)) or (len(self.chromoMinValues) == 0): #Comment_DB: error situation
            return None

        self.currentGeneration = []
        # Add Startchromo to the starting population once (if using preprocessor):
        if self.preprocessedInit == True: #Comment_DB: preprocessed initialization
            amount_of_startchromosoms = 0
            if use_start_3D:
                amount_of_startchromosoms +=1
                c = Chromosome() #Comment_DB: create chromosome object
                c.geneMinValues = self.chromoMinValues #Comment_DB: input in tape_EA
                c.geneMaxValues = self.chromoMaxValues #Comment_DB: input in tape_EA
                for i in range(len(self.startchromo3D)): #Comment_DB: startchromo is the startchromo in tape_EA
                    c.genes.append(self.startchromo3D[i]) #Comment_DB: append the initial startchromo from tape EA!

                c.fitness = None
                c.evalFunc = self.evalFunc
                self.currentGeneration.append(c) #Comment_DB: append the startchromo into the current generation
                c.initrange=self.initrange #Comment_DB: Interval for variation of preprocessed gene. As of now, 1 generation.
            if use_start_2D_edge:
                amount_of_startchromosoms += 1
                c = Chromosome() #Comment_DB: create chromosome object
                c.geneMinValues = self.chromoMinValues #Comment_DB: input in tape_EA
                c.geneMaxValues = self.chromoMaxValues #Comment_DB: input in tape_EA
                for i in range(len(self.startchromo2D_edge)): #Comment_DB: startchromo is the startchromo in tape_EA
                    c.genes.append(self.startchromo2D_edge[i]) #Comment_DB: append the initial startchromo from tape EA!

                c.fitness = None
                c.evalFunc = self.evalFunc
                self.currentGeneration.append(c) #Comment_DB: append the startchromo into the current generation
                c.initrange=self.initrange #Comment_DB: Interval for variation of preprocessed gene. As of now, 1 generation.
            if use_start_2D:
                amount_of_startchromosoms += 1
                c = Chromosome()  # Comment_DB: create chromosome object
                c.geneMinValues = self.chromoMinValues  # Comment_DB: input in tape_EA
                c.geneMaxValues = self.chromoMaxValues  # Comment_DB: input in tape_EA
                for i in range(len(self.startchromo2D)):  # Comment_DB: startchromo is the startchromo in tape_EA
                    c.genes.append(self.startchromo2D[i])  # Comment_DB: append the initial startchromo from tape EA!

                c.fitness = None
                c.evalFunc = self.evalFunc
                self.currentGeneration.append(c)  # Comment_DB: append the startchromo into the current generation
                c.initrange = self.initrange  # Comment_DB: Interval for variation of preprocessed gene. As of now, 1 generation.
            if use_best_Solution:
                amount_of_startchromosoms += 1
                c = Chromosome()  # Comment_DB: create chromosome object
                c.geneMinValues = self.chromoMinValues  # Comment_DB: input in tape_EA
                c.geneMaxValues = self.chromoMaxValues  # Comment_DB: input in tape_EA
                for i in range(len(self.bestFitIndividual)):  # Comment_DB: startchromo is the startchromo in tape_EA
                    c.genes.append(self.bestFitIndividual[i])  # Comment_DB: append the initial startchromo from tape EA!

                c.fitness = None
                c.evalFunc = self.evalFunc
                self.currentGeneration.append(c)  # Comment_DB: append the startchromo into the current generation
                c.initrange = self.initrange  # Comment_DB: Interval for variation of preprocessed gene. As of now, 1 generation.

            partial_pop_size = self.numChromosomes/amount_of_startchromosoms

            pop_size_range = partial_pop_size
            pop_size_range_3D,pop_size_range_2DE,pop_size_range_2D,pop_size_range_best = 0,0,0,0

            if use_start_3D:

                pop_size_range_3D = pop_size_range
                pop_size_range+=partial_pop_size

            if use_start_2D_edge:
                pop_size_range_2DE= pop_size_range
                pop_size_range+=partial_pop_size
            if use_start_2D:
                pop_size_range_2D= pop_size_range
                pop_size_range+=partial_pop_size
            if use_best_Solution:
                pop_size_range_best= pop_size_range
                pop_size_range+=partial_pop_size

        # Create the start population, if PreprocessedInit = True the population is created around Startchromo
        for i in range(self.numChromosomes): #Comment_DB: numChromosomes is the # of chromosomes (population) in each generation
            c = Chromosome() #Comment_DB: create chromosome object
            c.geneMinValues = self.chromoMinValues
            c.geneMaxValues = self.chromoMaxValues
            c.initrange = self.initrange
            c.p_randInit = self.p_randInit
            c.p_prepInit = self.p_prepInit
            # Using the preprocessed parameters (#Comment_DB: to generate the start population?)
            if self.preprocessedInit == True:

                if use_start_3D:
                   if i<pop_size_range_3D:c.preprocessedInit(self.generator, self.startchromo3D, self.useInteger)

                if use_start_2D_edge:
                    if i < pop_size_range_2DE and i >=pop_size_range_3D:c.preprocessedInit(self.generator, self.startchromo2D_edge, self.useInteger)
                if use_start_2D:
                    if i < pop_size_range_2D and i >=pop_size_range_3D and i >=pop_size_range_2DE:c.preprocessedInit(self.generator, self.startchromo2D, self.useInteger)
                if use_best_Solution:
                    if i < pop_size_range_best and i >=pop_size_range_3D and i >=pop_size_range_2DE and i >=pop_size_range_2D:c.preprocessedInit(self.generator, self.bestFitIndividual, self.useInteger)


                    """
                if use_start_2D & use_start_3D & use_start_2D_edge:
                    if i<(self.numChromosomes/3):c.preprocessedInit(self.generator, self.startchromo2D, self.useInteger)
                    elif i < (2*self.numChromosomes / 3):
                        c.preprocessedInit(self.generator, self.startchromo2D_edge, self.useInteger)
                    else:c.preprocessedInit(self.generator, self.startchromo3D, self.useInteger)

                elif use_start_2D & use_start_3D:
                    if i<(self.numChromosomes/2):c.preprocessedInit(self.generator, self.startchromo2D, self.useInteger)
                    else:c.preprocessedInit(self.generator, self.startchromo3D, self.useInteger)
                elif use_start_2D & use_start_2D_edge:
                    if i<(self.numChromosomes/2):c.preprocessedInit(self.generator, self.startchromo2D, self.useInteger)
                    else:c.preprocessedInit(self.generator, self.startchromo2D_edge, self.useInteger)
                elif use_start_3D & use_start_2D_edge:
                    if i<(self.numChromosomes/2):c.preprocessedInit(self.generator, self.startchromo3D, self.useInteger)
                    else:c.preprocessedInit(self.generator, self.startchromo2D_edge, self.useInteger)

                elif use_start_2D: c.preprocessedInit(self.generator, self.startchromo2D, self.useInteger)
                elif use_start_2D_edge: c.preprocessedInit(self.generator, self.startchromo2D_edge, self.useInteger)
                elif use_start_3D: c.preprocessedInit(self.generator, self.startchromo3D, self.useInteger)
                    """
            else:
                c.randomInit(self.generator, self.useInteger)
            c.evalFunc = self.evalFunc
            self.currentGeneration.append(c) #Comment_DB: add chromosomes to the generation (create the population)

        return 1

    def prepPopulation_read_in_population(self, start_2D, start_2D_edge, start_3D, loaded_individums_list_2D=None, loaded_individums_list_2D_E=None, loaded_individums_list_3D=None): #Comment_DB: For both random and preprocessed initialization
        """Radnomly initializes each chromosome according to the values in
        chromosMinValues and chromosMaxValues.
        """
        if (len(self.chromoMinValues) != len(self.chromoMaxValues)) or (len(self.chromoMinValues) == 0): #Comment_DB: error situation
            return None

        self.currentGeneration = []

        # Create the start population
        for i in range(self.numChromosomes): #Comment_DB: numChromosomes is the # of chromosomes (population) in each generation
            c = Chromosome() #Comment_DB: create chromosome object
            c.geneMinValues = self.chromoMinValues
            c.geneMaxValues = self.chromoMaxValues
            c.initrange = self.initrange
            c.p_randInit = self.p_randInit
            c.p_prepInit = self.p_prepInit
            # Using the preprocessed parameters (#Comment_DB: to generate the start population?)
            if self.preprocessedInit == True:

                if start_2D & start_3D & start_2D_edge:
                    if i<(self.numChromosomes/3):
                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_2D[i][j])
                    elif i < (2*self.numChromosomes / 3):
                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_3D[i - int(self.numChromosomes / 3)][j])
                    else:
                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_3D[i - int(2*self.numChromosomes / 3)][j])

                elif start_2D & start_3D:
                    if i < (self.numChromosomes / 2):
                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_2D[i][j])
                    else:

                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_3D[i - int(self.numChromosomes / 2)][j])

                elif start_2D & start_2D_edge:
                    if i < (self.numChromosomes / 2):
                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_2D[i][j])
                    else:

                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_2D_E[i - int(self.numChromosomes / 2)][j])

                elif start_3D & start_2D_edge:
                    if i<(self.numChromosomes/2):
                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_3D[i][j])
                    else:

                        for j in range(len(self.startchromo2D)):
                            c.genes.append(loaded_individums_list_2D_E[i-int(self.numChromosomes/2)][j])

                elif start_2D:
                    for j in range(len(self.startchromo2D)):
                        c.genes.append(loaded_individums_list_2D[i][j])

                elif start_2D_edge:
                    for j in range(len(self.startchromo2D)):
                        c.genes.append(loaded_individums_list_2D_E[i][j])
                elif start_3D:
                    for j in range(len(self.startchromo2D)):
                        c.genes.append(loaded_individums_list_3D[i][j])

            else:
                c.randomInit(self.generator, self.useInteger)

            c.evalFunc = self.evalFunc
            self.currentGeneration.append(c) #Comment_DB: add chromosomes to the generation (create the population)

        return 1
    def evaluate(self):
        """Evaluates each chromosome. Since fitness values are cached, don't
        hesistate to call many times. Also calculates sumFitness, avgFitness,
        maxFitness, minFitness, and finds bestFitIndividual, for your convienence.
        Be sure to assign an evalFunc
        """

        # Usually a generaion is already evaluated. Change of weighting factor at num_gen_set → new evaluation
        if self.generationNumber == self.num_gen_set2 - 1 or self.generationNumber == self.num_gen_set3 - 1 or self.generationNumber == self.num_gen_set4 - 1:
            self.currentGeneration[0].fitness = None  # Comment_DB: make sure getFitness() returns c.evaluate()

        self.sumFitness = 0.0
        self.avgFitness = 0.0

        for chromo in self.currentGeneration:
            chromo.current_generation_nr = self.generationNumber


        self.maxFitness = self.currentGeneration[0].getFitness()
        self.minFitness = self.currentGeneration[0].getFitness()
        self.bestFitIndividual = self.currentGeneration[0]


        for chromo in self.currentGeneration:
            if self.generationNumber == self.num_gen_set2 - 1 or self.generationNumber == self.num_gen_set3 - 1 or self.generationNumber == self.num_gen_set4 - 1:
                chromo.fitness = None  # Comment_DB: make sure getFitness() returns c.evaluate()
            #chromo.fitness = None  # Comment_DB: make sure getFitness() returns c.evaluate()

            ###Comment_DB: (variable weighting) sort the chromosomes again at set gens, since the new gammas change the fitnesses###



            f = chromo.getFitness()

            self.sumFitness = self.sumFitness + f
            if f > self.maxFitness:
                self.maxFitness = f
                self.bestFitIndividual = chromo

            elif f < self.minFitness:
                self.minFitness = f

        self.avgFitness = self.sumFitness / len(self.currentGeneration)  # Comment_DB: can be used if needed


    def calc_fitness_of_chromo(self, chromo):
        chromo.fitness = None  # Comment_DB: make sure getFitness() returns c.evaluate()
        f = chromo.getFitness()
        return f

    def mutate(self):
        """At probability mutationRate, mutates each gene of each chromosome. That
        is, each gene has a mutationRate chance of being randomly re-initialized.
        Right now, only mutate_Default is available for assignment to mutateFunc.
        """

        self.mutationCount = 0

        for i in range(self.replacementSize):
            self.nextGeneration[i] = self.mutateFunc(self.nextGeneration[i])
    def select(self):
        """Selects chromosomes from currentGeneration for placement into
        nextGeneration based on selectFunc.
        """

        self.nextGeneration = []
        for i in range(0, self.replacementSize, 2): #Comment_DB: Each iteration selects chromosome from replacementsize population twice already, thus every 2 indices!
            s1 = self.selectFunc()  # Comment_DB: elite[k-1] chromosome (0th return value in method)
            s2 = self.selectFunc()  # Comment_DB: another elite[k-1] chromosome (0th return value in method)
            s1.parent = (s1, s1)  # Comment_DB: Tuple
            s2.parent = (s2, s2)  # Comment_DB: Tuple

            s1.current_generation_nr = self.generationNumber
            s2.current_generation_nr = self.generationNumber

            self.nextGeneration.append(s1)
            self.nextGeneration.append(s2)  # Comment_DB: append the two elite[k-1] chromosomes into next generation. Do this replacementSize/2 times! (WheelPosition CHANGES! as selectFunc() is called!)


    def select_loop(self):
        s1 = self.selectFunc()  # Comment_DB: elite[k-1] chromosome (0th return value in method)
        s2 = self.selectFunc()  # Comment_DB: another elite[k-1] chromosome (0th return value in method)
        s1.parent = (s1, s1)  # Comment_DB: Tuple
        s2.parent = (s2, s2)  # Comment_DB: Tuple
        self.nextGeneration.append(s1)
        self.nextGeneration.append(s2)  # Comment_DB: append the two elite[k-1] chromosomes into next generation. Do this replacementSize/2 times! (WheelPosition CHANGES! as selectFunc() is called!)

    def crossover(self):
        """Performs crossover on pairs of chromos in nextGeneration with probability
        crossoverRate. Calls crossoverFunc, which must be set; current choices are
        crossover_OnePoint, crossover_TwoPoint and crossover_Uniform.
        """
        self.crossCount = 0

        for i in range(0, self.replacementSize, 2):
            (a, b) = self.crossoverFunc(self.nextGeneration[i], self.nextGeneration[i + 1])
            (self.nextGeneration[i], self.nextGeneration[i + 1]) = (a, b)

    def select_Roulette(self):
        """Perform Roulette (Monte Carlo) selection. Assign this function to
        selectFunc to use.
        In essence, we construct a big roulette wheel, with a slot for each
        individual. The size of each slot is proportional to the relative fitness
        of that individual. The wheel is then spun! whee! The more fit individuals
        have a greater chance of landing under the pointer. The individual that
        lands under the pointer is returned.
        """

        partialSum = 0.0
        # spin the wheel!!
        wheelPosition = self.generator.uniform(0, self.sumFitness)
        i = 0
        for chromo in self.currentGeneration:
            partialSum = partialSum + chromo.getFitness()
            if partialSum >= wheelPosition:
                return chromo
            i = i + 1
        return self.currentGeneration[-1]
    def select_Ranked(self):
        """Perform rank-based selection (using Roulette selection on a ranked list of chromosomes) Assign this function to
        selectFunc to use.
        """
        sortedGeneration = sorted(self.currentGeneration)
        genSize = len(self.currentGeneration)
        # formula for sum of natural numbers
        sumRanks = genSize * (genSize + 1) / 2
        wheelPosition = self.generator.uniform(0, sumRanks)
        pos = 1
        k = 1
        while (pos < wheelPosition):
            k = k + 1
            pos = pos + k
        return sortedGeneration[k - 1]
    def select_EliteRanked(self):
        """Ranked selection only from the top n individuals - the genes of
        the rest will die out. selectionSize is half of the populationSize
        by default. This can be overwritten, but make sure that your numbers
        add up, i.e. don't select 8 individuals from the top 4 parents...
        """
        self.currentGeneration = sorted(self.currentGeneration) #Comment_DB: sorted from worst fitness to best
        elites = self.currentGeneration[-self.selectionSize:] #Comment_DB: the last "selectionsize" chromosomes in the list!
        genSize = len(elites) #Comment_DB: the amount of best chromosomes in the generation
        # formula for sum of natural numbers
        sumRanks = genSize * (genSize + 1) / 2 #Comment_DB: 1+2+3+4+5... = (n*(n+1))/2 --> summing the ranks from 1 to genSize
        wheelPosition = self.generator.uniform(0, sumRanks) #Comment_DB: random number from 0 to sum of ranks (uniform is uniform distribution). Fixed in the while loop
        pos = 1
        k = 1
        while (pos < wheelPosition): #Comment_DB: counter k increases by 1 for each iteration, the higher the k, the higher the pos, thus the more likely the while loop will break
            k = k + 1
            pos = pos + k
        return elites[k - 1] #Comment_DB: Return the (k-1)th chromosome in the elites list!
    def crossover_Flat(self, chromo1, chromo2):
        """A crossover function that can be assigned to crossoverFunc.
        This takes two chromosomes and produces an offspring through linear
        combinations of the parents' genes.
        """
        newchromo1 = chromo1.copy()
        newchromo2 = chromo2.copy()

        for n in range(len(chromo1)):
            r1 = self.generator.uniform(0, 1)
            r2 = self.generator.uniform(0,1)
            if self.useInteger == 1:
                newchromo1[n] = int((r1 * chromo1.genes[n] + (1 - r1) * chromo2.genes[n]))
                newchromo2[n] = int((r2 * chromo1.genes[n] + (1 - r2) * chromo2.genes[n]))
            else:
                newchromo1[n] = int(r1 * chromo1.genes[n] + (1 - r1) * chromo2.genes[n])
                newchromo2[n] = int((r2 * chromo1.genes[n] + (1 - r2) * chromo2.genes[n]))
        newchromo1.fitness = None
        newchromo2.fitness = None
        return (newchromo1,newchromo2)
    def crossover_OnePoint(self, chromo1, chromo2):
        """A crossover function that can be assigned to crossoverFunc. This one
        takes two chromosomes, cuts them at some random point, and swaps the parts
        creating two new chromosomes, which are returned in a tuple. Note
        that there is only a crossoverRate chance of crossover happening.
        """
        prob = self.generator.random()
        if prob <= self.crossoverRate:
            self.crossoverCount = self.crossoverCount + 1
            cutPoint = self.generator.randint(0, len(chromo1) - 1)
            newchromo1 = chromo1[:cutPoint] + chromo2[cutPoint:]
            newchromo2 = chromo2[:cutPoint] + chromo1[cutPoint:]
            newchromo1.fitness = None
            newchromo2.fitness = None

            return (newchromo1, newchromo2)
        else:
            return (chromo1, chromo2)

        """A crossover function that can be assigned to crossoverFunc. This one
        takes two chromosomes, cuts them at two random points (creating three
        parts for each chromosomes), and swaps the parts around, creating two
        new chromosomes, which are returned in a tuple. Note
        that there is only a crossoverRate chance of crossover happening.
        """

        prob = self.generator.random()
        if prob <= self.crossoverRate:
            self.crossoverCount = self.crossoverCount + 1
            cutPoint1 = self.generator.randint(0, len(chromo1) - 1)
            cutPoint2 = self.generator.randint(1, len(chromo1))
            if cutPoint2 < cutPoint1:
                temp = cutPoint1
                cutPoint1 = cutPoint2
                cutPoint2 = temp

            newchromo1 = chromo1[:cutPoint1] + chromo2[cutPoint1:cutPoint2] + chromo1[cutPoint2:]
            newchromo2 = chromo2[:cutPoint1] + chromo1[cutPoint1:cutPoint2] + chromo2[cutPoint2:]

            return (newchromo1, newchromo2)
        else:
            return (chromo1, chromo2)
    def crossover_Uniform(self, chromo1, chromo2):
        """A crossover function that can be assigned to crossoverFunc. Creates
        two new chromosomes by flippinng a coin for each gene. If the coin is heads,
        the gene values in chromo1 and chromo2 are swapped (otherwise they are
        left alone). The two new chromosomes are returned in a tuple. Note
        that there is only a crossoverRate chance of crossover happening.
        """

        prob = self.generator.random()
        if prob <= self.crossoverRate:
            self.crossoverCount = self.crossoverCount + 1
            newchromo1 = chromo1.copy()
            newchromo2 = chromo2.copy()
            for i in range(len(chromo1)):
                # flip a coin...1 we switch, 0 we do nothing
                coin = self.generator.randint(0, 1)
                if coin == 1:
                    temp = newchromo1.genes[i]
                    newchromo1.genes[i] = int(newchromo2.genes[i])
                    newchromo2.genes[i] = int(temp)
            newchromo1.fitness = None
            newchromo2.fitness = None
            return (newchromo1, newchromo2)
        else:
            return (chromo1, chromo2)
    def mutate_Default(self, chromo):
        """Mutation function that can be assigned to mutateFunc. For each gene
        on each chromosome, there is a mutationRate chance that it will be
        randomly re-initialized. The chromosome is returned.
        """
        for i in range(len(chromo.genes)):
            prob = self.generator.random()
            if prob <= self.mutationRate:
                # then we mutate!
                self.mutationCount = self.mutationCount + 1
                f = 0
                if self.useInteger == 1:
                    f = self.generator.randint(self.chromoMinValues[i], self.chromoMaxValues[i])
                else:
                    f = self.generator.uniform(self.chromoMinValues[i], self.chromoMaxValues[i])
                chromo.genes[i] = int(f)
        return chromo
    def mutate_Uniform(self, chromo):
        """Mutation function that can be assigned to mutateFunc.
        Uniform mutation within the mutation radius.
        """
        for gen_nr in range(len(chromo.genes)):
            prob = self.generator.random() #Comment_DB: produces value between 0 (inclusive) and 1 (not inclusive)
            #print("mutation rate in mutate func:",self.mutationRate)
            if prob <= self.mutationRate:
                # then we mutate!
                f = 0
                chromorange = self.chromoMaxValues[gen_nr] - self.chromoMinValues[gen_nr]
                lowstop = max(self.chromoMinValues[gen_nr], chromo.genes[gen_nr] - self.mutationRange * chromorange)
                highstop = min(self.chromoMaxValues[gen_nr], chromo.genes[gen_nr] + self.mutationRange * chromorange)

                if lowstop > self.chromoMaxValues[gen_nr]: lowstop = self.chromoMaxValues[gen_nr] # Comment_DKu_Wenzel: Fehler: lowstop > highstop. Highstop min() ist max chromomax. Lowstop max() kann höher sein. Problem?
                if highstop < 0: highstop = 0

                if self.useInteger == 1:
                    f = self.generator.randint(int(lowstop), int(highstop))
                else:
                    f = self.generator.uniform(
                        self.chromoMinValues[gen_nr], self.chromoMaxValues[gen_nr])

                delta_f = int(f) - chromo.genes[gen_nr]

                chromo.genes[gen_nr] = int(f)

                mutate_compensate = False
                if mutate_compensate:
                    if delta_f > chromorange/36:  # 180°/36 = 5°, ab einer Änderung von 5° wird Kompensiert.
                        nr_of_bendpoints = (len(chromo.genes)-7)/3
                        if gen_nr <= nr_of_bendpoints*3-3:
                            list_alpha = []
                            for j in range(int(nr_of_bendpoints)):
                                list_alpha.append(1+j*3)

                            if (gen_nr in list_alpha):
                                chromo.genes[gen_nr + 3] = chromo.genes[gen_nr + 3] + delta_f
                            else:
                                chromo.genes[gen_nr + 3] = chromo.genes[gen_nr + 3] - delta_f

                            if chromo.genes[gen_nr + 3] < self.chromoMinValues[gen_nr + 3]:
                                chromo.genes[gen_nr + 3] = 1
                            elif chromo.genes[gen_nr + 3] > self.chromoMaxValues[gen_nr + 3]:
                                chromo.genes[gen_nr + 3] = self.chromoMaxValues[gen_nr + 3]
        return chromo
    def mutate_Gauss(self, chromo):
        """Mutation function that can be assigned to mutateFunc.
        Uniform mutation within the mutation radius.
        """
        for i in range(len(chromo.genes)):
            prob = self.generator.random()
            if prob <= self.mutationRate:
                # then we mutate!
                f = 0
                chromorange = self.chromoMaxValues[i] - self.chromoMinValues[i]
                variance = int(chromorange/16-self.generationNumber*(chromorange/16)/self.maxGenerations)
                f = int(gauss(chromo.genes[i],variance))

                if f < self.chromoMinValues[i]:
                    f = self.chromoMinValues[i]
                if f > self.chromoMaxValues[i]:
                    f = self.chromoMaxValues[i]
                chromo.genes[i] = int(f)
        return chromo
    def replace(self):

        """Replaces currentGeneration with nextGeneration according to the rules
        set forth in replaceFunc. Right now, replaceFunc can take the values of
        replace_SteadyState, replace_SteadyStateNoDuplicates and
        replace_Generational.
        """

        return self.replaceFunc()
    def replace_SteadyState(self):
        """Replacement function that can be assigned to replaceFunc.
        Takes the values in nextGeneration, sticks them into currentGeneration, sorts
        currentGeneration, and lops off enough of the least fit individuals
        to reduce the size of currentGeneration back to numChromosomes.
        """

        for chromo in self.nextGeneration:
            self.currentGeneration.append(chromo)
        #self.currentGeneration.sort(key = Chromosome.getFitness(chromo))
        self.currentGeneration.sort() #Comment_DB: sort by fitness
        self.currentGeneration.reverse()
        self.currentGeneration = self.currentGeneration[:self.numChromosomes]
        self.nextGeneration = []
    def replace_SteadyStateNoDuplicates(self):
        """Replacement function that can be assigned to replaceFunc. Same as
        replace_SteadyState, exccept that duplicate chromosomes are not inserted
        back into the currentGeneration.
        """
        # this one is like above, but no duplicates are allowed!
        for chromo in self.nextGeneration:
            flag = 0
            for chromo2 in self.currentGeneration:
                if chromo.isIdentical(chromo2):
                    flag = 1
            if flag == 0:
                self.currentGeneration.append(chromo)

        self.currentGeneration.sort() #Comment_DB: sort by fitness
        self.currentGeneration.reverse()
        self.currentGeneration = self.currentGeneration[:self.numChromosomes]

        self.nextGeneration = []
    def replace_SteadyStateNoSimilares(self):
        """Replacement function that can be assigned to replaceFunc. Same as
        replace_SteadyState, exccept that duplicate chromosomes are not inserted
        back into the currentGeneration.
        """
        # this one is like above, but no duplicates are allowed!
        for chromo in self.nextGeneration:
            flag = 0
            for chromo2 in self.currentGeneration:
                if chromo.isSimilar(chromo2):
                    flag = 1
            if flag == 0:
                self.currentGeneration.append(chromo)

        self.currentGeneration.sort()  # Comment_DB: sort by fitness
        self.currentGeneration.reverse()
        self.currentGeneration = self.currentGeneration[:self.numChromosomes]

        self.nextGeneration = []
    def replace_SteadyStateNoDuplicates_just_offsprings(self):
        """Replacement function that can be assigned to replaceFunc. Same as
        replace_SteadyState, exccept that duplicate chromosomes are not inserted
        back into the currentGeneration.
        """
        # this one is like above, but no duplicates are allowed!
        self.currentGeneration = []
        for chromo in self.nextGeneration:
            flag = 0
            for chromo2 in self.currentGeneration:
                if chromo.isIdentical(chromo2):
                    flag = 1
            if flag == 0:
                self.currentGeneration.append(chromo)

        self.currentGeneration.sort() #Comment_DB: sort by fitness
        self.currentGeneration.reverse()
        self.currentGeneration = self.currentGeneration[:self.numChromosomes]

        self.nextGeneration = []
    def replace_Generational(self):
        """Replacement function that can be assigned to replaceFunc. Wholesale
        replacement of currentGeneration with nextGeneration. assumes that
        replacementSize is equal to numChromosomes; otherwise, the
        currentGeneration will shrink in size to replacementSize in size.
        """
        self.currentGeneration = self.nextGeneration[:]
        self.nextGeneration = []