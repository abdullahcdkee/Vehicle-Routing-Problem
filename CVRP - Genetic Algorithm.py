###########################################################################################
###########################################################################################
###########################################################################################
################  Code by Abdullah Siddiqui (ms03586@st.habib.edu.pk) #####################
###########################################################################################
###########################################################################################
###########################################################################################

from scipy.spatial.distance import pdist, squareform
import numpy as np
import sys
import networkx as nx
import random
import matplotlib.pyplot as plt
import math
import numpy as np
import random
import matplotlib.pyplot as plt
import math


# Class that represents the genetic algorithm and associated functions
class EA:

	def __init__(self, nodelist, nodecoords, nodecap, size, journeys, max_cap, mutation_rate = 0.5):
		self.nodelist = nodelist
		self.nodecoords = nodecoords
		self.nodecap = nodecap
		self.population_size = size
		self.mutation_rate = mutation_rate
		self.population = [] 
		self.journeys = journeys
		self.max_cap = max_cap
		self.mutation_valid = 0
		self.mutation_invalid = 0
		self.initial_solution = []

		self.initialize_population()


	# Initialize population randomly while following capacity requirements
	def initialize_population(self):

		candidate = nodelist.copy()

		for i in range(self.population_size):
			route = []

			while (len(route)!= self.journeys):
				random.shuffle(candidate)
				cap = 0 
				nodesvisited = 0
				journey = []
				route = []

				for i in candidate:
					cap += self.nodecap[i-1]

					if cap > max_cap:
						route.append(journey)
						journey = []
						cap = self.nodecap[i-1]

					journey.append(i)
					nodesvisited += 1

				if journey != []:
					route.append(journey)

			self.population.append(route)

		self.initial_solution = self.population[1]
		self.initial_solution_fitness = self.individual_fitness(self.initial_solution)


	# Calculate euclidean distance between nodes
	def calculate_distance(self, node1, node2):
		distance = math.sqrt( (self.nodecoords[node1][0] - self.nodecoords[node2][0])**2 +  (self.nodecoords[node1][1] - self.nodecoords[node2][1])**2 )
		return distance


	# Calculate fitness for every member in the population
	# By default fitness is being maximized, for minimization the value's reciprocal is stored
	def fitness_func(self): 	
		self.fitness = []

		for i in range(len(self.population)):
			if (not self.is_valid(self.population[i])):
				self.population[i] = self.shuffle_route(self.population[i])

			individual_fitness = self.individual_fitness(self.population[i])
			self.fitness.append(1/individual_fitness)

		maximum = 0
		idx = 0

		for i in range(0,len(self.fitness)):
			if self.fitness[i] > maximum:
				maximum = self.fitness[i]
				idx = i

		self.max_solution = self.population[idx]
		self.max_fitness = maximum


	# Calculate fitness of an individual i.e. total distance travelled 
	def individual_fitness(self, candidate):
		individual_fitness = 0

		for j in range(len(candidate)):
			first_journey_node = candidate[j][0]
			individual_fitness += self.calculate_distance(0, first_journey_node)	

			for k in range(len(candidate[j]) - 1):
				idx1, idx2 = candidate[j][k], candidate[j][k+1]
				individual_fitness += self.calculate_distance(idx1, idx2)
				
			last_journey_node = candidate[j][-1]
			individual_fitness += self.calculate_distance(last_journey_node, 0) 

		return individual_fitness


	# Check if route is valid i.e. meets capacity requierments
	def is_valid(self, candidate):
		cap = 0
		for journey in candidate:
			for j in journey:
				cap += self.nodecap[j-1]

			if cap > self.max_cap:
				return False
			cap = 0

		return True


	# If route does not meet capacity requirements, shuffle it stochastically to ensure it does
	def shuffle_route(self, org_candidate):

		candidate = []
		for i in org_candidate:
			for j in i:
				candidate.append(j)

		random.shuffle(candidate)
		route = []

		while (len(route)!= self.journeys):
			random.shuffle(candidate)
			cap = 0 
			nodesvisited = 0
			journey = []
			route = []

			for i in candidate:
				cap += self.nodecap[i-1]

				if cap > max_cap:
					route.append(journey)
					journey = []
					cap = self.nodecap[i-1]

				journey.append(i)
				nodesvisited += 1

			if journey != []:
				route.append(journey)

		return route


	# Select either parents or survivors to proceed to the next stage through a variety of selection schemes
	# Num is the size of individuals to select and scheme defines the type of selection scheme
	def selection(self, num, scheme, survivor): 
		self.selected_candidates = []
		
		if scheme == "Random":
			
			for i in range(num):
				random_num = random.randrange(len(self.population))
				self.selected_candidates.append(self.population[random_num])
		

		elif scheme == "FPS":

			norm_fitness = []
			sum_fitness = sum(self.fitness)

			for i in self.fitness:
				norm_fitness.append(i/sum_fitness)

			cum_prob = []
			a = 0
			
			for k in norm_fitness:
				a += k
				cum_prob.append(a)
			
			np.array(cum_prob)
			
			for i in range(num):
				random_num = random.uniform(0,1)
				
				for j in range(len(cum_prob)-1):
					if (random_num > cum_prob[j] and random_num <= cum_prob[j+1]):
						self.selected_candidates.append(self.population[j+1])
						break
					elif (random_num <= cum_prob[0]):
						self.selected_candidates.append(self.population[0])
						break
		

		elif scheme == "Truncation":
			fitness_copy = np.copy(self.fitness)

			numbered_fitness_copy = list(enumerate(fitness_copy))
			numbered_fitness_copy.sort(key=lambda x:x[1])
			numbered_fitness_copy = numbered_fitness_copy[::-1]

			for i in range(num):
				index = numbered_fitness_copy[i][0]
				self.selected_candidates.append(self.population[index])
			

		elif scheme == "RBS":
			
			fitness_copy = np.copy(self.fitness)
			numbered_fitness_copy = list(enumerate(fitness_copy))
			numbered_fitness_copy.sort(key=lambda x:x[1])

			rank_dic = dict()
			rank_lst = []
			for i in range(1,len(self.fitness)+1):
				rank_dic[i] = numbered_fitness_copy[i-1]
				rank_lst.append(i)

			norm_rank_dic = dict()
			norm_lst = []  
			for j in range(1,len(self.fitness)+1):
				norm_rank_dic[rank_lst[j-1]/sum(rank_lst)] = rank_dic[j]
				norm_lst.append(rank_lst[j-1]/sum(rank_lst))
			
			cum_lst = []
			a = 0
			for p in norm_lst:
				a+=p
				cum_lst.append(a)
			
			for k in range(num):
				random_num2 = random.uniform(0,1)
				
				for m in range (len(cum_lst)-1):

					if (random_num2 > cum_lst[m] and random_num2 <= cum_lst[m+1]):
						fitness = norm_rank_dic[norm_lst[m+1]][1]
						index = norm_rank_dic[norm_lst[m+1]][0]
						self.selected_candidates.append(self.population[index])
						break
					elif random_num2 <= cum_lst[0]:
						fitness = norm_rank_dic[norm_lst[0]][1]
						index = norm_rank_dic[norm_lst[0]][0]
						self.selected_candidates.append(self.population[index])
						break
		

		elif scheme == 'BT':
			for i in range(num):
				cand_dic = dict()
				
				for j in range(2):
					rand = random.randrange(len(self.population))
					cand_dic[rand] = self.fitness[rand]
				
				max_index = max(cand_dic, key=cand_dic.get)
				self.selected_candidates.append(self.population[max_index])
	   

		else:
			print("Please enter a valid scheme")
		
		
		if survivor == True:
			self.population = np.array(self.selected_candidates)
		else:
			self.selected_candidates = np.array(self.selected_candidates)
   

	# Create new generation of off springs
	def crossover(self):
		selected_candidates_copy = np.copy(self.selected_candidates)
		selected_candidates_copy = selected_candidates_copy.tolist()
		temp_offsprings = []
		
		#Two point crossover
		for i in range(0,len(selected_candidates_copy)-1,2):
			parent1 = selected_candidates_copy[i]
			parent2 = selected_candidates_copy[i+1]

			child1 = self.create_child(parent1, parent2)
			child2 = self.create_child(parent2, parent1)
			
			offspring1 = np.array(child1)
			offspring2 = np.array(child2)
			temp_offsprings.append(offspring1)
			temp_offsprings.append(offspring2)


		self.offsprings = np.array(temp_offsprings)

   	
   	# Create child by performing two point crossover on parents
	def create_child(self,org_parent1, org_parent2):

		parent1 = []
		parent2 = []

		for i in org_parent1:
			for j in i:
				parent1.append(j)

		for i in org_parent2:
			for j in i:
				parent2.append(j)

		crossover_point1 = int(random.randint(0,len(parent1)))
		crossover_point2 = int(random.randint(0,len(parent1)))

		while crossover_point2 == crossover_point1:
			crossover_point2 = int(random.randint(0,len(parent1)))

		start = min(crossover_point1, crossover_point2)
		end = max(crossover_point1, crossover_point2)
		insertion = parent1[start:end]
		stored = []
		child = [None]*len(parent1)

		j = 0
		for i in range(len(parent1)):

			if i>=start and i<end:
				child[i] = insertion[i-start]
			else:
				while (parent2[j] in insertion or parent2[j] in stored or parent2[j] in child) and j<len(parent2):
					j += 1
				child[i] = parent2[j]
				stored.append(parent2[j])

		# Ensure that child created meets capacity requirements
		route = []

		candidate = child
		attempts = 0
		while (len(route)!= self.journeys):

			if (attempts>2):
				random.shuffle(candidate)

			cap = 0 
			journey = []
			route = []

			for i in candidate:
				cap += self.nodecap[i-1]

				if cap > max_cap:
					route.append(journey)
					journey = []
					cap = self.nodecap[i-1]

				journey.append(i)

			if journey != []:
				route.append(journey)

			attempts += 1

		org_child = route

		return org_child


	# mutate the offsprings by either swaps or inversion
	def mutation(self):
		for i in range(len(self.offsprings)):
			random_num = random.uniform(0,1)

			if random_num < self.mutation_rate:
				repeat = True

				while (repeat):
					repeat = False
					random_journey_idx = random.randint(0, self.journeys-1)
					while (len(self.offsprings[i][random_journey_idx])<2):
						random_journey_idx = random.randint(0, self.journeys-1)

					journey_length = len(self.offsprings[i][random_journey_idx])
					random_idx = random.randint(0, journey_length-1)
					random_idx2 = random.randint(0, journey_length-1)    
					counter = 0
					self.offsprings[i][random_journey_idx][random_idx], self.offsprings[i][random_journey_idx][random_idx2] =  self.offsprings[i][random_journey_idx][random_idx2], self.offsprings[i][random_journey_idx][random_idx]
					
					if (not self.is_valid(self.offsprings[i])):
						repeat = True

		self.population = np.concatenate((self.population,self.offsprings))



def solve_vrp(generations, iterations, parent_size, population_size, nodelist, nodecoords, nodecap, journeys, max_cap, mutation_rate = 0.5):

	total_BFS = []
	total_ASF = []
	best_solution_set = [None]*iterations
	best_fitness_set = [None]*iterations

	for i in range(generations):
		total_BFS.append([])
		total_ASF.append([])

	# Algorithm called for different iterations to average out the result values
	for i in range(iterations):
		vrp_prob = EA(nodelist, nodecoords, nodecap, population_size, journeys, max_cap, mutation_rate)
		print('iteration:', i)
		for j in range(generations):

			vrp_prob.fitness_func()

			max_generation_fitness = vrp_prob.max_fitness#

			if len(total_BFS[j])==0:
				total_BFS[j].append(max_generation_fitness)
			elif max_generation_fitness >= total_BFS[j][-1]:
				total_BFS[j].append(max_generation_fitness)
			else:
				total_BFS[j].append(total_BFS[j][-1])

			avg_generation_fitness = (np.average(vrp_prob.fitness) + sum(total_ASF[j])) / (len(total_ASF[j])+1)
			total_ASF[j].append(avg_generation_fitness)

			# Perform genetic algorithm operations
			vrp_prob.selection(parent_size, "BT", False) 
			vrp_prob.crossover()
			vrp_prob.mutation()
			vrp_prob.fitness_func()
			vrp_prob.selection(population_size, "Truncation", True) #survivor selection hence last argument is True
			
			# Calculate initial solution/fitness and best solution/fitness
			initial_solution = vrp_prob.initial_solution
			initial_fitness = vrp_prob.initial_solution_fitness
			best_solution_set[i] = vrp_prob.max_solution
			best_fitness_set[i] = 1/max_generation_fitness

	
	# Calculate best fitness and best solution	
	best_fitness = 1/(max(total_BFS[-1]))
	best_solution_set = list(best_solution_set)
	best_sol_idx = best_fitness_set.index(best_fitness)
	best_solution = list(best_solution_set[best_sol_idx])

	print(best_solution)
	print(best_fitness)

	empty = []
	for i in vrp_prob.max_solution:
	    if i not in empty:
	        empty.append(i)
	    else:
	        print('DUPLICATION')

	generations_list = []
	for i in range(1,generations+1):
	    generations_list.append(i)

	# Plot the best fitness and average fitness values
	BFS_gen = []
	ASF_gen = []
	gen_file = open("gen_BFS_ASF.txt","w")
	for i in range(0,generations):
	    BFS_gen.append(1/np.average(total_BFS[i]))
	    ASF_gen.append(1/np.average(total_ASF[i]))
	    gen_file.write(str(generations_list[i]))
	    gen_file.write(" | ")
	    gen_file.write(str(1/np.average(total_BFS[i])))
	    gen_file.write(" | ")
	    gen_file.write(str(1/np.average(total_ASF[i])))
	    gen_file.write("\n")

	BFS_iter_gen = []
	ASF_iter_gen = []
	BFS_iter = open("BFS_iteration_gen.txt","w")
	ASF_iter = open("ASF_iteration_gen.txt","w") 
	for i in range(0,generations):
	  temp_BFS = []
	  temp_ASF = []
	  BFS_iter.write(str(generations_list[i]))
	  ASF_iter.write(str(generations_list[i]))
	  for j in range(0,iterations):
	    temp_BFS.append(1/total_BFS[i][j])
	    temp_ASF.append(1/total_ASF[i][j])
	    BFS_iter.write(" | ")
	    ASF_iter.write(" | ")
	    BFS_iter.write(str(1/total_BFS[i][j]))
	    ASF_iter.write(str(1/total_ASF[i][j]))
	  BFS_iter_gen.append(temp_BFS)
	  ASF_iter_gen.append(temp_ASF)
	  BFS_iter.write(" | ")
	  BFS_iter.write(str(np.average(temp_BFS)))
	  BFS_iter.write("\n")
	  ASF_iter.write(" | ")
	  ASF_iter.write(str(np.average(temp_ASF)))
	  ASF_iter.write("\n")
	BFS_iter.close()
	ASF_iter.close()

	plt.close('all')
	plt.title('Fitness vs Generation (BT and Truncation)')
	plt.plot(generations_list, BFS_gen, label="BFS")
	plt.ylabel('Fitness')
	plt.xlabel('Generations')
	plt.plot(generations_list, ASF_gen, label="ASF")
	plt.legend(framealpha=1, frameon=True);
	plt.savefig('BTandTruncation.png')
	plt.show()

	return best_fitness

# datasets -> https://github.com/giulianoxt/vehicle-routing-aco/blob/master/augerat-a/A-n32-k5.vrp
# Read the data file
infile = open('A-n32-k5.txt', 'r')

# Read instance header
infile.readline().strip().split()
infile.readline().strip().split()
infile.readline().strip().split()
dim = int(infile.readline().strip().split()[-1])
infile.readline().strip().split()
max_cap = int(infile.readline().strip().split()[-1])
journeys = int(infile.readline().strip().split()[-1])
infile.readline().strip().split()

nodelist = []
nodecoords = []
nodecap = []

for i in range(0, dim):
	line = infile.readline()       
	
	if (line != ''):
		n,x,y = line.strip().split()
		if (float(n)==1):
			nodecoords.append([float(x), float(y)])
		else:
			nodelist.append(int(n)-1)
			nodecoords.append([float(x), float(y)])
	else:
		break

infile.readline()

for i in range(0, dim):
	line = infile.readline()       

	if (line != ''):
		n,c = line.strip().split()
		if (float(n)==1):
			pass
		else:
			nodecap.append(float(c))
	else:
		break

# Close input file
infile.close()


# Define parameter values
generations = 400
iterations = 1
parent_size = 850
population_size = 1500
mutation_rate = 0.5

solve_vrp(generations, iterations, parent_size, population_size, nodelist, nodecoords, nodecap, journeys, max_cap, mutation_rate)

# mutation = []
# fitness_set = []

# for i in range(0, 1000, 20):
# 	print(i/1000)
# 	fitness = solve_vrp(generations, iterations, parent_size, population_size, nodelist, nodecoords, nodecap, journeys, max_cap, i/1000)
# 	fitness_set.append(fitness)
# 	mutation.append(i/1000)


# plt.close('all')
# plt.title('Fitness vs Mutation Rate')
# plt.plot(mutation, fitness_set)
# plt.ylabel('Fitness')
# plt.xlabel('Mutation Rate')
# plt.show()

###########################################################################################
###########################################################################################
###########################################################################################
################  Code by Abdullah Siddiqui (ms03586@st.habib.edu.pk) #####################
###########################################################################################
###########################################################################################
###########################################################################################