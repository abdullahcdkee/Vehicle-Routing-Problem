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
import veroviz as vrv
import os
import pandas as pd

# Define directories and keys for the Veroviz library
CESIUM_DIR = os.environ['CESIUMDIR']
DATA_PROVIDER = 'ORS-online'
DATA_PROVIDER_ARGS = {'APIkey': os.environ['ORSKEY'], 'databaseName': None}

# Define boundary in which to generate the nodes
myBoundary = [[42.914949262700084, -78.90020370483398], 
              [42.871938424448466, -78.89556884765626], 
              [42.875083507773496, -78.82158279418947], 
              [42.91293772859624, -78.81729125976564]]

numCustomers = 26

vrv.createLeaflet(boundingRegion = myBoundary)

# Create single depot node:
nodesDF = vrv.generateNodes(
    numNodes         = 1,
    startNode        = 0,
    nodeType         = 'depot',
    leafletColor     = 'red',
    leafletIconType  = 'home',
    nodeDistrib      = 'uniformBB', 
    nodeDistribArgs  = { 'boundingRegion': myBoundary },
    dataProvider     = DATA_PROVIDER,
    dataProviderArgs = DATA_PROVIDER_ARGS,
    snapToRoad       = True)

# Create customer nodes:
nodesDF = vrv.generateNodes(
    initNodes       = nodesDF,
    numNodes        = numCustomers,
    startNode       = 1,
    nodeType        = 'customer',
    leafletColor    = 'blue',
    leafletIconType = 'star',
    nodeDistrib     = 'uniformBB', 
    nodeDistribArgs = { 'boundingRegion': myBoundary },
    dataProvider     = DATA_PROVIDER,
    dataProviderArgs = DATA_PROVIDER_ARGS,
    snapToRoad      = True)

# Create a map that shows the boundary and the nodes:
vrv.createLeaflet(nodes          = nodesDF, 
                  boundingRegion = myBoundary,
                  mapBackground  = 'Arcgis Roadmap')

[time, dist] = vrv.getTimeDist2D(nodes            = nodesDF, 
                                 matrixType       = 'all2all',
                                 routeType        = 'fastest', 
                                 dataProvider     = DATA_PROVIDER,
                                 dataProviderArgs = DATA_PROVIDER_ARGS)


# Class that represents the genetic algorithm and associated functions
class EA:
    def __init__(self, nodelist, nodecoords, nodecap, size, journeys, max_cap, dist, time = [], consider_time = False):
        self.nodelist = nodelist
        self.nodecoords = nodecoords
        self.nodecap = nodecap
        self.population_size = size
        self.mutation_rate = 0.5
        self.max_time = 10000 #seconds
        self.population = [] 
        self.dist = dist
        self.journeys = journeys
        self.max_cap = max_cap
        self.mutation_valid = 0
        self.mutation_invalid = 0
        self.initial_solution = []
        self.consider_time = consider_time

        self.initialize_population()


    # Initialize random chromosomes for the population
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


    # Calculate time of one trip based on the values calculated by the map interface
    def calculate_time(self, node1, node2):
        time_taken = time[(node1, node2)]
        return time_taken


    # Calculate time of entire route based on the values calculated by the map interface
    def calculate_total_time(self, candidate):

        if (time == []):
            return 0

        time_set = []
        total_time = 0
        for journey in candidate:
            total_time = 0
            total_time += self.calculate_time(0, journey[0])

            for j in range(len(journey)-1):
                total_time += self.calculate_time(journey[j], journey[j+1])
            total_time += self.calculate_time(journey[-1], 0)

            time_set.append(total_time)

        return max(time_set)


    # Calculate distance of entire route based on the values calculated by the map interface
    def calculate_distance(self, node1, node2):
        distance = dist[(node1, node2)]
        return distance


    # Calculate fitness for every member in the population
    # By default fitness is being maximized, for minimization the value's reciprocal is stored
    def fitness_func(self): 
        self.fitness = []

        for i in range(len(self.population)):
            if (not self.consider_time):
                if ( not self.is_valid(self.population[i]) ):
                    self.population[i] = self.shuffle_route(self.population[i])
            else:
                if ( (not self.is_valid(self.population[i])) or (not self.is_time_valid(self.population[i]))):
                    self.population[i] = self.shuffle_route(self.population[i])

            if (not self.consider_time):     
                individual_fitness = self.individual_fitness(self.population[i])
            else:
                individual_fitness = self.individual_fitness(self.population[i])  + ( 30*self.calculate_total_time(self.population[i]) )

            self.fitness.append(1/individual_fitness)

        maximum = 0
        idx = 0
        for i in range(0,len(self.fitness)):
            if self.fitness[i] > maximum:
                maximum = self.fitness[i]
                idx = i

        self.max_solution = self.population[idx]
        self.max_fitness = maximum


    # Calculate fitness of an individual i.e. total distance travelled and time taken weighted together
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


    # Check if route is valid i.e. meets time window requierments
    def is_time_valid(self, candidate):
        if (time == []):
            return True

        total_time = 0
        for journey in candidate:

            total_time = 0
            total_time += self.calculate_time(0, journey[0])

            for j in range(len(journey)-1):
                total_time += self.calculate_time(journey[j], journey[j+1])
            total_time += self.calculate_time(journey[-1], 0)

            if total_time > self.max_time:
                return False

        return True


    # If route does not meet capacity or time requirements, shuffle it stochastically to ensure it does
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
        repeat = False
        while (len(route)!= self.journeys or repeat):
            if (attempts>5):
                random.shuffle(candidate)

            if attempts>15:
                break

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

            if (self.consider_time):
                if (not self.is_time_valid(route)):
                    repeat = True

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
                    random_journey_idx = random.randint(0, len(self.offsprings[i])-1)

                    while (len(self.offsprings[i][random_journey_idx])<2):
                        random_journey_idx = random.randint(0, self.journeys-1)

                    journey_length = len(self.offsprings[i][random_journey_idx])
                    random_idx = random.randint(0, journey_length-1)
                    random_idx2 = random.randint(0, journey_length-1)    
                    self.offsprings[i][random_journey_idx][random_idx], self.offsprings[i][random_journey_idx][random_idx2] =  self.offsprings[i][random_journey_idx][random_idx2], self.offsprings[i][random_journey_idx][random_idx]
                    
                    if (not self.is_valid(self.offsprings[i])):
                        repeat = True

        self.population = np.concatenate((self.population,self.offsprings))


# function solves the vehicle routing problem
def solve_vrp(generations, iterations, parent_size, population_size, nodelist, nodecoords, nodecap, journeys, max_cap, dist, time = [], consider_time = False, max_time = 5000):
    total_BFS = []
    total_ASF = []
    best_solution_set = [None]*iterations
    best_fitness_set = [None]*iterations

    for i in range(generations):
        total_BFS.append([])
        total_ASF.append([])

    for i in range(iterations):
        vrp_prob = EA(nodelist, nodecoords, nodecap, population_size, journeys, max_cap, dist, time, consider_time)
        
        for j in range(generations):

            if ((consider_time) and (j==290)):
                vrp_prob.max_time = max_time

            vrp_prob.fitness_func()
            max_generation_fitness = vrp_prob.max_fitness

            if len(total_BFS[j])==0:
                total_BFS[j].append(max_generation_fitness)
            elif max_generation_fitness >= total_BFS[j][-1]:
                total_BFS[j].append(max_generation_fitness)
            else:
                total_BFS[j].append(total_BFS[j][-1])

            avg_generation_fitness = (np.average(vrp_prob.fitness) + sum(total_ASF[j])) / (len(total_ASF[j])+1)
            total_ASF[j].append(avg_generation_fitness)

            vrp_prob.selection(parent_size, "BT", False) 
            vrp_prob.crossover()
            vrp_prob.mutation()
            vrp_prob.fitness_func()
            vrp_prob.selection(population_size, "Truncation", True) #survivor selection hence last argument is True
            
            initial_solution = vrp_prob.initial_solution
            initial_fitness = vrp_prob.initial_solution_fitness
            best_solution_set[i] = vrp_prob.max_solution
            best_fitness_set[i] = 1/max_generation_fitness
    
    best_fitness = 1/(max(total_BFS[-1]))
    best_solution_set = list(best_solution_set)
    best_sol_idx = best_fitness_set.index(best_fitness)
    best_solution = list(best_solution_set[best_sol_idx])

    #best time
    best_time = vrp_prob.calculate_total_time(best_solution)
    #best distance
    best_dist = vrp_prob.individual_fitness(best_solution)

    print("Best Fitness: ",best_fitness)
    print("Best Time Fitness: ", best_time)
    print("Best Distance Fitness: ", best_dist)

    empty = []
    for i in vrp_prob.max_solution:
        if i not in empty:
            empty.append(i)

    generations_list = []
    for i in range(1,generations+1):
        generations_list.append(i)

    # plot the best and average fitness values
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
    # plt.show()

    initial_route = [0]
    for i in initial_solution:
        for j in i:
            initial_route.append(j)
        initial_route.append(0)

    best_route = [0]
    for i in best_solution:
        for j in i:
            best_route.append(j)
        best_route.append(0)

    return([initial_solution, best_solution, best_time])


# load data file to allocate capacity to the nodes and vehicles
infile = open('data2.txt', 'r')

# Read instance header
journeys, max_cap = infile.readline().strip().split()
journeys = int(journeys)
max_cap = int(max_cap)

nodelist = []
nodecoords = []
nodecap = []

N = 1000
for i in range(0, numCustomers+1):        
    line = infile.readline()

    if (line != ''):
        n,x,y,c = line.strip().split()
        if (float(c)==0):
            nodecoords.append([float(x), float(y)])
        else:
            nodelist.append(int(n)-1)
            nodecoords.append([float(x), float(y)])
            nodecap.append(float(c))

    else:
        break

total_capacity = sum(nodecap)
journeys = total_capacity//100 + 1

# Close input file
infile.close()


# Define parameter values
generations = 300
iterations = 1
parent_size = 150
population_size = 300

# Solve CVRP
random_route, dist_opt_route, dist_opt_sol_time = solve_vrp(generations, iterations, parent_size, population_size, nodelist, nodecoords, nodecap, journeys, max_cap, dist)
print("Dist Opt: ", dist_opt_route)
print()

# Specify the CVRP-TW to optimize itself with a time window of 5 minutes less than the time obtained in CVRP
max_time = dist_opt_sol_time - 300
# Solve the CVRP-TW problem
random_route, time_opt_route, time_opt_sol_time = solve_vrp(generations, iterations, parent_size, population_size, nodelist, nodecoords, nodecap, journeys, max_cap, dist, time, True, max_time)
print("Time Opt: ", time_opt_route)
print()


# Visualize the improved routes via Cesium and Folium
def vrvSolver(nodesDF, dist, time, truck_route):
    import pandas as pd

    network = []
    count = 0
    truck = 'truck'
    for path in truck_route:
        truck = truck + str(count)
        route = {truck:[0]}

        for node in path:
            route[truck].append(node)
        route[truck].append(0)
        count += 1
        network.append(route)

    # Define configuration of the different vehicles used to deliver the goods  
    configs = {'truck0': {
                    'vehicleModels': ['veroviz/models/ub_truck.gltf'],
                    'leafletColor': 'blue',
                    'cesiumColor': 'Cesium.Color.BLUE',
                    'packageModel': 'veroviz/models/box_blue.gltf',
                    'modelScale': 100,
                    'minPxSize': 75 }, 
                'truck01': {
                    'vehicleModels': ['veroviz/models/car_green.gltf'],
                    'leafletColor': 'green',
                    'cesiumColor': 'Cesium.Color.YELLOWGREEN',
                    'packageModel': 'veroviz/models/rectangle_green.gltf',
                    'modelScale': 100,
                    'minPxSize': 75 }, 
                'truck012': {
                    'vehicleModels': ['veroviz/models/car_red.gltf'],
                    'leafletColor': 'red',
                    'cesiumColor': 'Cesium.Color.RED',
                    'packageModel': 'veroviz/models/rectangle_red.gltf',
                    'modelScale': 100,
                    'minPxSize': 75 }, 
                'truck0123': {
                    'vehicleModels': ['veroviz/models/car_blue.gltf'],
                    'leafletColor': 'purple',
                    'cesiumColor': 'Cesium.Color.LIGHTSKYBLUE',
                    'packageModel': 'veroviz/models/rectangle_blue.gltf',
                    'modelScale': 100,
                    'minPxSize': 75 }, 
                'truck01234': {
                    'vehicleModels': ['veroviz/models/ub_truck.gltf'],
                    'leafletColor': 'yellow',
                    'cesiumColor': 'Cesium.Color.YELLOW',
                    'packageModel': 'veroviz/models/box_blue.gltf',
                    'modelScale': 100,
                    'minPxSize': 75 }
              }
        
    serviceTime = 30 # seconds    
    
    # Initialize an empty "assignments" dataframe.  
    assignmentsDF = vrv.initDataframe('assignments')
    

    time_set = []

    for route in network: 
        startTime = 0   # line moves outside the loop
        for vehicle in route:
            for i in list(range(0, len(route[vehicle])-1)):
                startNode = route[vehicle][i]
                endNode   = route[vehicle][i+1]

                startLat  = nodesDF[nodesDF['id'] == startNode]['lat'].values[0]
                startLon  = nodesDF[nodesDF['id'] == startNode]['lon'].values[0]
                endLat    = nodesDF[nodesDF['id'] == endNode]['lat'].values[0]
                endLon    = nodesDF[nodesDF['id'] == endNode]['lon'].values[0]

                if ((vehicle == 'drone') and (startNode == 0)):
                    # Use the 3D model of a drone carrying a package
                    myModel = configs[vehicle]['vehicleModels'][1]
                else:
                    # Use the 3D model of either a delivery truck or an empty drone
                    myModel = configs[vehicle]['vehicleModels'][0]

                if (vehicle == 'truck' or vehicle == 'truck0' or vehicle == 'truck01' or vehicle == 'truck012' or vehicle == 'truck0123' or vehicle == 'truck01234'):
                    # Get turn-by-turn navigation for the truck, as it travels
                    # from the startNode to the endNode:
                    shapepointsDF = vrv.getShapepoints2D(
                        # odID           = odID,
                        objectID         = vehicle, 
                        modelFile        = myModel,
                        modelScale       = configs[vehicle]['modelScale'], 
                        modelMinPxSize   = configs[vehicle]['minPxSize'], 
                        startTimeSec     = startTime,
                        startLoc         = [startLat, startLon],
                        endLoc           = [endLat, endLon],
                        routeType        = 'fastest',
                        leafletColor     = configs[vehicle]['leafletColor'], 
                        cesiumColor      = configs[vehicle]['cesiumColor'], 
                        dataProvider     = DATA_PROVIDER,
                        dataProviderArgs = DATA_PROVIDER_ARGS) 
                else:
                    # Get a 3D flight profile for the drone:
                    shapepointsDF = vrv.getShapepoints3D(
                        objectID           = vehicle, 
                        modelFile          = myModel,
                        modelScale         = configs[vehicle]['modelScale'], 
                        modelMinPxSize     = configs[vehicle]['minPxSize'], 
                        startTimeSec       = startTime,
                        startLoc           = [startLat, startLon],
                        endLoc             = [endLat, endLon],
                        takeoffSpeedMPS    = 5,               
                        cruiseSpeedMPS     = 20,              
                        landSpeedMPS       = 3,               
                        cruiseAltMetersAGL = 100,              
                        routeType          = 'square',
                        cesiumColor        = configs[vehicle]['cesiumColor']) 

                # Update the assignments dataframe:
                assignmentsDF = pd.concat([assignmentsDF, shapepointsDF], ignore_index=True, sort=False)

                # Update the time
                startTime = max(shapepointsDF['endTimeSec'])

                # Add loitering for service
                assignmentsDF = vrv.addStaticAssignment(
                    initAssignments      = assignmentsDF, 
                    objectID             = vehicle, 
                    modelFile            = myModel, 
                    modelScale           = configs[vehicle]['modelScale'], 
                    modelMinPxSize       = configs[vehicle]['minPxSize'], 
                    loc                  = [endLat, endLon],
                    startTimeSec    = startTime,
                    endTimeSec = startTime + serviceTime)

                # Update the time again
                startTime = startTime + serviceTime

                # Add a package at all non-depot nodes:
                if (endNode != 0):
                    assignmentsDF = vrv.addStaticAssignment(
                        initAssignments      = assignmentsDF, 
                        objectID             = 'package %d' % endNode,
                        modelFile            = configs[vehicle]['packageModel'], 
                        modelScale           = 100, 
                        modelMinPxSize       = 35, 
                        loc                  = [endLat, endLon],
                        startTimeSec    = startTime,
                        endTimeSec = -1)

        time_set.append(startTime)

    return assignmentsDF

# Visualize CVRP
print("Distance Optimal Route:")
assignmentsDF = vrvSolver(nodesDF, dist, time, dist_opt_route)

# Visualize CVRP-TW
print("Time Optimal Route:")
solutionDF = vrvSolver(nodesDF, dist, time, time_opt_route)

print("Routes Calculated")

# Create a Leaflet map showing the nodes and the routes.
random_map = vrv.createLeaflet(nodes=nodesDF, arcs=assignmentsDF)
optimal_map = vrv.createLeaflet(nodes=nodesDF, arcs=solutionDF)

random_map.save('distance_map.html')
optimal_map.save('time_map.html')

# Create a Cesium movie showing the nodes, routes, and package deliveries.
vrv.createCesium(
    assignments = assignmentsDF, 
    nodes       = nodesDF, 
    startTime   = '10:00:00', 
    cesiumDir   = os.environ['CESIUMDIR'],
    problemDir  = 'Opt_Dist')

vrv.createCesium(
    assignments = solutionDF, 
    nodes       = nodesDF, 
    startTime   = '10:00:00', 
    cesiumDir   = os.environ['CESIUMDIR'],
    problemDir  = 'Opt_Time')


###########################################################################################
###########################################################################################
###########################################################################################
################  Code by Abdullah Siddiqui (ms03586@st.habib.edu.pk) #####################
###########################################################################################
###########################################################################################
###########################################################################################

