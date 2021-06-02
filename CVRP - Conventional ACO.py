###########################################################################################
###########################################################################################
###########################################################################################
################  Code by Abdullah Siddiqui (ms03586@st.habib.edu.pk) #####################
###########################################################################################
###########################################################################################
###########################################################################################

import getopt
import math
import random
import numpy
from functools import reduce
import sys
from scipy.spatial.distance import pdist, squareform
import numpy as np
import sys
import networkx as nx
import random
import matplotlib.pyplot as plt
import math
import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import math


# obtain the route traversed by one ant i.e. one possible solution while considering capacity requirements
def single_ant_route(nodelist, graph_edges, nodecap, pheromone_trail, max_cap, alpha, beta):
        nodes = nodelist.copy()
        route = []

        while len(nodes)!=0 : 
                trip = [] 

                # ant initialized at a random node      
                node = int(np.random.choice(nodes))     
                ant_cap = max_cap - nodecap[node]
                trip.append(int(node))
                nodes.remove(int(node))

                while len(nodes) != 0:
                        # define transition probabilities from node i to j 
                        transition_probabilities = []
                        for vertex in nodes:
                                transition_probabilities.append(  (( pheromone_trail[ int(min(vertex, node)), int(max(vertex, node)) ] )**alpha) * (Q/graph_edges[ int(min(vertex, node)), int(max(vertex, node)) ])**beta )  

                        transition_probabilities = transition_probabilities/np.sum(transition_probabilities) 

                        cumsum = []
                        a = 0
                        for k in transition_probabilities:
                                a += k 
                                cumsum.append(a)

                        # select node depending on transition probabilities
                        random_num = random.uniform(0,1)

                        if (len(cumsum) == 1):
                                node = nodes[0] 

                        elif len(cumsum) == 2:
                            if random_num < 0.5:
                                node = nodes[0]
                            else:
                                node = nodes[1]

                        else:
                            for i in range(len(cumsum) - 1):
                                if (random_num > cumsum[i]  and  random_num <= cumsum[i+1]):
                                        node = nodes[i+1]
                                        break
                                elif (random_num <= cumsum[0]):
                                        node = nodes[0]
                                        break

                        ant_cap = ant_cap - nodecap[node]

                        if (ant_cap>0):
                            if node not in trip:
                                trip.append(int(node))
                            if (int(node) in nodes):
                                    nodes.remove(int(node))
                        else:
                            break

                route.append(trip)

        return route


# obtain fitness value for the specified route i.e. total distance traversed by an ant to explore all nodes
def fitness_function(route, graph_edges):
        fitness = 0
        complete_route = [1]

        for trip in route:
                for node in trip:
                        complete_route.append(int(node))
                complete_route.append(1)

        for i in range(len(complete_route)-1):
                fitness += graph_edges[min(int(complete_route[i]), int(complete_route[i+1])), max(int(complete_route[i]), int(complete_route[i+1])) ]

        return fitness


# define dictionaries detailing distances between nodes and the pheromone trail between nodes
def create_graph(nodelist, nodecoords, nodecap):
        graph_edges = {}
        pheromone_trail = {}
        for i in nodelist:
                for j in nodelist: 
                        distance = np.sqrt( (nodecoords[i][0] - nodecoords[j][0])**2 + (nodecoords[i][1] - nodecoords[j][1])**2 ) 
                        graph_edges[ int(min( int(i), int(j) )), int(max(int(i), int(j))) ] = distance

                        if i!=j:
                                pheromone_trail[ int(min(int(i), int(j) )), int(max(int(i), int(j) )) ] = 1

        return graph_edges, pheromone_trail


# update pheromone trail based on standard formula
def update_pheromone_trail(pheromone_trail, routes, best_solution, rho, Q, sigma):
        fitness = []
        for value in routes:
                fitness.append(value[1])

        avg_fitness = sum(fitness)/len(routes)

        # standard pheromone update equation factoring in the evaporation rate
        for key, value in pheromone_trail.items():
                pheromone_trail[key] = (rho + Q/avg_fitness)*value

        routes.sort(key = lambda x:x [1])

        if (best_solution != []):

                if (routes[0][1] < best_solution[1]):
                        best_solution = routes[0]

                # consider effect of elitist ants on updating pheromone value
                for trip in best_solution[0]:
                        for i in range(len(trip) - 1):
                                pheromone_trail[ (min(trip[i], trip[i+1]), max(trip[i], trip[i+1])) ] += sigma/best_solution[1]

        else:

                best_solution = routes[0]

        # incorporate effect of elitist ant strategy on pheromone trail
        for i in range(sigma):
                trips = routes[i][0]
                fitness = routes[i][1]

                for trip in trips:
                        for j in range(len(trip) - 1):
                                pheromone_trail[ (min(trip[j], trip[j+1]), max(trip[j], trip[j+1])) ] += (sigma - (j+1)) / (fitness**(j+1))                             

        return best_solution


# solve the vrp problem using aco
def aco_solver(nodelist, nodecoords, nodecap, max_cap, iterations, alpha, beta, Q, sigma, rho, Q, ants):
        graph_edges, pheromone_trail = create_graph(nodelist, nodecoords, nodecap)
        best_solution = []
        BFS = []
        AFS = []

        # run algorithm for specified number of iterations to find optimal solution
        for i in range(iterations):
                routes = []
                avg_fitness = []

                # create a colony of ants i.e. possible solutions in each iteration
                for j in range(ants):
                        route = single_ant_route(nodelist, graph_edges, nodecap, pheromone_trail, max_cap, alpha, beta)
                        fitness = fitness_function(route, graph_edges)
                        routes.append([route, fitness])
                        avg_fitness.append(fitness)

                average_fitness = sum(avg_fitness)/len(avg_fitness)
                best_solution = update_pheromone_trail(pheromone_trail, routes, best_solution, rho, Q, sigma)
                BFS.append(best_solution[1])
                AFS.append(average_fitness)
                print("Best: ", best_solution[1], " Avg: ", average_fitness)

        return best_solution, BFS, AFS


# Read data file
infile = open('A-n37-k5.txt', 'r')

# Read instance header
(infile.readline().strip().split())
(infile.readline().strip().split())
(infile.readline().strip().split())
dim = int(infile.readline().strip().split()[-1])
infile.readline().strip().split()
max_cap = int(infile.readline().strip().split()[-1])
journeys = int(infile.readline().strip().split()[-1])
infile.readline().strip().split()

nodecoords = {}
nodecap = {}
nodelist = []

for i in range(0, dim):
        line = infile.readline()       
        if (line != ''):
                n,x,y = line.strip().split()

                if (n!=1):
                        nodelist.append(int(n))

                nodecoords[int(n)] = [float(x), float(y)]
        else:
                break

infile.readline()

for i in range(0, dim):
        line = infile.readline()       
        if (line != ''):
                n,c = line.strip().split()
                nodecap[int(n)] = int(c)
        else:
                break

# Define parameter values
alpha = 5
beta = 8
sigma = 6
rho = 0.9 
Q = 10
iterations = 1000
ants = dim

# solve vrp using aco
best_sol, BFS, AFS = aco_solver(nodelist, nodecoords, nodecap, max_cap, iterations, alpha, beta, Q, sigma, rho, Q, ants)

# plot the results obtained
plt.title('Fitness vs Iteration')
plt.plot([i for i in range(iterations)], BFS, label="BFS")
plt.ylabel('Fitness')
plt.xlabel('Iterations')

plt.title('Fitness vs Iteration')
plt.plot([i for i in range(iterations)], AFS, label="AFS")
plt.ylabel('Fitness')
plt.xlabel('Iterations')
plt.legend(framealpha=1, frameon=True);
plt.show()

###########################################################################################
###########################################################################################
###########################################################################################
################  Code by Abdullah Siddiqui (ms03586@st.habib.edu.pk) #####################
###########################################################################################
###########################################################################################
###########################################################################################

        