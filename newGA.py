from __future__ import division
from random import randint
import numpy as np
import Tree_Harry  
import pdb
import time

def t1():
    # Quick run test method to allow everything to be done in a single call.     
    
    # t = the normal return from the run method. 
    t = run(5, 0, 512, 20, 500000, 5, 2.5,50)
    # a = list of all nodes in the tree.
    a = (t[0].getLeavesList([]))
    # b = a sorted list of all nodes in the tree. 
    b = sorted(a, key = lambda x: x.genome, reverse = True)
    return [t,a,b]

def run(dims, minGenomeVal, maxGenomeVal, startPopulation, iterations, noOfResults, varienceThreshold, revisitAmount, startRes = None):
    """
    Create a number of individuals (i.e. a population).

    dims                : the number of dimensions
    minGenomeVal        : the minimum value of the genome space
    maxGenomeVal        : the maximum value of the genome space 
    startPopilation     : the number of individuals initially added before the evolutions
    iterations          : the number of genetic iterations to perform
    noOfResults         : determine how many results are top be returned
    varienceThreshold   : the varience threshold to determine when to change resolution
    revisitAmount       : the number of revisits needed before a resolution change can occur. 
    startRes            : if you want the resoltuion to start at anything but the minimim. 

    return              : the tree object, top n fitness genomes, top n varience cell values
    """    
    
    #Get the start time
    start_time = time.time()    
    
    # initialise the tree object
    tree = population(startPopulation, dims, minGenomeVal, maxGenomeVal, startRes)
    
    # perform the defined number of genetic iterations        
    for i in range(0, iterations):
        if i % 10000 == 0:
            print "{0:.0f}%".format(float(i)/float(iterations) * 100)
        
        tree = evolve(tree, dims, minGenomeVal, maxGenomeVal, varienceThreshold, revisitAmount)
    
    # get the top n Fitness / varience results   
    topFitness          = getTopFitness(tree, noOfResults)
    topVarienceResults  = getTopVarience(tree, noOfResults)
    
    # Check the same genome isn't in the tree more than once.     
    duplicates(tree.getLeavesList([]))
    
    # get the run time    
    end_time = ("--- %s seconds ---" % (time.time() - start_time))    
    
    # Return the tree object, top n fitness genomes, top n varience cell values
    return [tree, topFitness, topVarienceResults, end_time]

    
def individual(length, minGenomeVal, maxGenomeVal):
    """
    
    Generate a random genome.
    
    length          : the number of dimensions
    minGenomeVal    : the minimum value of the genome space
    maxGenomeVal    : the maximum value of the genome space 
    
    return          : a random genome
    
    """
    
    #Account for indexing starting at 0
    maxGenomeVal = maxGenomeVal - 1
    
    return [ randint(minGenomeVal,maxGenomeVal) for x in xrange(length) ]

def population(count, dim, minGenomeVal, maxGenomeVal, startRes):
    """
    Create a number of individuals (i.e. a population).

    count           : the number of individuals in the starting population
    dim             : the number of dimensions
    minGenomeVal    : the minimum value of the genome space
    maxGenomeVal    : the maximum value of the genome space 

    return          : the tree object
    """
    # Setup the initial feature space
    initialResolution = [maxGenomeVal]*dim
    
    # Account for indexing starting at 0
    maxGenomeVal = maxGenomeVal - 1

    # Initialise the tree.
    tree = Tree_Harry.Tree(res = initialResolution, carg = [None] * 3, child = [], parent = None, blockLocation = None, maxReso=None)
    
    # Half the entire features spaces resolution    
    halvedResolution = [x/2 for x in initialResolution]   
    
    if startRes:
        halvedResolution = [startRes for x in initialResolution]
    
    tree.changeResolution(halvedResolution)
    
    # Add the initial population
    for x in range(0, count):
        genome        = individual(dim, minGenomeVal, maxGenomeVal)
        initgenomeFitness = MICHALEWICZ_getFitness(genome)
        tree.addNode(Tree_Harry.Node(genome, fitness = initgenomeFitness))
        
    return tree


def evolve(inputTree, dim, minGenomeVal, maxGenomeVal, varienceThreshold, revisitAmount):
    """
    Evolve the population, either adding a new random individual, mutating an individual or crossingover an invididual.

    inputTree           : the population tree
    dim                 : the number of dimensions
    minGenomeVal        : the minimum value of the genome space
    maxGenomeVal        : the maximum value of the genome space   
    varienceThreshold   : the varience threshold to determine when to change resolution    
    revisitAmount       : the number of revisits needed before a resolution change can occur        

    return              : the tree object
    """ 
    
    #Account for indexing starting at 0
    maxGenomeVal    = maxGenomeVal - 1    
    
    tree            = inputTree
    affectedCell    = None
    
    #Change this to randint(0, 2) to include the random adding of individuals.
    randomNum       = randint(1, 2)
    
    if randomNum == 0:
        # randomly add other individuals to promote genetic diversity. ** NOT CURRENTLY USED **
        genome           = individual(dim, minGenomeVal, maxGenomeVal)
        rndgenomeFitness = MICHALEWICZ_getFitness(genome)
        tree.addNode(Tree_Harry.Node(genome, fitness = rndgenomeFitness))
        affectedCell = tree.searchNode(genome)        
        #print "** Rnd Node Added **"        
        
    elif randomNum == 1:
        # mutate an individual
        mutated = False       
        
        # Make sure the mutation is within one cells width but also doesnt go outside the feature space. If so, try again. 
        while mutated == False:       
            
            randomCell          = tree.getRandom()
            # Error found here 14/08/17, randomCell.cargo[0].genome was modifying the object, not just retrieving the value. list() was added to solve the problem. 
            randomNodeGenome    = list(randomCell.cargo[0].genome)
            nodesCurrentRes     = randomCell.resolution        
            
            dimToChange         = randint(0, dim - 1)
            rndChangeAmount     = randint(1, nodesCurrentRes[0])
            posOrMinus          = randint(0, 1)
            
            if posOrMinus == 0:
                if randomNodeGenome[dimToChange] - rndChangeAmount >= minGenomeVal:                    
                    randomNodeGenome[dimToChange] = randomNodeGenome[dimToChange] - rndChangeAmount                   
                    mutated = True
            
            else:
                if randomNodeGenome[dimToChange] + rndChangeAmount <= maxGenomeVal - 1:                    
                    randomNodeGenome[dimToChange] = randomNodeGenome[dimToChange] + rndChangeAmount                    
                    mutated = True
   
        genomeFitness = MICHALEWICZ_getFitness(randomNodeGenome)        
        tree.addNode(Tree_Harry.Node(randomNodeGenome, fitness = genomeFitness))                
        affectedCell = tree.searchNode(randomNodeGenome)          
        #print "** Rnd Node Mutated **"
   
    else:     
        
        # crossover parents to create children
        parent1      = tree.getRandom().cargo[0].genome
        parent2      = tree.getRandom().cargo[0].genome
        
        # Make sure they are not the same node, if so, pick a new parent. 
        while parent1 == parent2:
            parent2   =  tree.getRandom().cargo[0].genome
        
        # Create the new node at the mid-point
        parents      = [parent1, parent2]
        child        = [sum(e)/len(e) for e in zip(*parents)]
        child        = [int(i) for i in child]
        childFitness = MICHALEWICZ_getFitness(child)        
        
        tree.addNode(Tree_Harry.Node(child, fitness = childFitness))    
        affectedCell = tree.searchNode(child)
        #print "** Rnd Nodes Crossover **"
             
        
    #See if a resolution change is necessasary
    if affectedCell is not None:
        if isinstance (affectedCell.cargo[0], Tree_Harry.Node):
            if affectedCell.cargo[1] > varienceThreshold: 
                if len(affectedCell.cargo[2]) > revisitAmount:
                    tree.findAndChangeRes(affectedCell.cargo[0].genome)   
            
    return tree
 
def duplicates(inputList):
    """
    Checks for duplicate genomes in the tree. Used for error testing but kept in for a final error check at the end of a run. 
    
    inputList   : the list generated from tree.getLeavesList().
    
    return      : returns if there are no duplicates, otherwise enter the python debugger. 
    """    
    
    uniq = []
    a = inputList  
    
    for x in a:
        if x.genome not in uniq:
            uniq.append(x.genome)
        else:
            print x.genome
            pdb.set_trace()   
    return
    
def getTopVarience(inTree, noOfResults):
    """
    Returns the top n cell varience values. 
    
    inputTree   : the tree object.
    noOfResults : how many top results to return. 
    
    return      : the new list containing the top n cell varience values. 
    """  
    
    tree = inTree
    treeList = tree.getTheLeaves([])
    
    newlist = sorted(treeList, key=lambda x: x[1], reverse=True)[:noOfResults]
    newlist = [[x[0], x[1]] for x in newlist]    
    
    return newlist

    
def getTopFitness(inTree, noOfResults):
    """
    Returns the top n cell fitness values. 
    
    inputTree   : the tree object.
    noOfResults : how many top results to return. 
    
    return      : the new list containing the top n cell fitness values. 
    """  
    
    tree = inTree
    treeList = tree.getTheLeaves([])
    
    for x in treeList:
        assert isinstance(x[0], Tree_Harry.Node)
    
    newlist = []
    for x in treeList:
        theGenome = x[0].genome
        newlist.append([theGenome, MICHALEWICZ_getFitness(theGenome)])

    newlist = sorted(newlist, key=lambda x: x[1], reverse=False)[:noOfResults]  
    return newlist
    
    
def MICHALEWICZ_getFitness(x, m = 10):
    """
    Returns a fitness value based on the Michalewicz formula.  
    
    x           : the genome
    m           : changes the look of the feature space. Usually left at 10. 
    
    return      : a fitness value. 
    """  
        
    m       = m
    d       = len(x)
    theSum  = 0
    
    for i in range (0, d):
        xi      = x[i]
        new     = np.sin(xi) * (np.sin((i+1) * xi**2/np.pi))**(2*m) 
        theSum  = theSum + new

    return -(theSum)