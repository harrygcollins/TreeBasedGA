 # -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 17:36:36 2017

@author: Harry


To initialise the tree, in setup use:

    @res = the maximum size of the feature space, I've used [32,32,32] for a 3-dimension example. 
        
    yourTree =  Tree(res = [32,32,32], carg = None, child = [], parent = None, blockLocation = None, maxReso=None)
    
    yourTree.changeResolution([1,1,1])
    

To add a node:

    @param = A new instance of a Node class with the genome of the node as the param to that.         
    yourTree.addNode(Node([25,25,25]))        
    

@Return Node Class or None
@Param genome to search.
    
    yourTree.searchNode([5,5,5])

"""

import numpy as np
import random as rnd
import pdb

def getBlockLocation(listInput, resolution):
    assert (len(listInput) == len(resolution))    
    block = []
    for index in range(len(listInput)):
        if listInput[index] == 0 or resolution[index] ==0:
            block.append(int((int((listInput[index]) / 1))))
        else:
            block.append(int((int((listInput[index]) / (resolution[index])))))
    
    for bb in block:
        if bb > 1:
            pdb.set_trace()            
    return block

def getInBlockCoord(listInput, resolution):
    assert (len(listInput) == len(resolution))
    coord = []
    for index in range(len(listInput)):
        if int(listInput[index]) == 0 or int(resolution[index]) ==0:
            coord.append(int((int((listInput[index]) % 1))))
        else:
            coord.append(int((int((listInput[index]) % (resolution[index])))))
    return coord            
    
def getIntFromBinary(inList):
    binaryString = ""
    for i in inList:
        binaryString = binaryString + str(i)
    return int(binaryString, 2)

class Node(object):
    def __init__(self, genome, fitness = 100, blockLocation = None, remainder = None):
        self.genome = genome
        self.blockLocation = blockLocation
        self.remainder = remainder
        self.fitness = fitness
        
    def __repr__(self):
        return str(self.genome)
    
  
def getVarience(inputList):
    tempList = []
    
    for x in inputList:
        tempList.append(x.fitness)
    
    # ddof = 1 gives the unbiased varience.
    return float(np.var(tempList, ddof=1))
  
class Tree(object):
                                 
    def __init__(self, res, carg = None, child = [], parent = None, blockLocation = None, maxReso = None):

        """
        Cargo:
                [0]     =   Node Object.
                [1]     =   Cell Varience (float).
                [2]     =   List of all nodes that have been in that cell.                
        """        
           
        self.resolution = res
        self.cargo = carg 
        self.children = child
        self.parent = parent
        
        if blockLocation is None:
            self.blockLocation = [0 for x in self.resolution]
        else:
            self.blockLocation = blockLocation        
        
        if maxReso is None:
            self.maxResolution = res
        else:
            self.maxResolution = maxReso
    
    def __repr__(self):
        return str(self.cargo[0])
    
    
    # returns a random node from the Tree.
    def getRandom(self):
        if self.children == []:
            return self            
        else:
            return rnd.choice(self.children).getRandom()

    
    def addNode(self, newNode):     
        
        assert isinstance(newNode, Node)        
        
        
        if newNode.remainder is None:            
            newNode.remainder = getInBlockCoord(newNode.genome, self.resolution)

        newNode.blockLocation = getBlockLocation(newNode.remainder, self.resolution) 
        newNode.remainder = getInBlockCoord(newNode.remainder, self.resolution)
        
        #First Insertion
        if(self.cargo[0] == None): 
            #print "Tree:    First Insertion"

            self.cargo[0] = newNode
            self.cargo[1] = 0.0
            self.cargo[2] = [newNode]
            self.blockLocation = newNode.blockLocation
            assert isinstance(self.cargo[0], Node)
            self.moveNodesDown()
            
            return

        #Root Case, first covering first insertion and then when the resolution isnt changed.        
        elif((self.parent is None) and (getBlockLocation(newNode.genome, self.resolution) == self.blockLocation)):

            
            if (not isinstance(self.cargo[0], str)):
                if(self.cargo[0].blockLocation == newNode.blockLocation):
                    if(newNode.fitness > self.cargo[0].fitness):
                        #print "Tree:    Node Added"
                        self.cargo[0] = newNode
                        self.cargo[2].append(newNode)
                        self.cargo[1] = getVarience(self.cargo[2])
                        self.blockLocation = newNode.blockLocation
                        self.moveNodesDown()
                        assert isinstance(self.cargo[0], Node)
                        return
                    else:
                        #print "Tree:    Revisit, fitness lower"
                        self.cargo[2].append(newNode)
                        self.cargo[1] = getVarience(self.cargo[2])
                        assert isinstance(self.cargo[0], Node)
                        return
            
        #Parent Case
        if(self.cargo[0] == "Parent"):
            
            tempRes = [x / 2 for x in self.resolution]
            
            newNode.blockLocation   = getBlockLocation(newNode.remainder, tempRes)  
            newNode.remainder       = getInBlockCoord(newNode.remainder, tempRes)          
            
            #Search if there is already something in a cell it should be in
            for c in self.children:

                if(c.blockLocation == newNode.blockLocation):
                    if(c.cargo[0] == "Parent"):
                        c.addNode(newNode)
                        return 
                        
                    else:
                        #REVISIT, keep the fittest individual and update varience. 
                        if newNode.genome != c.cargo[0].genome:
                            
                            if(newNode.fitness > c.cargo[0].fitness):
                                #print "Tree:    Revisit, fitness higher"
                                c.cargo[0] = newNode
                                c.cargo[2].append(newNode)
                                c.cargo[1] = getVarience(c.cargo[2])
                                c.moveNodesDown()
                            else:
                                #print "Tree:    Revisit, fitness lower"
                                c.cargo[2].append(newNode)
                                c.cargo[1] = getVarience(c.cargo[2])
                            assert isinstance(c.cargo[0], Node)
                            return
                            
                        else:
                            # Genome already in the tree, discard.
                            #print " ########################################################## "
                            #print str(c.cargo[0].genome) + " :: " + str(newNode.genome)
                            #print " ########################################################## "
                            return
            
            #If there are no children at all, it is the first and should be added. 
            #(new very fist child of this branch)
            if (self.children == []):
                #print "Tree:    Node Added"
                self.children.append(Tree(res = tempRes, 
                                      carg = [newNode,0.0,[newNode]],
                                      parent = self, 
                                      child = [], 
                                      blockLocation = newNode.blockLocation, 
                                      maxReso = self.maxResolution))
                
                self.children[-1].cargo[0].blockLocation = self.children[-1].blockLocation                         
                self.children[-1].moveNodesDown()      
                #assert isinstance(self.children[-1].cargo[0], Node)
                return
                
            
            #If the parent already has children, but not in the block you are in.    
            if self.children != []: 
                
                #If children[0] is a parent, not at correct resolution
                if(self.children[0].cargo[0] is "Parent"):
                    self.children.append(Tree(res = tempRes, 
                                              carg = ["Parent",0.0, []], 
                                              parent = self, 
                                              child = [], 
                                              blockLocation = newNode.blockLocation, 
                                              maxReso = self.maxResolution))                      
                    
                    self.children[-1].addNode(newNode)
                    return
                
                else:
                    #The resolution is correct, it is in a new block so just add to the end of the children list. (new cell)
                    #print "Tree:    New Node Added"
                    self.children.append(Tree(res = tempRes, 
                                              carg = [newNode, 0.0, [newNode]], 
                                              parent = self, 
                                              child = [], 
                                              blockLocation = newNode.blockLocation, 
                                              maxReso = self.maxResolution))
                                              
                    self.children[-1].cargo[0].blockLocation = self.children[-1].blockLocation 
                    self.children[-1].moveNodesDown()  
                    assert isinstance(self.children[-1].cargo[0], Node)          
                    return
        

                            
    def moveNodesDown(self):
        
        if isinstance(self.cargo[0], Node) and self.resolution > self.maxResolution and self.resolution[0] / 2 >= 1:

            tempRes = [x / 2 for x in self.resolution]
            
            self.cargo[0].blockLocation = getBlockLocation(self.cargo[0].remainder, tempRes)              
            self.cargo[0].remainder = getInBlockCoord(self.cargo[0].remainder, tempRes)         
            
            self.children.append(Tree(res = tempRes, 
                                      carg = [self.cargo[0], 0.0, [self.cargo[0]]], 
                                      parent = self, 
                                      child = [], 
                                      blockLocation = self.cargo[0].blockLocation,
                                      maxReso = self.maxResolution))
            self.cargo[0] = "Parent"
            self.cargo[1] = 0.0
            self.cargo[2] = []
            assert isinstance(self.children[-1].cargo[0], Node)
            self.children[-1].moveNodesDown()
        else: 

            for b in self.children:
                b.moveNodesDown()
            return 
            
    def searchNode (self, coord):
        if len(self.children) == 0:
            if self.cargo[0].genome == coord:
                return self
        else:
            branchSum = []
            for b in self.children:
                ret = b.searchNode(coord)
                if ret != None:
                    branchSum.append(ret)
                    return branchSum[0]
                    
    def findAndChangeRes(self, coord):
        if len(self.children) == 0:
            if self.cargo[0].genome == coord:
                newRes = [x/2 for x in self.resolution]
                self.changeResolution(newRes)     
                return
        else:                              
            for b in self.children:
                b.findAndChangeRes(coord)  
            
            return
        
  
    def changeResolution(self, newResolution):  
        
        # ---- Deals with the Root Node --- #  
        
        # Check root for a greater resolution.
        if self.resolution > newResolution:
            self.maxResolution = newResolution
            
            # If no values have been added yet. 
            if self.cargo[0] is None:               
                self.maxResolution = newResolution
            
            # If there is a single node at the root.
            elif isinstance(self.cargo[0], Node):
                self.moveNodesDown()

        # --- Deals with the children if the node is a Parent --- #
            
            # If the root is a parent, look at the resolution of the children.
            elif self.cargo[0] is "Parent":
                
                for child in self.children:
                   if child.cargo[0] is "Parent":
                       child.changeResolution(newResolution)
                           
                   elif child.resolution > newResolution:
                       tempRes = [x / 2 for x in child.resolution]
        
                       child.cargo[0].blockLocation = getBlockLocation(child.cargo[0].remainder, tempRes)
                       child.cargo[0].remainder = getInBlockCoord(child.cargo[0].remainder, tempRes)   
                       
                       child.children.append(Tree(res = tempRes, 
                                                  carg = [child.cargo[0], 0.0, [child.cargo[0]]], 
                                                  parent = child, 
                                                  child = [], 
                                                  blockLocation = child.cargo[0].blockLocation,  
                                                  maxReso = child.maxResolution))
                       
                       child.cargo[0] = "Parent"
                       child.cargo[1] = 0.0
                       child.cargo[2] = []
                   
                       if child.children[-1].resolution > newResolution:
                           child.changeResolution(newResolution)

            # If the cargo does not match one of these, throw an exception.
            else:
                raise ValueError('The Resolution could not change, error somewhere around here')
            
            # If it gets this far, Resolution has changed. 
            return

        # Resolution already finer than the new resolution, Resolution Change not needed.   
        else:
            return

    
    def getTheLeaves (self, theList):
        #returns the leaves of a given branch. For leaves of the tree, specify data
        if len(self.children) == 0:
            assert isinstance(self.cargo[0], Node)
            theList.append([self.cargo[0], self.cargo[1]])
            last = theList[-1]
            assert isinstance(last[0], Node)
            return theList
        else:
            theList = theList
            for b in self.children:
                theList = b.getTheLeaves(theList)
                last = theList[-1]
                assert isinstance(last[0], Node)
            return theList

    def getLeavesList (self, theList): #must start with an empty list as a parameter.
        if len(self.children) == 0 and isinstance(self.cargo[0],Node) and self.cargo[0].genome != None:                
            return self.cargo[0]
        
        else:
            for b in self.children:
                x = b.getLeavesList(theList)
                if isinstance(x, Node):
                    theList.append(x)

            return theList