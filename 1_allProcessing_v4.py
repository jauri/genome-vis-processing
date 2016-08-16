import json
import sys
import Queue
from sets import Set

#------------ Functions to create initial data types
def createBlock(idx, score, numberAligned, segList, maxLength):
    return {'id': idx, 'score': score, 'numberAligned': numberAligned, 'segList': segList, "maxLength": maxLength}

#--Returns species id, unless it hasn't been encountered yet
def getSpeciesId(speciesName):
    if speciesName in speciesToIdHM:
        return speciesToIdHM[speciesName]
    else:
        id = len(speciesToIdHM) 
        species.append( {"id": id, "name": speciesName, "contigList": [] , "maxNT": 0} )
        speciesToIdHM[speciesName] = id
        return id

#--Returns contig id, unless it hasn't been encountered yet 
def getContigId(contigName, speciesId):
    if contigName in contigsToIdHM:
        return contigsToIdHM[contigName]
    else:
        id = len(contigsToIdHM) 
        contigs.append( {"id": id, "name": contigName, "segmentList": [], "orderedBlocks": [] , "contigLength": 0, "speciesId": speciesId} )
        contigsToIdHM[contigName] = id
        return id

#-- Given a segment read from the file, add to the contig in the appropriate order
#----- so we have an ordered list of segments and blocks from alignment file
def addSegmentToContig(segment, contigId):
    
    segId = segment["id"]
    segStart = segment["start"]
    segBlock = segment["block"]
    segLength = segment["size"] 

    insertSegmentAtIdx = 0
    segArray = contigs[contigId]["segmentList"]
    orderedBlocks = contigs[contigId]["orderedBlocks"] 

    found = False 
    for i in range( 0, len(segArray) ):
        if segStart < segments[ segArray[i] ]["start"]:
            insertSegmentAtIdx = i 
            found = True
            break
            
    if( not found ):
        segArray.append(segId) 
        orderedBlocks.append(segBlock)
    else:
        segArray.insert(insertSegmentAtIdx, segId ) 
        orderedBlocks.insert( insertSegmentAtIdx, segBlock )

    contigs[contigId]["segmentList"] = segArray #reassign this list back to the contig 
    contigs[contigId]["orderedBlocks"] = orderedBlocks

#############################################################################
###########----------------------BEGIN----------------------------###########
#############################################################################

#------input parameters
runName = "ecoli"
mafFileName = "inputData/EcoliMugsyOut2.maf"
f = open( mafFileName)

#------data strutures to output
blocks = []  #stores the aligned blocks of segments
segments = [] #just the segment info for each genome
contigs = [] #which contigs are present for the strains
species = [] #which species are represented 
allJson = {} #this is where all data will go

#------- to help with lookups
speciesToIdHM = {}
contigsToIdHM = {} 

#------loop through file, this is the id assigned to blocks and segments 
blockNumber = 0
segmentNumber = 0

#------ Create 0th block, which is empty and represents the start
blocks.append( createBlock(0, 0, 0, {}, 0) )


#-------READ FILE-----------------------------------------------------------------------------

#look throug maf file
# format : some line start with a, other lines start with s
# a lines are the aligned blocks
# s lines are the segments in the block, they occur after an a line
for line in f: 

    lineType = line[:1] # is it an 'a' line or an 's' line? 

    if lineType == "a":

        if( blockNumber > 0 ): #add blocks after all the segments within have been read 
            blocks.append( createBlock(blockNumber, score, numberAligned, segList, maxLength) )

        #increment count.  Note, there is no 0th block from data, but we invented one as a source 
        blockNumber = blockNumber+1 
        # print blockNumber

        #here is the information we need for each block 
        segList = []; #this will store the ids of the segments in the block
        maxLength = 0 #what is the length of the longest segment in this block

        #block info 
        tokens = line.split(" ") #split the line
        score = int(tokens[1].split("=")[1]) #score value for block
        numberAligned = int(tokens[3].split("=")[1]) #how many segments are in the block
        

    #these lines always follow an 'a' line 
    if lineType == "s":
        
        tokens = line.split()

        speciesName = int(tokens[1].split('.')[0]) #which species is this from (maf file is missing part of species name)
        contig = tokens[1].split('.')[1] #which contig is it from?  
        speciesId = getSpeciesId(speciesName) #assign a unique id to the species, or retrieve from hm
        contigId = getContigId(contig, speciesId) #assign a unique id to the contig, or retrieve from hm

        # addToSpecies(speciesId, speciesName, contigId, segmentNumber);
        # addToContigs(contigId, contig, speciesId, segmentNumber);

        start = int(tokens[2]) #where does this segment start from within its contig?
        size = int(tokens[3]) #how long is this segment (not counting gaps )

        #is it reverse compliment? 
        reverseComplimentStr = tokens[4] 
        revComp = False
        if reverseComplimentStr.find("+") == -1:
            revComp = True
        featureLen = int(tokens[5]) #if so, need to know the length of this contig

        if( revComp ):
            start = featureLen - start - size;  #reverse compliment coordinates start at end of contig

        if( size > maxLength ): 
            maxLength = size #keep track of which segment is longest 

        segList.append( {"id": segmentNumber} )  #store the id of the segment 

        #store info about the segment 
        segments.append( {"id": segmentNumber, "speciesId": speciesId, "contigId":contigId, "start": start, "size": size, "reverseCompliment": revComp, "block": (blockNumber)} );

        #handle species and contig info
        addSegmentToContig( segments[segmentNumber], contigId )
        contigs[contigId]["contigLength"] = featureLen;

        segmentNumber = segmentNumber+1


#add last one not covered in loop 
blocks.append( createBlock(blockNumber, score, numberAligned, segList, maxLength) )

#and an ending block
blockNumber = blockNumber+1 
blocks.append( createBlock((blockNumber), 0, len(species), {}, 0) )  

#fix the starting block  to have the right number of species
blocks[0]["numberAligned"] = len(species) 

#-------Resolve species and contig relationships

#add contig to correct species, with longest contig first
for i in range(0, len(contigs) ):
    speciesId = contigs[i]["speciesId"]

    insertIdx = 0
    contigList = species[speciesId]["contigList"] 

    contigLen = contigs[i]["contigLength"]

    #put it in the right order (longest first)
    for j in range(0, len( contigList) ):
        otherLen = contigs[ contigList[j] ]["contigLength"]
        if contigLen > otherLen: 
            insertIdx = j
            break

    contigList.insert( insertIdx, contigs[i]["id"] )
    species[speciesId]["contigList"] = contigList 
    species[speciesId]["maxNT"] = species[speciesId]["maxNT"] + contigs[i]["contigLength"] 


#--------- Order blocks of aligned segments for each species 
# This gives each species preferred order for the aligned blocks 

for i in range(0, len(species) ):
    contigList = species[i]["contigList"]
    orderedBlocks = [0]  
    for j in reversed(contigList): #range(0, len( contigList) ):  
        orderedBlocks = orderedBlocks + contigs[ j ]["orderedBlocks"]  

    orderedBlocks.append( blockNumber )
    species[i]["orderedBlocks"] = orderedBlocks

##DEBUGGING
# print "ordered blocks: "
# for i in range(0, len(species) ):
#     print "species: %d" % i
#     print species[i]["orderedBlocks"]



# done reading the files! 


##############################################################
##############################################################
###----------------- ORDERING AND POSITIONING -------------###
##############################################################
##############################################################

genomesAndAlignments = {}
maxLengths = [] 

startingPositions = []  #the starting positions for the blocks
blockSize = [] #the size of each block
blockDepth = [] # not using this so much right now, but the 'depth' from start of the block
overlapsWith = [] #which blocks overlap with each block, not using so much yet
globalSegmentDepth = [] # not using right now

lengthOfCoordinateSystem = 0 #given all the blocks and their lengths, how long is the total set of blocks in nt?

##------INITIALIZE -------
#populate the array storing blocks and strains
for i in range(0, len(species) ):
    genomesAndAlignments[ species[i]["id"] ] = species[i]["orderedBlocks"] 
    for j in range(0, len( genomesAndAlignments[ species[i]["id"] ] ) ) :
        genomesAndAlignments[ species[i]["id"] ][j] = int( genomesAndAlignments[ species[i]["id"] ][j] )

numBlocks = len(blocks)
numGenomes = len(species)

#populate the max length array- the length of the longest segment in each block
for i in range(0, len(blocks) ):
    maxLengths.append( int(blocks[i]["maxLength"]) )
    startingPositions.append(0)
    blockSize.append(0)
    blockDepth.append(0)
    overlapsWith.append([])
    globalSegmentDepth.append(0)


###-------- CREATE GRAPH -------- 

#here is where I store the initial ordered sequence as a graph in adjacency list format
# note, there will be cycles in this graph 
adjLists = []
for i in range(0, numBlocks):
    arr = []
    adjLists.append(arr)

allEdges = [] #this is used to get the sorted list of edges by weight for the greedy algorithm 


# add edge and sum up weight
def addEdge(g, f, t):
    found = False 
    for i in range( len(g[f]) ):
        if( g[f][i]["to"] == t ):
             g[f][i]["w"] = g[f][i]["w"]+1
             found = True
    if not found: 
        g[f].append({"to": t, "w": 1, "from": f})

#add edge and overwrite weight for a new weight
def addEdgeWithWeight(g, f, t, w, c):
    found = False 
    for i in range( len(g[f]) ):
        if( g[f][i]["to"] == t ):
             g[f][i]["w"] = w
             found = True
    if not found: 
        g[f].append({"to": t, "w": w, "from": f, "cost":c})

#remove a specific edge (eg when it is in a cycle)
def removeEdge( g, f, t):
    for i in range( len(g[f]) ):
        if( g[f][i]["to"] == t ):
            g[f].pop(i)


root = 0 # the 0th block is the root 

#populate the adjacecy list with initial data
for key in genomesAndAlignments:
    arr = genomesAndAlignments[key]

    for i in range(len(arr)-1):
        addEdge( adjLists, arr[i], arr[i+1] )

##DEBUGGING 
# print len(adjLists)

#populate the all edges dataset and the block size (really the number of edges into and out of node)
for i in range(0, len(adjLists) ):
    for j in range( 0, len (adjLists[i]) ):
        allEdges.append( { "from": i, "to": adjLists[i][j]["to"], "w": adjLists[i][j]["w"] , "cost": 0})
        blockSize[ i ] = blockSize[ i ] + adjLists[i][j]["w"]


#note, if this works can improve performance by doing lookup on segs by species id 
#and have null value where no species segment is present
#compute the costs of the edges
#    if an edge connects blocks where most segments are 'far away', adds to the 'cost'
for i in range(0, len( allEdges) ):
    e = allEdges[i]

    parentSegments = blocks[ e["from"] ]["segList"]
    childSegments = blocks[ e["to"] ]["segList"]
    cost = 0

    for j in range(0, len(parentSegments)):
        pSeg = parentSegments[j]["id"]
        pSpecies = segments[pSeg]["speciesId"]

        for k in range(0, len(childSegments)):
            cSeg = childSegments[k]["id"]
            cSpecies = segments[cSeg]["speciesId"]
            if( cSpecies == pSpecies ):
                cost = cost + abs( segments[cSeg]["start"] - (segments[pSeg]["start"]+segments[pSeg]["size"])  )#child start minus parent end
    e["cost"] = cost             


# print "ADDED COSTS "
# print allEdges


### DEBUGGING
#print the starting point
# print "------- inputs"
# print adjLists
# print allEdges


#this is used in the sorting for the greedy algorithm
def compareWeights(item1, item2):
    if item1["w"] > item2["w"]:
        return -1
    elif item1["w"] < item2["w"]:
        return 1
    else:
        return 0

def compareWeightsAndCosts(item1, item2):
    if item1["w"] > item2["w"]:
        return -1
    elif item1["w"] < item2["w"]:
        return 1
    else:
        if item1["cost"] < item2["cost"]:
            return -1
        elif item1["cost"] > item2["cost"]:
            return 1
        else:
            return 0
    return 0 #should never get here... 

def compareCosts(item1, item2):
    if item1["w"] < item2["w"]:
        return -1
    elif item1["w"] > item2["w"]:
        return 1
    else:
        return 0


#get the sorted array, so we can greedily add edges by weight
sortedArr = sorted(allEdges, cmp=compareWeightsAndCosts)

#DEBUGGING
# print "------- sorted"
# print sortedArr

#recursive helper function
def cyclicUtil(i, visited, recStack, g):
    # print "idx= %d" % i 
    # print "g.len= %d" % len(g)
    if( not visited[i] ):
        visited[i] = True
        recStack[i] = True

        arr = g[i]
        for j in range(0, len(arr)):
            # print "       j= %d" % j 
            # print "       arr.len= %d" % len(arr)
            if( not visited[ arr[j]["to"] ] and cyclicUtil( arr[j]["to"], visited, recStack, g) ):
                return True
            elif( recStack[ arr[j]["to"] ]):
                return True

    recStack[i] = False
    return False

#recursive main function to check for cycles
def isCyclic(g):
    visited = []
    recStack = []
    for i in range(0, len(g)):
        visited.append(False)
        recStack.append(False)

    for i in range(0, len(g)):
        if( cyclicUtil(i, visited, recStack, g) ):
            return True

    return False



#now build up the graph using greedy approach selected edges with most weight
#but don't add if it would create a cycle

#initialize
maxNoCycle = [] # holds the adjacency list of edges with no cycles 
for i in range(0, numBlocks):
    arr = []
    maxNoCycle.append(arr)

#loop through edges in order of weights
possibleOrphanList = [] # holds edges removed in order to create cycle-less graph

##--------  GREEDY GRAPH NO CYCLE GENERATION
#loop through edges in order of weights
for i in range( len(sortedArr) ):
    e = sortedArr[i]

    addEdgeWithWeight(maxNoCycle, e["from"], e["to"], e["w"], e["cost"])

    # print "ADDED %d, %d" % (e["from"], e["to"] )

    if( isCyclic(maxNoCycle) ):
        removeEdge(maxNoCycle, e["from"], e["to"] ) 
        # print "REMOVED f: %d, t: %d, w: %d" % (e["from"], e["to"], e["w"])
        possibleOrphanList.append( {"to": e["to"], "from": e["from"], "w": e["w"], "cost": e["cost"]} ) 
    #     for i in range(0, e["w"]):
    #         print "%d %d" % (1000, e["to"] )
    # else:
    #     for i in range(0, e["w"]):
    #         print "%d %d" % (e["from"], e["to"] )

# print "%d %d" %( 0, 1000)

### Are any 'orphans' created from this cycle-less graph
## orphans are blocks with no parents
actualOrphans = Set() 
allRemoved = Set() 
# blackList = [False]*len(blocks) 
for i in range(0, len(possibleOrphanList)):
    foundParent = False
    parent = 0
    allRemoved.add( possibleOrphanList[i]["to"] )#removing duplicates ... 

    for j in range(0, len(maxNoCycle) ):
        arr = maxNoCycle[j]  #children

        for k in range(0, len(arr)):
            e = arr[k]

            if e["to"] == possibleOrphanList[i]["to"] :
                foundParent = True

    if not foundParent:
        actualOrphans.add(possibleOrphanList[i]["to"] )
        # blackList[ possibleOrphanList[i]["to"] ] = True










## DEBUGGING
# print "------- no cycle max weights"
# print maxNoCycle

## DEBUGGING 
# print "-------- orphans "
# print actualOrphans

#add orphans to list

#### DOESN"T WORK
orphansFinal = [] 
while len(actualOrphans) > 0:
    orphansFinal.append( actualOrphans.pop() )


for i in range(0, len(orphansFinal)):
    #get orphan's first child
    #childId = maxNoCycle[orphansFinal[i]][0]["to"] #children
    foundParent = False
    newParent = 0
    edges = maxNoCycle[orphansFinal[i]]
    for j in range(0, len(edges)): 
        e = edges[j] #this is the orphan's first edge out

        #find another parent for the selected child
        for k in range(0, len(maxNoCycle)):
            children = maxNoCycle[k]
            for l in range(0, len(children) ):
                childId = children[l]["to"]

                if( childId == e["to"] and children[l]["from"] != e["from"] ):
                    newParent = children[l]["from"]
                    foundParent = True
                    break

    # if( foundParent ):
    #     print "%d %d" % ( newParent, orphansFinal[i] )

## create black list from orphans
# queue = Queue.Queue() 
# for i in range(0, len(orphansFinal) ):
#     queue.put(orphansFinal[i])

# visited = [0]*len(blocks)

# while not queue.empty():
#     parent = queue.get() 
#     visited[parent] = True 

#     children = maxNoCycle[parent] #children
#     for i in range(0, len(children) ):
#         childId = children[i]["to"] 
#         blackList[childId] = True
#         if( not visited[ childId ] ):
#             queue.put(children[i]["to"]) 

# print blackList 


##########-------FIRST PASS POSITIONING 
root = 0
queue = Queue.Queue() 
queue.put(root)

visited = [0]*len(blocks)
positioned = [0]*len(blocks)
positioningCost = [0]*len(blocks)


while not queue.empty():
    parent = queue.get() 
    visited[parent] = visited[parent]+1

    children = maxNoCycle[parent] #children
    l = maxLengths[parent] #max length of this segment
    startingPosition = startingPositions[parent]

    minDistToChild = sys.maxint

    sortedChildren = sorted( children, cmp=compareCosts)

    runningEnd =startingPosition + l 
    for j in range(0, len(sortedChildren)):
        childId = sortedChildren[j]["to"] 
        childStart = startingPositions[ childId ]

        if positioned[childId] == 1:
            if( positioningCost[childId] > sortedChildren[j]["cost"] ): #the cost of the former placement is greater
                startingPositions[childId] = runningEnd
                positioningCost[childId] = sortedChildren[j]["cost"]
                queue.put(childId)
                positioned[childId] = 1
            elif( sortedChildren[j]["cost"] == 0):
                 if( childStart < runningEnd ):
                    startingPositions[ childId ] = runningEnd  #shift child over
                    queue.put( children[j]["to"] )
                    positioned[childId] = 1
                    positioningCost[childId] = 0
                    # print "shove over child: %d pushed by: %d from: %d to: %d  visited: %d" %( childId, b, childStart, startingPositions[ children[j]["to"] ], visited[childId] )

        else:
            startingPositions[childId] = runningEnd
            positioningCost[childId] = sortedChildren[j]["cost"]
            queue.put(childId)
            positioned[childId] = 1



    # if( len(sortedChildren) > 0):
    #     childId = sortedChildren[j]["to"] 
    #     childStart = startingPositions[ childId ]

    # for j in range(0, len(sortedChildren)):
    #     childId = sortedChildren[j]["to"] 
    #     childStart = startingPositions[ childId ]



    # for j in range(0, len(children) ):
    #     childId = children[j]["to"] 
    #     childStart = startingPositions[ childId ]

    #     if( positioned[childId] ):
    #         print "already positioned"
    #     else:


    #     #if child's starting position should change, add to queuw
    #     if( childStart < startingPosition + l ):
    #         startingPositions[ childId ] = startingPosition + l   #shift child over
    #         queue.put( children[j]["to"] )
    #         # print "shove over child: %d pushed by: %d from: %d to: %d  visited: %d" %( childId, b, childStart, startingPositions[ children[j]["to"] ], visited[childId] )
    #     elif( visited[childId] == 0 ):
    #         queue.put( children[j]["to"] )

        # distFromEnd = startingPositions[childId] - (startingPositions[b]+l) 
        # if(distFromEnd< minDistToChild):
        #     minDistToChild = distFromEnd
          
        #check to make sure parent is in range of children
    # if( minDistToChild > 0 ):
    #     startingPositions[b] = startingPositions[childId] - l
        # queue.put(b)

####CHECK FOR ORPHANS ...

# #find all the nodes which start at 0, and add them to a new 'root queue' and position behind child
# rootList = []
# for i in range(0, len(maxNoCycle) ):
#     arr = maxNoCycle[i] #children
#     l = maxLengths[i] #max length of this segment
#     if( startingPositions[ i ] == 0 ):
#         rootList.append( i )
#         earliestChildPosition = 9999999999999999  
            
#         for j in range(0, len(arr)):
#             e = arr[j]
#             if( startingPositions[ e["to"] ] <  earliestChildPosition ) :
#                 earliestChildPosition = startingPositions[ e["to"] ]

#         ##OLD METHOD
#         # if earliestChildPosition != 9999999999999999:
#         #     moveTo = earliestChildPosition - maxLengths[i ]
#         #     if( moveTo < 0 ):
#         #         startingPositions[i] = 0  #don't put it at negative value, always start the parent at 0
#         #     else:
#         #         startingPositions[ i ] = moveTo 
#         #     # print "MOVED i = %d and pos = %d with max length = %d" % (i, startingPositions[ i ], maxLengths[i] )

#         ##NEW METHOD
#         if earliestChildPosition != 9999999999999999:
#             moveTo = earliestChildPosition #- maxLengths[i ]
#             if( moveTo < 0 ):
#                 startingPositions[i] = 0  #don't put it at negative value, always start the parent at 0
#             else:
#                 startingPositions[ i ] = moveTo 
#             print "MOVED i = %d and pos = %d with max length = %d" % (i, startingPositions[ i ], maxLengths[i] )

#     # print "visited: ---- "
#     # print visited
# print "first pass outcome: "
# for i in range(0, len(blocks)):
#     print "block %d, starts at %d" %( i, startingPositions[i] )


### RESOLVE ORPHANS- place next to their earliest child
# rootList = []
# for i in range(0, len(orphansFinal)):
#     rootList.append(orphansFinal[i])
#     arr = maxNoCycle[ orphansFinal[i] ] #get children
#     earliestChildPosition = 9999999999999999

#     for j in range(0, len(arr)):
#         if ( startingPositions[ e["to"] ] <  earliestChildPosition ) :
#             earliestChildPosition = startingPositions[ e["to"] ]#put it right on top of its child 

#     if earliestChildPosition != 9999999999999999:
#         moveTo = earliestChildPosition #- maxLengths[i ]
#         if( moveTo < 0 ):
#             startingPositions[i] = 0  #don't put it at negative value, always start the parent at 0
#         else:
#             startingPositions[ i ] = moveTo 
#         print "MOVED i = %d and pos = %d with max length = %d" % (i, startingPositions[ i ], maxLengths[i] )


##### RESOLVE ORPHANS- place after their parents  JUST ACTUAL ORPHANS 
# rootList = []
# for i in range(0, len(orphansFinal)):
#     rootList.append(orphansFinal[i])

#     parentList = []
#     parentPosition = 0
#     for j in range(0, len( possibleOrphanList) ): 
#         if possibleOrphanList[j]["to"] == orphansFinal[i]:
#         #     parentList.append(possibleOrphanList[j]["from"] )
#             parentEnd = startingPositions[ possibleOrphanList[j]["from"] ] + maxLengths[ possibleOrphanList[j]["from"] ]
#             if parentEnd > parentPosition:
#                 parentPosition =  parentEnd

#     startingPositions[orphansFinal[i]] = parentEnd


rootList = []
rootList.append(0)
# allRemovedNextPass = Set() 
while len(allRemoved) > 0:
    block = allRemoved.pop() 

    # if( block == 94 )
    #     print "block 94 is at %d " % startingPositions[block]

    if( startingPositions[block] == 0 ):
        rootList.append(block)

        parentList = []
        parentPosition = 0#sys.maxint
        minCost = sys.maxint
        for j in range(0, len( possibleOrphanList) ): 
            if possibleOrphanList[j]["to"] == block:
                if possibleOrphanList[j]["cost"] < minCost:
                    parentPosition = startingPositions[ possibleOrphanList[j]["from"] ] + maxLengths[ possibleOrphanList[j]["from"] ]
                    minCost = possibleOrphanList[j]["cost"]
                #parentList.append(possibleOrphanList[j]["from"] )
                # parentEnd = startingPositions[ possibleOrphanList[j]["from"] ] + maxLengths[ possibleOrphanList[j]["from"] ]
                # if parentEnd > parentPosition:
                #     parentPosition =  parentEnd

        startingPositions[block] = parentPosition
        # print "moving %d to behind parent %d" 



 # for i in range(0, len( possibleOrphanList) ): 
 #    if( startingPositions[ possibleOrphanList[j]["to"] ]):
 #        parentEnd = startingPositions[ possibleOrphanList[j]["from"] ] + maxLengths[ possibleOrphanList[j]["from"] ]
 #        if parentEnd > parentPosition:
 #            parentPosition =  parentEnd


    # arr = maxNoCycle[ orphansFinal[i] ] #get children
    # earliestChildPosition = 9999999999999999

    # for j in range(0, len(arr)):
    #     if ( startingPositions[ e["to"] ] <  earliestChildPosition ) :
    #         earliestChildPosition = startingPositions[ e["to"] ]#put it right on top of its child 

    # if earliestChildPosition != 9999999999999999:
    #     moveTo = earliestChildPosition #- maxLengths[i ]
    #     if( moveTo < 0 ):
    #         startingPositions[i] = 0  #don't put it at negative value, always start the parent at 0
    #     else:
    #         startingPositions[ i ] = moveTo 
    #     print "MOVED i = %d and pos = %d with max length = %d" % (i, startingPositions[ i ], maxLengths[i] )

# print "orphan resolution: "
# for i in range(0, len(blocks)):
#     print "block %d, starts at %d" %( i, startingPositions[i] )

## SECOND METHOD TO RESOLVE ORPHANS 
# for i in range(0, len(rootlist)):
#     orphansFinal


# ###### SECOND PASS 
queue = Queue.Queue() 
for i in range(0, len(rootList)): #add multiple roots this time 
    queue.put(rootList[i])

visited = [0]*len(blocks)
positioned = [0]*len(blocks)
positioningCost = [0]*len(blocks)


while not queue.empty():
    parent = queue.get() 
    visited[parent] = visited[parent]+1

    children = maxNoCycle[parent] #children
    l = maxLengths[parent] #max length of this segment
    startingPosition = startingPositions[parent]

    minDistToChild = sys.maxint

    sortedChildren = sorted( children, cmp=compareCosts)

    runningEnd =startingPosition + l 
    for j in range(0, len(sortedChildren)):
        childId = sortedChildren[j]["to"] 
        childStart = startingPositions[ childId ]

        if positioned[childId] == 1:
            if( positioningCost[childId] > sortedChildren[j]["cost"] ): #the cost of the former placement is greater
                startingPositions[childId] = runningEnd
                positioningCost[childId] = sortedChildren[j]["cost"]
                queue.put(childId)
                positioned[childId] = 1
            elif sortedChildren[j]["cost"] == 0:
                 if( childStart < startingPosition + l ):
                    startingPositions[ childId ] = startingPosition + l   #shift child over
                    queue.put( children[j]["to"] )
                    positioned[childId] = 1
                    positioningCost[childId] = 0
                    # print "shove over child: %d pushed by: %d from: %d to: %d  visited: %d" %( childId, b, childStart, startingPositions[ children[j]["to"] ], visited[childId] )

        else:
            startingPositions[childId] = runningEnd
            positioningCost[childId] = sortedChildren[j]["cost"]
            queue.put(childId)
            positioned[childId] = 1



# queue = Queue.Queue() 
# for i in range(0, len(rootList)): #add multiple roots this time 
#     queue.put(rootList[i])

# visited = [0]*len(blocks)

# while not queue.empty():
#     b = queue.get() 
#     visited[b] = visited[b]+1

#     children = maxNoCycle[b] #children
#     l = maxLengths[b] #max length of this segment
#     startingPosition = startingPositions[b]

#     minDistToChild = 9999999999999999

#     for j in range(0, len(children) ):
#         childId = children[j]["to"] 
#         childStart = startingPositions[ childId ]

#         #if child's starting position should change, add to queuw
#         if( childStart < startingPosition + l ):
#             startingPositions[ childId ] = startingPosition + l   #shift child over
#             queue.put( children[j]["to"] )
#             # print "shove over child: %d pushed by: %d from: %d to: %d  visited: %d" %( childId, b, childStart, startingPositions[ children[j]["to"] ], visited[childId] )
#         elif( visited[childId] == 0 ):
#             queue.put( children[j]["to"] )

#         # distFromEnd = startingPositions[childId] - (startingPositions[b]+l) 
#         # if(distFromEnd< minDistToChild):
#         #     minDistToChild = distFromEnd
          
#         #check to make sure parent is in range of children
#     # if( minDistToChild > 0 ):
#     #     startingPositions[b] = startingPositions[childId] - l
#         ## queue.put(b)


# print "final positions : "
# for i in range(0, len(blocks)):
#     print "block %d, starts at %d" %( i, startingPositions[i] )



###### ADDED THIS TO COMPUTE CONSENSUS ORDERS
## lets start with consensus order between two species
# blockPlaced = [False]*len(blocks) 

# newIdx = len(species)
# species.append( {"id": newIdx, "name": "consensus1", "contigList": [] , "maxNT": 0, "orderedBlocks": [] } )

# for i in range(0, len( species[0]["orderedBlocks"] )):
#     species[newIdx]["orderedBlocks"].append( species[0]["orderedBlocks"][i] )
#     blockPlaced[ species[0]["orderedBlocks"][i] ] = True

# for i in range( 1, len( species[1]["orderedBlocks"] )-1 ):
#     block = species[1]["orderedBlocks"][i]
#     parent = species[1]["orderedBlocks"][i-1]
#     child = species[1]["orderedBlocks"][i+1] 

#     if( !blockPlaced[block] ){
#         parentIdx = len(blocks)
#         childIdx = len(blocks)

#         if( blockPlaced[parentIdx] || blockPlaced[childIdx] )

#         for()

#     }

#     species[newIdx]["orderedBlocks"].append( species[0]["orderedBlocks"][i] )
#     blockPlaced[ species[0]["orderedBlocks"][i] ] = 1

# print species[newIdx]["orderedBlocks"]



## output the original dag








for i in range( 0, len(overlapsWith) ):
    for j in range( i, len(overlapsWith) ):
        if( i != j ):
            if( startingPositions[j] >= startingPositions[i] ):
                if( startingPositions[j] < startingPositions[i] + maxLengths[i] ):
                    overlapsWith[i].append(j)
                    overlapsWith[j].append(i)
            elif( startingPositions[j] <= startingPositions[i] ):
                if( not( startingPositions[j]+ maxLengths[j] < startingPositions[i]) ):
                    overlapsWith[i].append(j)
                    overlapsWith[j].append(i)




for i in range( 0, len(overlapsWith) ):
    for j in range(0, len(overlapsWith[i] ) ) :
        globalSegmentDepth[i] = globalSegmentDepth[i]+ blockSize[j]


# print "input maxLengths:"
# print maxLengths
# 
# print "------number of lines in each block:"
# print blockSize
# print "------depth of each block from start"
# print blockDepth
# print "------- overlaps with:"
# print overlapsWith
# print "------- segmentDepth"
# print globalSegmentDepth 




# #add data to the input data structure
for i in range(0, len(blocks) ):
    blocks[i]["startingPosition"] = startingPositions[i]
    # blocks[i]["overlapsWith"] = overlapsWith[i]
    if( blocks[i]["startingPosition"] + blocks[i]["maxLength"] > lengthOfCoordinateSystem ):
        lengthOfCoordinateSystem = blocks[i]["startingPosition"] + blocks[i]["maxLength"] 

# print "length of coordinate system = %d" % lengthOfCoordinateSystem

# print "ordered blocks and their start position: "
# print "ordered blocks: "
# for i in range(0, len(species) ):
#     print "species: %d" % i
#     for j in range(0, len(species[i]["orderedBlocks"]) ):
#         print "id: %d start: %d end: %d" %( species[i]["orderedBlocks"][j], startingPositions[ species[i]["orderedBlocks"][j] ], startingPositions[ species[i]["orderedBlocks"][j] ]+ maxLengths[ species[i]["orderedBlocks"][j] ] )




######## CONSENSUS order
placed = []
blockIdToOrder = []
for i in range(0, len(blocks)):
    placed.append(False)
    blockIdToOrder.append(-1)


orderedBlocks1 = species[0]["orderedBlocks"]
orderedBlocks2 = species[1]["orderedBlocks"]
orderedBlocks3 = species[2]["orderedBlocks"]
consensusOrder = [] 
 
for i in range(0, len(orderedBlocks1)):
    placed[ orderedBlocks1[i] ] = True
    blockIdToOrder[ orderedBlocks1[i]] = i
    consensusOrder.append(orderedBlocks1[i])


#### number 2
for i in range(1, len(orderedBlocks2)):
    thisBlock = orderedBlocks2[i]
    if not placed[thisBlock] :
        parent = orderedBlocks2[i-1]
        child = orderedBlocks2[i+1]

        reorder = False
        spliceIdx = len(consensusOrder)

        if placed[parent] and placed[child]:
            if blockIdToOrder[parent] < blockIdToOrder[child]:
                spliceIdx = blockIdToOrder[parent]+1
                consensusOrder.insert( spliceIdx, thisBlock )
                reorder = True
            elif blockIdToOrder[parent] > blockIdToOrder[child]: #invert
                spliceIdx = blockIdToOrder[parent]  #place  
                consensusOrder.insert( spliceIdx, thisBlock )
                reorder = True             
        elif placed[parent]:
            spliceIdx = blockIdToOrder[parent]+1
            consensusOrder.insert( spliceIdx, thisBlock )
            reorder = True
        elif placed[child]:
            spliceIdx = blockIdToOrder[child]  #place  
            consensusOrder.insert( spliceIdx, thisBlock )
            reorder = True      


        if reorder:
            for i in range(0, len(consensusOrder)):
                blockId = consensusOrder[i]
                placed[blockId] = True
                blockIdToOrder[blockId] = i 

### number 3
for i in range(1, len(orderedBlocks3)):
    thisBlock = orderedBlocks3[i]
    if not placed[thisBlock] :
        parent = orderedBlocks3[i-1]
        child = orderedBlocks3[i+1]

        reorder = False
        spliceIdx = len(consensusOrder)

        if placed[parent] and placed[child]:
            if blockIdToOrder[parent] < blockIdToOrder[child]:
                spliceIdx = blockIdToOrder[parent]+1
                consensusOrder.insert( spliceIdx, thisBlock )
                reorder = True
            elif blockIdToOrder[parent] > blockIdToOrder[child]: #invert
                spliceIdx = blockIdToOrder[parent]  #place  
                consensusOrder.insert( spliceIdx, thisBlock )
                reorder = True  
        elif placed[parent]:
            spliceIdx = blockIdToOrder[parent]+1
            consensusOrder.insert( spliceIdx, thisBlock )
            reorder = True
        elif placed[child]:
            spliceIdx = blockIdToOrder[child]  #place  
            consensusOrder.insert( spliceIdx, thisBlock )
            reorder = True              

        if reorder:
            for i in range(0, len(consensusOrder)):
                blockId = consensusOrder[i]
                placed[blockId] = True
                blockIdToOrder[blockId] = i 

# print consensusOrder 
maxNT = 0
for i in range(0, len(consensusOrder)):
    maxNT = maxNT + blocks[ consensusOrder[i] ]["maxLength"]



id = len(species)
species.append( {"id": id, "name": "consensus1", "contigList": [-1] , "maxNT": maxNT, "orderedBlocks": consensusOrder} )



allJson = {}
allJson["blocks"] = blocks
allJson["segments"] = segments
allJson["contigs"] = contigs
allJson["species"] = species
allJson["lengthOfCoordinateSystem"] = lengthOfCoordinateSystem
allJson["unplaced"] = len(consensusOrder)-len(blocks)


data_string = json.dumps(allJson, sort_keys=True, indent=2)
print data_string


        
