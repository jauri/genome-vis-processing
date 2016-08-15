import json
import Queue
from sets import Set

#build a block
def createBlock(idx, score, numberAligned, segList, maxLength):
    return {'id': idx, 'score': score, 'numberAligned': numberAligned, 'segList': segList, "maxLength": maxLength}

def getSpeciesId(speciesName):
    if speciesName in speciesToIdHM:
        return speciesToIdHM[speciesName]
    else:
        id = len(speciesToIdHM) 
        species.append( {"id": id, "name": speciesName, "contigList": [] , "maxNT": 0} )
        speciesToIdHM[speciesName] = id
        return id

def getContigId(contigName, speciesId):
    if contigName in contigsToIdHM:
        return contigsToIdHM[contigName]
    else:
        id = len(contigsToIdHM) 
        contigs.append( {"id": id, "name": contigName, "segmentList": [], "orderedBlocks": [] , "contigLength": 0, "speciesId": speciesId} )
        contigsToIdHM[contigName] = id
        return id

#add segment to the contig in order of its start position 
def addSegmentToContig(segment, contigId):
    #first add segment to contig
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
            insertSegmentAtIdx = i #i plus one, because you want it to become the 'i+1'
            # print "segStart = %d, other start %d "  % (segStart, segments[ segArray[i] ]["start"])

            found = True
            break
        # print "or don't insert, append %d segStart = %d, other start %d "  % (segBlock, segStart, segments[ segArray[i] ]["start"])
    

        
    if( not found ):
        segArray.append(segId) 
        orderedBlocks.append(segBlock)
    else:
        segArray.insert(insertSegmentAtIdx, segId ) 
        orderedBlocks.insert( insertSegmentAtIdx, segBlock )

        #why is it here and not in the loop?  
        #so that if it needs to be first, because it is lowest start value, it'll have value 0
    contigs[contigId]["segmentList"] = segArray #reassign this list back to the contig 

    #then produce ordered list of blocks 
    contigs[contigId]["orderedBlocks"] = orderedBlocks

    #debugging... 
    # if( segBlock == 24 ):
    #     print "positioned between segments: %d and %d at idx %d " %( segArray[insertSegmentAtIdx-1], segArray[insertSegmentAtIdx+1], insertSegmentAtIdx )
    #     print segments[segArray[insertSegmentAtIdx-1]]
    #     print segments[segArray[insertSegmentAtIdx]]
    #     print segments[segArray[insertSegmentAtIdx+1]]

    # if( segment["speciesId"] == 2 and len(orderedBlocks) > 2 ):
    #     print "%d positioned between blocks: %d and %d at idx %d array length = %d" %( orderedBlocks[ insertSegmentAtIdx ], orderedBlocks[insertSegmentAtIdx-1], orderedBlocks[insertSegmentAtIdx+1], insertSegmentAtIdx, len(orderedBlocks))
    #     print "%d start: %d" % (orderedBlocks[ insertSegmentAtIdx-1 ], segments[segArray[insertSegmentAtIdx-1]]["start"])
    #     print "%d start: %d" % (orderedBlocks[ insertSegmentAtIdx ], segments[segArray[insertSegmentAtIdx]]["start"])
    #     print "%d start: %d" % (orderedBlocks[ insertSegmentAtIdx+1 ], segments[segArray[insertSegmentAtIdx+1]]["start"])
    # print "current list: "
    # print orderedBlocks


#not using anymore
# def addToSpecies( speciesId, speciesName, contigId, segmentNumber ):
#   if len(species) >= speciesId:
#       species.append( {"speciesId": speciesId, "speciesName": speciesName, "contigs":[], "segments": [] })
#       species[speciesId]["contigs"].append(contigId)
#       species[speciesId]["segments"].append(segmentNumber)
#   else:
#       species[speciesId]["segments"].append(segmentNumber) 

#       found = False
#       for i in range(0, len(species[speciesId]["contigs"]) ) :
#           if species[speciesId]["contigs"] == contigId:
#               found = True
#       if not found:
#           species[speciesId]["contigs"].append(contigId)

# def addToContigs( contigId, contigName, speciesId, segmentNumber ):
#   if len(contigs) >= contigId:
#       contigs.append( {"contigId": contigId, "contigName": contigName, "species":speciesId, "segments": [] } )
#       contigs[contigId]["segments"].append(segmentNumber)
#   else:
#       contig[contigId]["segments"].append(segmentNumber)


###########-----------------BEGIN! 

#------input parameters
runName = "ecoli"
mafFileName = "inputData/EcoliMugsyOut2.maf"
f = open( mafFileName)

#------data strutures to output
blocks = []  #stores the aligned blocks of segments
segments = [] #just the segment info for each genome
contigs = [] #which contigs are present for the strains
species = [] #which species are represented 
allJson = {}

#------- to help with lookups
speciesToIdHM = {}
contigsToIdHM = {} 

#------loop through file, this is the id assigned to blocks and segments 
blockNumber = 0
segmentNumber = 0

#------ Create 0th block, which is empty and represents the start
blocks.append( createBlock(0, 0, 0, {}, 0) )

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

blockNumber = blockNumber+1 
#and an ending block
blocks.append( createBlock((blockNumber), 0, len(species), {}, 0) )  

#fix the starting block 
blocks[0]["numberAligned"] = len(species) 

#add contig to species, ordered by length of contig
for i in range(0, len(contigs) ):
    speciesId = contigs[i]["speciesId"]
    # print "species %d" % speciesId 

    insertIdx = 0
    contigList = species[speciesId]["contigList"] 
    # print "contig list so far: "
    # print contigList 

    contigLen = contigs[i]["contigLength"]

    for j in range(0, len( contigList) ):
        otherLen = contigs[ contigList[j] ]["contigLength"]
        if contigLen > otherLen: #reverse sort, longest first
            insertIdx = j
            break

    contigList.insert( insertIdx, contigs[i]["id"] )
    species[speciesId]["contigList"] = contigList 
    species[speciesId]["maxNT"] = species[speciesId]["maxNT"] + contigs[i]["contigLength"] 
    # print "contig list after add: "
    # print contigList 



#create ordered list of blocks for the species
#longest contig first 
for i in range(0, len(species) ):
    contigList = species[i]["contigList"]
    orderedBlocks = [0]  
    for j in reversed(contigList): #range(0, len( contigList) ):  
        orderedBlocks = orderedBlocks + contigs[ j ]["orderedBlocks"]  

        # orderedBlocks = orderedBlocks + contigs[ contigList[j] ]["orderedBlocks"]  
    # orderedBlocks = orderedBlocks + contigs[ contigList[0] ]["orderedBlocks"]
    orderedBlocks.append( blockNumber )
    species[i]["orderedBlocks"] = orderedBlocks

# print "ordered blocks: "
# for i in range(0, len(species) ):
#     print "species: %d" % i
#     print species[i]["orderedBlocks"]

# print "ordered blocks: "
# for i in range(0, len(species) ):
#     print "species: %d" % i
    




# allJson["blocks"] = blocks
# allJson["segments"] = segments
# allJson["contigs"] = contigs
# allJson["species"] = species

# data_string = json.dumps(allJson, sort_keys=True, indent=2)
# print  data_string


##############################################################
##############################################################

genomesAndAlignments = {}
maxLengths = [] 

startingPositions = []
blockSize = []
blockDepth = []
overlapsWith = []
globalSegmentDepth = []

lengthOfCoordinateSystem = 0

#populate the array storing blocks and strains
for i in range(0, len(species) ):
    genomesAndAlignments[ species[i]["id"] ] = species[i]["orderedBlocks"] 
    for j in range(0, len( genomesAndAlignments[ species[i]["id"] ] ) ) :
        genomesAndAlignments[ species[i]["id"] ][j] = int( genomesAndAlignments[ species[i]["id"] ][j] )

numBlocks = len(blocks)
numGenomes = len(species)

#populate the max length array
for i in range(0, len(blocks) ):
    maxLengths.append( int(blocks[i]["maxLength"]) )
    startingPositions.append(0)
    blockSize.append(0)
    blockDepth.append(0)
    overlapsWith.append([])
    globalSegmentDepth.append(0)


#here is where I store the initial data as an adjacency list
#this is the unfiltered graph
adjLists = []
for i in range(0, numBlocks):
    arr = []
    adjLists.append(arr)

allEdges = [] #this is used to get the sorted list of edges by weight for the greedy algo

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
def addEdgeWithWeight(g, f, t, w):
    found = False 
    for i in range( len(g[f]) ):
        if( g[f][i]["to"] == t ):
             g[f][i]["w"] = w
             found = True
    if not found: 
        g[f].append({"to": t, "w": w, "from": f})

#remove a specific edge (eg when it is in a cycle)
def removeEdge( g, f, t):
    for i in range( len(g[f]) ):
        if( g[f][i]["to"] == t ):
            g[f].pop(i)


root = 0

#populate the adjacecy list with initial data
for key in genomesAndAlignments:
    arr = genomesAndAlignments[key]

    for i in range(len(arr)-1):
        addEdge( adjLists, arr[i], arr[i+1] )

# print len(adjLists)

#populate the all edges dataset and the block size (really the number of edges into and out of node)
for i in range(0, len(adjLists) ):
    # print i
    # print "len adjLists[i]= %d" % len (adjLists[i]) 
    # print "len block size= %d" % len(blockSize)
    for j in range( 0, len (adjLists[i]) ):
        allEdges.append( { "from": i, "to": adjLists[i][j]["to"], "w": adjLists[i][j]["w"] })
        blockSize[ i ] = blockSize[ i ] + adjLists[i][j]["w"]

# blockSize[len(blockSize)-1] = 5  #force the last one to be 5, because has no outgoing edges

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

#get the sorted array
sortedArr = sorted(allEdges, cmp=compareWeights)

#print it just to check
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
maxNoCycle = []
possibleOrphanList = [] 
for i in range(0, numBlocks):
    arr = []
    maxNoCycle.append(arr)

#loop through edges in order of weights
for i in range( len(sortedArr) ):
    e = sortedArr[i]

    addEdgeWithWeight(maxNoCycle, e["from"], e["to"], e["w"])

    # print "ADDED %d, %d" % (e["from"], e["to"] )

    if( isCyclic(maxNoCycle) ):
        removeEdge(maxNoCycle, e["from"], e["to"] ) #empty right now
        # print "REMOVED f: %d, t: %d, w: %d" % (e["from"], e["to"], e["w"])
        possibleOrphanList.append( {"to": e["to"], "from": e["from"], "w": e["w"]} )


#loop through possible orphans and determine if they really are orphans
# actualOrphans = []
# for i in range(0, len(possibleOrphanList) ):
#     orphan = possibleOrphanList[i]

#     for j in range(0, len(maxNoCycle) ):
#         e = 

# for i in range(0, len(maxNoCycle) ):
#     arr = maxNoCycle[i] #children
#     for j in range(0, len(arr)):
#         e = arr[j]

#         foundParent = False 
#         for k in range(0, len(possibleOrphanList)):
#             if e["to"] == possibleOrphanList[k]:
#                 foundParent = True

#         if not foundParent:
#             actualOrphans


### START HERE!!!
actualOrphans = Set() 
for i in range(0, len(possibleOrphanList)):
    foundParent = False
    parent = 0

    for j in range(0, len(maxNoCycle) ):
        arr = maxNoCycle[j]  #children

        for k in range(0, len(arr)):
            e = arr[k]

            if e["to"] == possibleOrphanList[i]["to"] :
                foundParent = True

    if not foundParent:
        actualOrphans.add(possibleOrphanList[i]["to"] )



# print "------- no cycle max weights"
# print maxNoCycle


# print "-------- orphans "
# print actualOrphans

#add orphans to list
orphansFinal = [] 
while len(actualOrphans) > 0:
    orphansFinal.append( actualOrphans.pop() )


# print isCyclic( maxNoCycle )


# print "------- positioning segments"
# def recursivePosition(blockId, parentEndPosition):
#     if parentEndPosition > startingPositions[blockId]:
#         startingPositions[blockId] = parentEndPosition

#     arr = maxNoCycle[i] #children
#     l = maxLengths[i] #max length of this segment
#     # for i in range(0, len( arr)):
#     #     recursivePosition( arr[j]["to"], l = )


#only works because 0th is the root 
# blockDepth[0] = 0
# for i in range(0, len(maxNoCycle) ):
#     arr = maxNoCycle[i] #children
#     l = maxLengths[i] #max length of this segment
#     for j in range(0, len(arr)):
#         e = arr[j]
#         if( startingPositions[ e["to"] ] < maxLengths[ e["from"] ] + startingPositions[ e["from"] ] ):
#             startingPositions[ e["to"] ] = maxLengths[ e["from"] ] + startingPositions[ e["from"] ]         
#         if( blockDepth[ e["to"] ] < blockDepth[ e["from"] ] + 1 ):
#             blockDepth[ e["to"] ] = blockDepth[ e["from"] ] + 1


##########FIRST PASS POSITIONING 
root = 0
queue = Queue.Queue() 
queue.put(root)

visited = [0]*len(blocks)


while not queue.empty():
    b = queue.get() 
    visited[b] = visited[b]+1

    children = maxNoCycle[b] #children
    l = maxLengths[b] #max length of this segment
    startingPosition = startingPositions[b]

    minDistToChild = 9999999999999999

    for j in range(0, len(children) ):
        childId = children[j]["to"] 
        childStart = startingPositions[ childId ]

        #if child's starting position should change, add to queuw
        if( childStart < startingPosition + l ):
            startingPositions[ childId ] = startingPosition + l   #shift child over
            queue.put( children[j]["to"] )
            # print "shove over child: %d pushed by: %d from: %d to: %d  visited: %d" %( childId, b, childStart, startingPositions[ children[j]["to"] ], visited[childId] )
        elif( visited[childId] == 0 ):
            queue.put( children[j]["to"] )

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


##### RESOLVE ORPHANS- place after their parents 
rootList = []
for i in range(0, len(orphansFinal)):
    rootList.append(orphansFinal[i])

    parentList = []
    parentPosition = 0
    for j in range(0, len( possibleOrphanList) ): 
        if possibleOrphanList[j]["to"] == orphansFinal[i]:
        #     parentList.append(possibleOrphanList[j]["from"] )
            parentEnd = startingPositions[ possibleOrphanList[j]["from"] ] + maxLengths[ possibleOrphanList[j]["from"] ]
            if parentEnd > parentPosition:
                parentPosition =  parentEnd

    startingPositions[orphansFinal[i]] = parentEnd





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

while not queue.empty():
    b = queue.get() 
    visited[b] = visited[b]+1

    children = maxNoCycle[b] #children
    l = maxLengths[b] #max length of this segment
    startingPosition = startingPositions[b]

    minDistToChild = 9999999999999999

    for j in range(0, len(children) ):
        childId = children[j]["to"] 
        childStart = startingPositions[ childId ]

        #if child's starting position should change, add to queuw
        if( childStart < startingPosition + l ):
            startingPositions[ childId ] = startingPosition + l   #shift child over
            queue.put( children[j]["to"] )
            # print "shove over child: %d pushed by: %d from: %d to: %d  visited: %d" %( childId, b, childStart, startingPositions[ children[j]["to"] ], visited[childId] )
        elif( visited[childId] == 0 ):
            queue.put( children[j]["to"] )

        # distFromEnd = startingPositions[childId] - (startingPositions[b]+l) 
        # if(distFromEnd< minDistToChild):
        #     minDistToChild = distFromEnd
          
        #check to make sure parent is in range of children
    # if( minDistToChild > 0 ):
    #     startingPositions[b] = startingPositions[childId] - l
        ## queue.put(b)


# print "final positions : "
# for i in range(0, len(blocks)):
#     print "block %d, starts at %d" %( i, startingPositions[i] )



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


allJson = {}
allJson["blocks"] = blocks
allJson["segments"] = segments
allJson["contigs"] = contigs
allJson["species"] = species
allJson["lengthOfCoordinateSystem"] = lengthOfCoordinateSystem


data_string = json.dumps(allJson, sort_keys=True, indent=2)
print data_string


        
