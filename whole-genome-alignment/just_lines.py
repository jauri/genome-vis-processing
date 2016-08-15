import json
import sys
import Queue
from sets import Set

#------------ Functions to create initial data types
def createBlock(idx, score, numberAligned, segList, maxLength):
    return {'id': idx, 'score': score, 'numberAligned': numberAligned, 'segList': segList, "maxLength": maxLength, "compressedNextTo": 0, "compressedStart": 0, "neighborOrder": 0, "isGapBlock": False}

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
        contigs.append( {"id": id, "name": contigName, "segmentList": [], "orderedBlocks": [] , "contigLength": 0, "speciesId": speciesId, "order": 0} )
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
mafFileName = "inputData/seven_ecoli.maf"#"twenty_five_ecoli.maf"#"inputData/seven_ecoli.maf" #"inputData/aBaumannii.maf" #"inputData/EcoliMugsyOut2.maf"
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

        speciesName = tokens[1].split('.')[0] #which species is this from (maf file is missing part of species name)
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


#-------- give contigs an order value, for coloring by its order within species
for i in range(0, len(species) ):
    contigList = species[i]["contigList"]

    for j in range(0, len(contigList) ):
        contigs[ contigList[j] ]["order"] = j;

#--------- Order blocks of aligned segments for each species 
# This gives each species preferred order for the aligned blocks 

for i in range(0, len(species) ):
    contigList = species[i]["contigList"]

    #create contig order

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

#### REMOVED_1 #######

#### END REMOVED_1 ###########

#####  REMOVED_2 #########

##### END REMOVED_2 #########


######## CONSENSUS order
placed = []
blockIdToOrder = []
for i in range(0, len(blocks)):
    placed.append(False)
    blockIdToOrder.append(-1)

consensusOrder = [] 
orderedBlocks1 = species[0]["orderedBlocks"]
for i in range(0, len(orderedBlocks1) ):
    placed[ orderedBlocks1[i] ] = True
    blockIdToOrder[ orderedBlocks1[i]] = i
    consensusOrder.append(orderedBlocks1[i])


for s in range(1, len(species)):
    orderedBlocks_s = species[s]["orderedBlocks"]

    for i in range(1, len(orderedBlocks_s)):
        thisBlock = orderedBlocks_s[i]
        if not placed[thisBlock] :
            parent = orderedBlocks_s[i-1]
            child = orderedBlocks_s[i+1]

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
# for i in range(1, len(orderedBlocks3)):
#     thisBlock = orderedBlocks3[i]
#     if not placed[thisBlock] :
#         parent = orderedBlocks3[i-1]
#         child = orderedBlocks3[i+1]

#         reorder = False
#         spliceIdx = len(consensusOrder)

#         if placed[parent] and placed[child]:
#             if blockIdToOrder[parent] < blockIdToOrder[child]:
#                 spliceIdx = blockIdToOrder[parent]+1
#                 consensusOrder.insert( spliceIdx, thisBlock )
#                 reorder = True
#             elif blockIdToOrder[parent] > blockIdToOrder[child]: #invert
#                 spliceIdx = blockIdToOrder[parent]  #place  
#                 consensusOrder.insert( spliceIdx, thisBlock )
#                 reorder = True  
#         elif placed[parent]:
#             spliceIdx = blockIdToOrder[parent]+1
#             consensusOrder.insert( spliceIdx, thisBlock )
#             reorder = True
#         elif placed[child]:
#             spliceIdx = blockIdToOrder[child]  #place  
#             consensusOrder.insert( spliceIdx, thisBlock )
#             reorder = True              

#         if reorder:
#             for i in range(0, len(consensusOrder)):
#                 blockId = consensusOrder[i]
#                 placed[blockId] = True
#                 blockIdToOrder[blockId] = i 

# print consensusOrder 
maxNT = 0
for i in range(0, len(consensusOrder)):
    maxNT = maxNT + blocks[ consensusOrder[i] ]["maxLength"]


lengthOfCoordinateSystem = maxNT
# length of the whole coordinate system 
# for i in range(0, len(blocks) ):
#     blocks[i]["startingPosition"] = startingPositions[i]
#     if( blocks[i]["startingPosition"] + blocks[i]["maxLength"] > lengthOfCoordinateSystem ):
#         lengthOfCoordinateSystem = blocks[i]["startingPosition"] + blocks[i]["maxLength"] 

### DEBUGGING 
# print "length of coordinate system = %d" % lengthOfCoordinateSystem


####### COMPRESSED COORDINATE 
lastNeighbor = [0]*len(species) #for each species, as you go through the blocks, determine the proper 'neighbor'
lastPosition = [0]*len(species) #keeping a running log of the coordinates for the blocks, as you position them next to their neighbor
countNeighbors = [0]*len(blocks) #keep track of how many blocks are directly adjacent to another block, for y positioning

speciesForTimeLine = []; 
blocksForTimeLine = []; 

for i in range(0, len(species) ):
    sp = species[i]
    speciesForTimeLine.append( {"id": sp["id"], "name": sp["name"], "contigList": sp["contigList"], "timeline": []} )
    # print "species: "
    # print speciesForTimeLine[i]



#layout blocks in x 
for i in range(0, len(consensusOrder) ): #loop through the blocks in order
    # print consensusOrder[i]

    b = blocks[ consensusOrder[i] ] #get a block
    segs = b["segList"]  # get its segments- this will tell you which species are in this block

    greatestLastPosition = 0; 
    greatestLastPositionIdx = 0;

    #for each species, find its latest position in x
    #determine which species 'pushes' this block over the furthest
    #the block should fall after the last position of all species in the block
    for s in range(0, len(segs) ): #loop through these segments 
        # print "   seg:  %d" % segs[s]["id"]
        sp = segments[ segs[s]["id"] ]["speciesId"] #get the species
        # print "        species: %d " % sp 
        if lastPosition[ sp ] > greatestLastPosition: #look up the 'last position' of this species in the running tally
            greatestLastPosition = lastPosition[ sp ] #if it is greatest, note this 
            greatestLastPositionIdx = sp  

    #at this point we have stored the 'neighbor' for this block and its position in 'time '
   
    # print " block: %d has max length %d and starts at: %d next to: %d" % ( consensusOrder[i] ,  blocks[ consensusOrder[i] ]["maxLength"] , greatestLastPosition , lastNeighbor[greatestLastPositionIdx] )

    #old storage method 
    blocks[consensusOrder[i]]["compressedStart"] = greatestLastPosition;
    blocks[consensusOrder[i]]["compressedNextTo"] = lastNeighbor[greatestLastPositionIdx]
    blocks[consensusOrder[i]]["neighborOrder"] = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ]
    countNeighbors[ lastNeighbor[greatestLastPositionIdx] ] = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ]  + 1; #increment

    #new storage method - for the species 
    for s in range(0, len(segs) ): #for each species present 
        # print "   seg:  %d" % segs[s]["id"]
        sp = segments[ segs[s]["id"] ]["speciesId"] #get species id

        timeline = speciesForTimeLine[ sp ]["timeline"] #get its 'timeline' 
        timeline.append( {"block": consensusOrder[i], "start": blocks[consensusOrder[i]]["compressedStart"], "size": blocks[consensusOrder[i]]["maxLength"], "slot":  blocks[consensusOrder[i]]["neighborOrder"], "speciesOrder": s, "speciesInBlock": len(segs) } );
        speciesForTimeLine[sp]["timeline"] = timeline
        

    # print "idx: %d neighbor: %d order: %d" % (consensusOrder[i], lastNeighbor[greatestLastPositionIdx], blocks[consensusOrder[i]]["neighborOrder"]) 

    #update


    for s in range(0, len(segs) ):
        sp = segments[ segs[s]["id"] ]["speciesId"]
        lastPosition[sp] = greatestLastPosition + blocks[ consensusOrder[i] ]["maxLength"] 
        lastNeighbor[sp] = consensusOrder[i]


compressedMaxLength = 0; 
for s in range(0, len(species)):
    if lastPosition[s] > compressedMaxLength:
        compressedMaxLength = lastPosition[s]



print speciesForTimeLine[0]

# #### add the empty blocks
# speciesPresent = [0]*len(species) #keep track of which species have been encountered for each set of neighbor blocks
# currentNeighbor = blocks[ consensusOrder[0] ]["compressedNextTo"]; #which neighbor is 'current '
# # count = 0 
# for i in range(0, len(consensusOrder)): #loop through the consensus order
#     b = blocks[ consensusOrder[i] ] #get a block
#     segs = b["segList"] #get the inner segments 
#     # count = count + 1

#     if currentNeighbor != b["compressedNextTo"] : #check- is this a new neighbor set? 

#         # if so, resolve it 
#         absentSpecies= [] # which ones are missing 
#         for j in range(0, len(speciesPresent)):
#             if speciesPresent[j] == 0:
#                 absentSpecies.append(j) 

#         #if any are missing 
#         if( len(absentSpecies) > 0 ):
#             idx = len(blocks) #make a new id for it
#             segList = [] #and new set of segments  
#             for k in range(0, len(absentSpecies) ): # for each species 
#                 s = {"id": len(segments), "speciesId": absentSpecies[k], "contigId": -1, "size": 100, "start": 0, "reverseCompliment": False};
#                 segList.append( {"id": s["id"]} ) #add empty segment
#                 segments.append(s) #include in list...?  Maybe...? =/

#             #create a new block to hold these absent species     
#             order = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ]
#             newBlock = {'id': idx, 'score': 0, 'numberAligned': len(absentSpecies), 'segList': segList, "maxLength": 100, "compressedNextTo": currentNeighbor, "compressedStart": 0, "neighborOrder": order+1, "isGapBlock": True}
#             consensusOrder.insert( i, newBlock["id"] )
#             blocks.append(newBlock)

#             countNeighbors[ lastNeighbor[greatestLastPositionIdx] ] = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ] + 1
#             print "idx: %d neighbor: %d order: %d " % (idx, currentNeighbor, order) 
#             # print "################## new block:"
#             # print newBlock


#         currentNeighbor = b["compressedNextTo"]
#         count = 0 
#         for sp in range(0, len(speciesPresent)):
#             speciesPresent[sp] = 0 

#         #i = i +1 


#     for s in range(0, len(segs) ):
#         # print "   seg:  %d" % segs[s]["id"]
#         sp = segments[ segs[s]["id"] ]["speciesId"]
#         speciesPresent[sp] = 1


# print consensusOrder

id = len(species)
species.append( {"id": id, "name": "consensus1", "contigList": [-1] , "maxNT": maxNT, "orderedBlocks": consensusOrder} )

# print "CONSENSUS:"
# print consensusOrder


allJson = {}
allJson["blocks"] = blocks
allJson["segments"] = segments
allJson["contigs"] = contigs
allJson["species"] = species
allJson["lengthOfCoordinateSystem"] = lengthOfCoordinateSystem
allJson["unplaced"] = len(consensusOrder)-len(blocks)
allJson["idOfConsensusSpecies"] = len(species)-1
allJson["compressedMaxLength"] = compressedMaxLength

data_string = json.dumps(allJson, sort_keys=True, indent=2)
# print data_string


        
