import json
import sys
import Queue
from sets import Set

#------------ Functions to create initial data types
def createBlock(idx, score, numberAligned, segList, maxLength):
    return {'id': idx, 'score': score, 'numberAligned': numberAligned, 'segList': segList, "maxLength": maxLength, "compressedNextTo": 0, "compressedStart": 0, "neighborOrder": 0, "isGapBlock": False, "neighborInFrontList": [], "neighborsBehindList": []}

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

fixedWidthMode = False

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

        if fixedWidthMode:
            size = 100000

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
lastNeighbor = [0]*len(species)
lastPosition = [0]*len(species)
countNeighbors = [0]*len(blocks)
lastPositionFixedWidth = [0]*len(blocks)

for i in range(0, len(consensusOrder) ):
    # print consensusOrder[i]

    b = blocks[ consensusOrder[i] ]
    segs = b["segList"]

    greatestLastPosition = 0;
    greatestLastPositionIdx = 0;
    # greatestLastFixedWidthPosition = 0;
    # greatestLastFixedWidthPositionIdx = 0; 
    for s in range(0, len(segs) ):
        # print "   seg:  %d" % segs[s]["id"]
        sp = segments[ segs[s]["id"] ]["speciesId"]
        # print "        species: %d " % sp 
        if lastPosition[ sp ] > greatestLastPosition:
            greatestLastPosition = lastPosition[ sp ]
            greatestLastPositionIdx = sp 

        #see whether to add to neighbors behind- avoid duplicates
        present = False
        neighborsBehindList = blocks[ consensusOrder[i] ]["neighborsBehindList"]
        for n in range(0, len( neighborsBehindList )):
            if neighborsBehindList[n] == lastNeighbor[sp] :
                present = True 
        
        if not present:
            neighborsBehindList.append(lastNeighbor[sp])
            blocks[ consensusOrder[i] ]["neighborsBehindList"] = neighborsBehindList
        # if lastPositionFixedWidth[ sp ] > greatestLastFixedWidthPosition:
        #     greatestLastFixedWidthPosition = lastPositionFixedWidth[ sp ]
        #     greatestLastFixedWidthPositionIdx = sp

    # print " block: %d has max length %d and starts at: %d next to: %d" % ( consensusOrder[i] ,  blocks[ consensusOrder[i] ]["maxLength"] , greatestLastPosition , lastNeighbor[greatestLastPositionIdx] )

    blocks[consensusOrder[i]]["compressedStart"] = greatestLastPosition;
    blocks[consensusOrder[i]]["compressedNextTo"] = lastNeighbor[greatestLastPositionIdx]
    blocks[consensusOrder[i]]["neighborOrder"] = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ]
    # blocks[consensusOrder[i]]["compressedFixedWidth"] = greatestLastFixedWidthPosition; 
    # blocks[consensusOrder[i]]["compressedNextToFixedWidth"] = greatestLastFixedWidthPositionIdx; 
    
    neighborFrontList = blocks[ lastNeighbor[greatestLastPositionIdx] ]["neighborInFrontList"]
    neighborFrontList.append( consensusOrder[i] )
    blocks[ lastNeighbor[greatestLastPositionIdx] ]["neighborInFrontList"] = neighborFrontList
    countNeighbors[ lastNeighbor[greatestLastPositionIdx] ] = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ]  + 1; #increment
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









#### REMOVED! 
# ##  NEW METHOD 
# blockCount = len(blocks)
# for i in range (0, len(consensusOrder)):
#     b = blocks[ consensusOrder[i] ] #get a block

#     frontNeighbors = b["neighborInFrontList"]

#     speciesPresent = [0]*len(species) #keep track of which species have been encountered for each set of neighbor blocks
#     segList = [] 
#     absentSpecies = []

#     if not blocks[consensusOrder[i]]["isGapBlock"]:

#         for j in range(0, len(frontNeighbors) ): #for all neighbors 

#             segs = blocks[frontNeighbors[j]]["segList"] #get the inner segments

#             for s in range(0, len(segs) ): #mark which species are present 
#                 # print "   seg:  %d" % segs[s]["id"]
#                 sp = segments[ segs[s]["id"] ]["speciesId"]
#                 speciesPresent[sp] = 1

#         for k in range(0, len(speciesPresent)): #find absences 
#             if speciesPresent[k] == 0:
#                 absentSpecies.append(k)

#         # print len(absentSpecies)

#         for k in range(0, len(absentSpecies) ): # for each absent species , make a dummy segment 
#             s = {"id": len(segments), "speciesId": absentSpecies[k], "contigId": -1, "size": 1000, "start": 0, "reverseCompliment": False};
#             segList.append( {"id": s["id"]} ) #add empty segment
#             segments.append(s) #include in list...?  Maybe...? =/

#         if len(absentSpecies) > 0 : 
#             # print idx
#             newBlock = {"id": blockCount, 'score': 0, 'numberAligned': len(absentSpecies), 'segList': segList, 'maxLength': 1000,'compressedNextTo': blocks[consensusOrder[i]]["id"], "compressedStart": (blocks[consensusOrder[i]]["compressedStart"]+blocks[consensusOrder[i]]["maxLength"]), "neighborOrder": len(frontNeighbors)+1, "isGapBlock": True, "neighborInFrontList": [] }
#             # print newBlock
#             blocks.append(newBlock)
#             frontNeighbors.append(newBlock["id"])
#             blocks[consensusOrder[i]]["neighborInFrontList"] = frontNeighbors
#             consensusOrder.insert( i+1 , blockCount )
#             blockCount = blockCount+1

#         # i = i+1 #skip



## Position blocks into horizontal slots- like the algorithm
##### Right now done in an arbitrary order, no optimization
slots = []
for i in range(0, len(species)+1):
    slots.append([])

count = 0
for i in range(0, len(consensusOrder)):

    b = blocks[ consensusOrder[i] ]
    bStart = b["compressedStart"]
    bStop = b["compressedStart"] + b["maxLength"]

    count = count + 1
    #find an available slot
    assignedSlot = -1
    j = 0
    blocksAssigned = [] 
    found = False
    for j in range(0, len(slots)):

        blocksAssigned = slots[j]
        
        overlaps = False
        for k in range(0, len(blocksAssigned) ):
            other = blocks[ blocksAssigned[k] ]
            otherStart = other["compressedStart"]
            otherStop = other["compressedStart"] + other["maxLength"]

            if bStart <= otherStart and bStop >= otherStart:
                overlaps = True
            if bStart >= otherStart and bStart <= otherStop:
                overlaps = True


        if not overlaps and not found:
            assignedSlot = j
 
            if j % 2 == 0:
                assignedSlot = len(slots)/2-j/2
            else:
                assignedSlot = len(slots)/2+j/2

            # print assignedSlot

            blocksAssigned.append( consensusOrder[i] )
            found = True
            blocks[consensusOrder[i]]["slot"] = assignedSlot

        blocksAssigned = [] 






#### POSITION SEGMENTS WITHIN THE BLOCK
previous = consensusOrder[1]  ## skip the 0th, it has no segments 
previousSegs = blocks[ previous ]["segList"]

for i in range(1, len(consensusOrder)):

    #classify segments in all blocks
    b = blocks[ consensusOrder[i] ]
    segs = b["segList"]
    neighborsBehindList = b["neighborsBehindList"]

    # print "considering block %d " % b["id"]

    # if consensusOrder[i] == 14:
        # print "block 14: "

    #for each segment, find the neighbor behind
    for s in range(0, len(segs)):  #for all segs 
        latestNeighbor = -1
        latestNeighborPosition = 0
        latestNeighborSlot = -1 
        latestNeighborIdx = 0
        sp = segments[segs[s]["id"]]["speciesId"]
        found = False
        inSlot = -1

        # print "--block: %d" % consensusOrder[i]
        # print "---- segment: %d and species: %s" %( segs[s]["id"], species[sp]["name"] )

        for n in range(0, len(neighborsBehindList)): #look at previous neighbors 

            otherSegs = blocks[ neighborsBehindList[n] ] ["segList"]
            for os in range(0, len(otherSegs)):
                spOs = segments[ otherSegs[os]["id"] ]["speciesId"]

                # print "------------- neighbor: %d and segment: %d species: %s" %( neighborsBehindList[n], otherSegs[os]["id"], species[spOs]["name"] )

                if spOs == sp:
                    # found = True
                    # print "-------------MATCH!-----------"
                    # inSlot = blocks[ neighborsBehindList[n] ]["slot"]
                    if blocks[ neighborsBehindList[n] ]["compressedStart"] + blocks[neighborsBehindList[n]]["maxLength"] > latestNeighborPosition:
                        latestNeighborPosition =  blocks[ neighborsBehindList[n] ]["compressedStart"] + blocks[neighborsBehindList[n]]["maxLength"]
                        latestNeighborSlot = blocks[ neighborsBehindList[n] ]["slot"]#inSlot
                        latestNeighbor = blocks[ neighborsBehindList[n] ]
                        found = True
                        latestNeighborIdx = os
                        # print "-------------IS LATEST-----------"
        #done looping through back neighbors 

        if found:

            if latestNeighborSlot == b["slot"]:
                # print "       static"
                b["segList"][s]["class"] = "static";
                b["segList"][s]["prevIdx"]  = latestNeighborIdx
                b["segList"][s]["prevSlot"] = latestNeighborSlot
                b["segList"][s]["latestNeighborPosition"] = latestNeighborPosition
                b["segList"][s]["latestNeighborBlock"] = latestNeighbor["id"]
            if latestNeighborSlot < b["slot"]:
                # print "      lower in neighbor, rising now"
                b["segList"][s]["class"] = "descending"#"rising";
                b["segList"][s]["prevIdx"]  = latestNeighborIdx
                b["segList"][s]["prevSlot"] = latestNeighborSlot
                b["segList"][s]["latestNeighborPosition"] = latestNeighborPosition
                b["segList"][s]["latestNeighborBlock"] = latestNeighbor["id"]
            if latestNeighborSlot > b["slot"]:
                # print "       higher in neighbor, descending now"
                b["segList"][s]["class"] = "rising"#"descending"
                b["segList"][s]["prevIdx"]  = latestNeighborIdx
                b["segList"][s]["prevSlot"] = latestNeighborSlot
                b["segList"][s]["latestNeighborPosition"] = latestNeighborPosition
                b["segList"][s]["latestNeighborBlock"] = latestNeighbor["id"]      

            # print "segment %s in block %d with slot %d was found in previous neighbor %d in slot %d at idx %d " %( species[segments[segs[s]["id"]]["speciesId"]]["name"], b["id"], b["slot"], latestNeighbor["id"], latestNeighborSlot, latestNeighborIdx )




    blocks[consensusOrder[i]] = b
            # if not found:
            #     print ""


#done with a seg classification, now order them to minimize line crossings 
for i in range(2, len(consensusOrder)):

    b = blocks[ consensusOrder[i] ]
    segs = b["segList"]
    neighborsBehindList = b["neighborsBehindList"]

    reorderedSegs = [] 
    for s in range(0, len(segs)):  #for all segs 
        insertionIdx = len(reorderedSegs)
        found = False
        sp = segments[segs[s]["id"]]["speciesId"]


        ### Update idx within the slot 
        neighborBlock = blocks[ b["segList"][s]["latestNeighborBlock"] ]
        # print neighborBlock["id"]
        otherSegs = neighborBlock["segList"]
        for os in range(0, len(otherSegs)):
            spOs = segments[ otherSegs[os]["id"] ]["speciesId"]

            # print "------------- neighbor: %d and segment: %d species: %s" %( neighborsBehindList[n], otherSegs[os]["id"], species[spOs]["name"] )
            # print spOs 
            if spOs == sp:
                b["segList"][s]["prevIdx"] = os
                # print "SAME .  idx is %d" % os

        # print b["segList"][s]["prevIdx"]

        for rs in range(0, len(reorderedSegs)):
            # print b["segList"][s]
            # print reorderedSegs[rs]
            if not found:
                # print "------- looking at: block %d and segment: %s which was at idx: %d .  Compre to segment %s which was at idx %d " %( b["id"], species[ segments[ b["segList"][s]["id"] ]["speciesId"]]["name"], b["segList"][s]["prevIdx"], species[ segments[ reorderedSegs[rs]["id"] ]["speciesId"]]["name"], reorderedSegs[rs]["prevIdx"] ) 
                #descending: 
                if b["segList"][s]["class"] == "descending" and reorderedSegs[rs]["class"] != "descending": #stop the moment it is no longer descending
                    insertionIdx = rs
                    found = True 
                if b["segList"][s]["class"] == "descending" and reorderedSegs[rs]["class"] == "descending" and b["segList"][s]["prevSlot"] < reorderedSegs[rs]["prevSlot"]:
                    insertionIdx = rs
                    found = True 
                if b["segList"][s]["class"] == "descending" and reorderedSegs[rs]["class"] == "descending" and b["segList"][s]["prevSlot"] == reorderedSegs[rs]["prevSlot"] and b["segList"][s]["latestNeighborPosition"] < reorderedSegs[rs]["latestNeighborPosition"]:
                    insertionIdx = rs
                    found = True
                if b["segList"][s]["class"] == "descending" and reorderedSegs[rs]["class"] == "descending" and b["segList"][s]["prevSlot"] == reorderedSegs[rs]["prevSlot"] and b["segList"][s]["latestNeighborPosition"] == reorderedSegs[rs]["latestNeighborPosition"] and b["segList"][s]["prevIdx"] < reorderedSegs[rs]["prevIdx"]:
                    insertionIdx = rs
                    found = True
                
                #static 
                if b["segList"][s]["class"] == "static" and reorderedSegs[rs]["class"] == "rising":
                    insertionIdx = rs
                    found = True
                # if b["segList"][s]["class"] == "static" and reorderedSegs[rs]["class"] == "static" and b["segList"][s]["latestNeighborPosition"] < reorderedSegs[rs]["latestNeighborPosition"]:
                #     insertionIdx = rs
                #     found = True 
                if b["segList"][s]["class"] == "static" and reorderedSegs[rs]["class"] == "static" and  b["segList"][s]["latestNeighborPosition"] == reorderedSegs[rs]["latestNeighborPosition"] and b["segList"][s]["prevIdx"] < reorderedSegs[rs]["prevIdx"]:
                    insertionIdx = rs
                    found = True 

                #rising
                if b["segList"][s]["class"] == "rising" and reorderedSegs[rs]["class"] == "rising" and b["segList"][s]["prevSlot"] < reorderedSegs[rs]["prevSlot"]:
                    insertionIdx = rs
                    found = True
                if b["segList"][s]["class"] == "rising" and reorderedSegs[rs]["class"] == "rising" and b["segList"][s]["prevSlot"] == reorderedSegs[rs]["prevSlot"] and b["segList"][s]["latestNeighborPosition"] < reorderedSegs[rs]["latestNeighborPosition"]:
                    insertionIdx = rs
                    found = True   
                if b["segList"][s]["class"] == "rising" and reorderedSegs[rs]["class"] == "rising" and b["segList"][s]["prevSlot"] == reorderedSegs[rs]["prevSlot"] and b["segList"][s]["latestNeighborPosition"] == reorderedSegs[rs]["latestNeighborPosition"] and b["segList"][s]["prevIdx"] < reorderedSegs[rs]["prevIdx"]:
                    insertionIdx = rs
                    found = True 


        # print "FINAL: block %d  inserting  %s  at index %d in class %s" % ( b["id"], species[ segments[ b["segList"][s]["id"] ]["speciesId"]]["name"] , insertionIdx, b["segList"][s]["class"])
        reorderedSegs.insert( insertionIdx, b["segList"][s] )

    # print "before: " 
    # print segs 
    # print "after: "
    # print reorderedSegs
    # print "--------- "

    b["segList"] = reorderedSegs
    blocks[consensusOrder[i]] = b


#now remove whitespace


#first assign blocks to an initial y coordinate- based on slot
for i in range(0, len(blocks)):
    blocks[i]["yAssigned"] = blocks[i]["slot"]*len(species)*10
    # print blocks[i]["yAssigned"]


#shift y coordinate to line up segments within slots
minYShift = 0
for i in range(0, len(consensusOrder)):
    b = blocks[ consensusOrder[i] ]
    segs = b["segList"]
    neighborsBehindList = b["neighborsBehindList"]

    maxShift = len(species) - len(segs) + 1
    yShift = []
    if len(segs) < len(species):
        for s in range(0, len(segs)):
            if segs[s]["class"] == "static":
                yDiff =  segs[s]["prevIdx"] - s
                if yDiff < maxShift and yDiff >= 0:
                    yShift.append(yDiff)
        
    # if yDiff < maxShift: 
    #     yShift = yDiff  
                # if yDiff < minYShift:
                    # minYShift = yDiff
    consensusShift = 99999999999
    for j in range(0, len(yShift)):
        if yShift[j] < consensusShift:
            consensusShift = yShift[j]
    if consensusShift == 99999999999:
        consensusShift = 0

    # print "maxShift: %d and consensusShift %d " %( maxShift, consensusShift)
    b["yAssigned"] =  b["yAssigned"] + consensusShift*10



#first create 'slot sets' of blocks:
### contiguous blocks in a slot that share a subset of lines 










    # reorderedSegs.append( b["segList"][0] )
    # insertionIdx = 0
    # for s in range(1, len(segs)): 
    #     for rs in range(0, len(reorderedSegs)):
    #         if b["segList"][s]["class"] == "descending": # if this one is descending
    #             if reorderedSegs[rs]["class"] == "descending": #and the other is descending
    #                 if reorderedSegs[rs]      #if they are in the same slot

    #                 if reorderedSegs[rs]["prevIdx"] > b["segList"][s]["prevIdx"]: #insert in front
    #                     reorderedSegs.insert( rs+1 , b["segList"][s] )


    # for s in range(1, len(segs)): 
    #     for rs in range(0, len(reorderedSegs)):
    #         if b["segList"][s]["class"] == "descending":
    #             if reorderedSegs[rs]["class"] == "descending":
    #                 if reorderedSegs[rs]["prevIdx"] > b["segList"][s]["prevIdx"]: #insert in front
    #                     reorderedSegs.insert( rs+1 , b["segList"][s] )
    #                 else:
    #                     reorderedSegs.insert( rs , b["segList"][s] ) #insert behind
    #             else:  # if reorderedSegs[rs]["class"] == "static" or reorderedSegs[rs]["class"] == "rising": #insert behind
    #                 reorderedSegs.insert( rs , b["segList"][s] ) #insert behind
    #         if b["segList"][s]["class"] == "rising":
    #             if reorderedSegs[rs]["class"] == "rising":
    #                 if reorderedSegs[rs]["prevIdx"] < b["segList"][s]["prevIdx"]: #insert in front
    #                     reorderedSegs.insert( rs+1 , b["segList"][s] )
    #                 else:
    #                     reorderedSegs.insert( rs , b["segList"][s] ) #insert behind
    #             #if static or descending, keep going... 
    #         if b["segList"][s]["class"] == "static":
    #             if reorderedSegs[rs]["class"] == "rising":#insert behind rising
    #                 reorderedSegs.insert( rs , b["segList"][s] ) #insert behind
    #             if reorderedSegs[rs]["class"] == "static":
    #                 if reorderedSegs[rs]["prevIdx"] > b["segList"][s]["prevIdx"]: #insert in front
    #                     reorderedSegs.insert( rs , b["segList"][s] )






    # print "***** outcome: slot: %d and block: %d" % (assignedSlot , consensusOrder[i] )       #

    #add to the slot

# for i in range(0, len(slots)):
#     print "slot: %d" % i
#     for j in range(0, len(slots[i])):
#         print "       %d " % slots[i][j]



#add the consensus species for storing the consensus order
id = len(species)
species.append( {"id": id, "name": "consensus1", "contigList": [-1] , "maxNT": maxNT, "orderedBlocks": consensusOrder} )


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
print data_string


        
