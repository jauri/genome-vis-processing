import json
import sys
import Queue
from sets import Set
import pprint
import math

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

#---------------------------------------------------------------------------------------------
#-------READ FILE-----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

#look through maf file
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


#---------------------------------------------------------------------------------------------
#-------CONSENSUS ORDER GENERATION -----------------------------------------------------------
#---------------------------------------------------------------------------------------------

## the current method uses neighbors to position segments
## basically for a given unplaced segment, find a parent (behind) or child (in front)
## insert and repeat until all blocks have been placed

## init data structure.... 
placed = []
blockIdToOrder = []
for i in range(0, len(blocks)):
    placed.append(False) #everyone begins unplaced 
    blockIdToOrder.append(-1) #everyone has no position 

## start with a 'seed' species, which sets up the order 
consensusOrder = [] 
orderedBlocks1 = species[0]["orderedBlocks"]
for i in range(0, len(orderedBlocks1) ):
    placed[ orderedBlocks1[i] ] = True
    blockIdToOrder[ orderedBlocks1[i]] = i
    consensusOrder.append(orderedBlocks1[i])

## loop through all other species, not the reference 
for s in range(1, len(species)):
    orderedBlocks_s = species[s]["orderedBlocks"]  #get its ordered blocks 

    for i in range(1, len(orderedBlocks_s)): #loop through all of its blocks  
                                            #... note, later may need a 'while not done' statement
        thisBlock = orderedBlocks_s[i] 

        ## I tried to remove small blocks here, but it caused problems later.. can try something else... 
        # isTooSmall = False
        # if blocks[thisBlock]["maxLength"] < 10:
        #     isTooSmall = True
        if not placed[thisBlock]: #and not isTooSmall: #removed this 
            parent = orderedBlocks_s[i-1] #parent is behind
            child = orderedBlocks_s[i+1] #child is in front

            reorder = False 
            spliceIdx = len(consensusOrder) #find a place to 'splice' it in

            if placed[parent] and placed[child]: #are both parent and child placed? 
                if blockIdToOrder[parent] < blockIdToOrder[child]: #and in the right order 
                    spliceIdx = blockIdToOrder[parent]+1
                    consensusOrder.insert( spliceIdx, thisBlock ) #if so, put it in front of parent 
                    reorder = True
                elif blockIdToOrder[parent] > blockIdToOrder[child]: #if not the right order, invert
                    spliceIdx = blockIdToOrder[parent]  #place in front of parent 
                    consensusOrder.insert( spliceIdx, thisBlock ) 
                    reorder = True             
            elif placed[parent]: #if just parent placed 
                spliceIdx = blockIdToOrder[parent]+1 #put in front of parent (can't judge inversions now)
                consensusOrder.insert( spliceIdx, thisBlock ) 
                reorder = True 
            elif placed[child]: #if just child placed, put behind child
                spliceIdx = blockIdToOrder[child]  #place  
                consensusOrder.insert( spliceIdx, thisBlock )
                reorder = True      

                #maybe later could look for a downstream or upstream block... 

            if reorder:
                for i in range(0, len(consensusOrder)):
                    blockId = consensusOrder[i]
                    placed[blockId] = True
                    blockIdToOrder[blockId] = i 
            ## note: no else.  so, I'm not sure how everyone is getting placed.... 

### compute a linear 'length' of all the blocks together 
maxNT = 0
for i in range(0, len(consensusOrder)):
    maxNT = maxNT + blocks[ consensusOrder[i] ]["maxLength"]


lengthOfCoordinateSystem = maxNT

### DEBUGGING 
# print "length of coordinate system = %d" % lengthOfCoordinateSystem

#---------------------------------------------------------------------------------------------
#-------COMPRESSED COORDINATE GENERATION -----------------------------------------------------
#---------------------------------------------------------------------------------------------
### the next step is to 'stack' the blocks, so instead of a long list, it is more like a DAG 


lastNeighbor = [0]*len(species) ##lookup the 'most recent' neighbor for a species , running tally
lastPosition = [0]*len(species) ##lookup where this most recent neighbor is 
countNeighbors = [0]*len(blocks) #a count of neighbors- am I using this? 
lastPositionFixedWidth = [0]*len(blocks) ## i think i am not using this... 

for i in range(0, len(consensusOrder) ):

    b = blocks[ consensusOrder[i] ]
    segs = b["segList"]

    greatestLastPosition = 0;
    greatestLastPositionIdx = 0;

    for s in range(0, len(segs) ):
        sp = segments[ segs[s]["id"] ]["speciesId"]
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

    # print " block: %d has max length %d and starts at: %d next to: %d" % ( consensusOrder[i] ,  blocks[ consensusOrder[i] ]["maxLength"] , greatestLastPosition , lastNeighbor[greatestLastPositionIdx] )

    blocks[consensusOrder[i]]["compressedStart"] = greatestLastPosition;
    blocks[consensusOrder[i]]["compressedNextTo"] = lastNeighbor[greatestLastPositionIdx]
    blocks[consensusOrder[i]]["neighborOrder"] = countNeighbors[ lastNeighbor[greatestLastPositionIdx] ]

    
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


## compute how long the compressed coordinate system is, 
compressedMaxLength = 0; 
for s in range(0, len(species)):
    if lastPosition[s] > compressedMaxLength:
        compressedMaxLength = lastPosition[s]



#---------------------------------------------------------------------------------------------
#-------ASSIGN BLOCKS TO A HORIZONTAL SLOT ---------------------------------------------------
#---------------------------------------------------------------------------------------------
## Position blocks into horizontal slots- like the algorithm from the paper Design considerations for optimizing storylines
##### this is needed to find the gap blocks and initially position segments
#### later we will need to instead 'read in' sets of randomly generated slots for each block and assess whether it violates or passes

## initialize slots 
slots = []
for i in range(0, len(species)+2):
    slots.append([])

#assign all blocks to a slot 
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

        ## this is to ensure that blocks are centered in the middle, not the top 
        if j % 2 == 0:
            jIdx = int( math.floor( len(slots)/2.0)+math.floor(j/2.0) )
        else:
            jIdx = int( math.floor(len(slots)/2.0)-math.ceil(j/2.0) )                  

        #grab blocks already assigned at this index     
        blocksAssigned = slots[jIdx]
        
        #see if there are overlaps with any blocks in this slot 
        overlaps = False
        for k in range(0, len(blocksAssigned) ):
            other = blocks[ blocksAssigned[k] ]
            otherStart = other["compressedStart"]
            otherStop = other["compressedStart"] + other["maxLength"]

            if bStart <= otherStart and bStop >= otherStart:
                overlaps = True
            if bStart >= otherStart and bStart <= otherStop:
                overlaps = True


        #if not, assign it 
        if not overlaps and not found:
            assignedSlot = jIdx

            blocksAssigned.append( consensusOrder[i] )
            found = True
            blocks[consensusOrder[i]]["slot"] = assignedSlot
            slots[jIdx] = blocksAssigned

        #clear array for next loop 
        blocksAssigned = [] 

#done assigning slots


#---------------------------------------------------------------------------------------------
#-------CLASSIFY SEGMENTS WITHIN THE BLOCK AND MAKE GAP BLOCKS -------------------------------
#---------------------------------------------------------------------------------------------

#here is where each segment finds its previous neighbor, and uses this to classify itself
# as coming from above or below or being static
# this is needed for organizing lines within a block

previous = consensusOrder[1]  ## skip the 0th, it has no segments 
previousSegs = blocks[ previous ]["segList"]
blockCount = len(blocks)
segmentCount = len(segments)
done = False 
i = 0
gapBlocks = {}
gapBlocks["null"] = {}

while not done:
    idx = i
    #classify segments in all blocks
    b = blocks[ consensusOrder[idx] ]
    segs = b["segList"]
    neighborsBehindList = b["neighborsBehindList"]

    for s in range(0, len(segs)):  #for all segs 
        latestNeighbor = -1
        latestNeighborPosition = 0
        latestNeighborSlot = -1 
        latestNeighborIdx = 0
        sp = segments[segs[s]["id"]]["speciesId"]
        foundLatestNeighbor = False
        inSlot = -1

        neighborsBehindList = b["neighborsBehindList"]

        for n in range(0, len(neighborsBehindList)): #look at previous neighbors 

            otherSegs = blocks[ neighborsBehindList[n] ] ["segList"]
            for os in range(0, len(otherSegs)):
                spOs = segments[ otherSegs[os]["id"] ]["speciesId"]


                if spOs == sp:
                    if blocks[ neighborsBehindList[n] ]["compressedStart"] + blocks[neighborsBehindList[n]]["maxLength"] > latestNeighborPosition:
                        latestNeighborPosition =  blocks[ neighborsBehindList[n] ]["compressedStart"] + blocks[neighborsBehindList[n]]["maxLength"]
                        latestNeighborSlot = blocks[ neighborsBehindList[n] ]["slot"]#inSlot
                        latestNeighbor = blocks[ neighborsBehindList[n] ]
                        foundLatestNeighbor = True
                        latestNeighborIdx = os
                        # print "-------------IS LATEST-----------"
        #done looping through back neighbors 


        ## CREATE GAP BLOCK if latest neighbor is far away 
        if latestNeighborPosition != b["compressedStart"]:

            strNeigh = "%d"%latestNeighbor["id"]

            gapSize = b["compressedStart"] - latestNeighborPosition 
            # if gapSize > 10:
            #     needNew = True  ## not doing this right now.... 
            needNew = True 
            if strNeigh in gapBlocks: #have we already created a gap block for this prev neighbor?              
                neighborGapSize = gapBlocks[strNeigh]["maxLength"] #maybe, but need to see if it covers our gap

                distance = gapSize - neighborGapSize  
                if distance > 10: #the gap that was found wasn't big enough to cover our gap
                    needNew = True #so we need a new block still 

                    #but, we should alter the existing gap block so that the future gap block will be in front 
                    gapBlocks[strNeigh]["neighborInFrontList"].append( blockCount )

                    #add a segment to the short gapBlock
                    segId = segmentCount  
                    block = gapBlocks[strNeigh]["id"] 
                    speciesId = sp
                    start = -1 #just for gap segments
                    size = gapBlocks[strNeigh]["maxLength"] #it'll fill this gap block.. 
                    reverseCompliment = False 
                    contigId = -1 #just for gap segments
                    prevIdx = latestNeighborIdx #the segment is in front of the latest neighbor
                    prevSlot = latestNeighborSlot #which has this slot 
                    theClass = "static" #classify the segment, to avoid line crosses later 
                    if latestNeighborSlot < gapBlocks[strNeigh]["slot"]:
                        theClass = "descending"
                    if latestNeighborSlot > gapBlocks[strNeigh]["slot"]:
                        theClass = "rising"

                    #create new seg and add it to gap block and segment list 
                    newSeg = {"id": segId, "block": block, "contigId": contigId, "speciesId": speciesId, "start": start, "size": size, "reverseCompliment": reverseCompliment}
                    gapBlocks[strNeigh]["segList"].append( {"id": newSeg["id"], "prevIdx": prevIdx, "prevSlot": prevSlot, "latestNeighborPosition": latestNeighborPosition, "latestNeighborBlock": latestNeighbor["id"], "class": theClass} )
                    segments.append( newSeg ) 
                    segmentCount = segmentCount+1
                    gapBlocks[strNeigh]["numberAligned"] = gapBlock["numberAligned"] +1 

                    #then update the latest neighbor to point to this gap block
                    #so that the next gap added- which will fill the long gap, consideres this gap to be 'behind'
                    latestNeighbor = gapBlocks[ strNeigh ] 
                    latestNeighborPosition = gapBlocks[strNeigh]["compressedStart"] + gapBlocks[strNeigh]["maxLength"]
                    latestNeighborSlot = gapBlocks[strNeigh]["slot"]
                    latestNeighborIdx = len(gapBlocks[strNeigh]["segList"])-1

                else:
                    needNew = False #otherwise, we don't need to make a new block, can just retrieve it 
            #if we need a new block- either because gap block wasn't big enough or no gap block as been created 
            if needNew: 
            
                strNeigh = "%d"%latestNeighbor["id"] #reset just in case 

                #create gap block  
                id = blockCount
                score = 0
                numberAligned = 0 # 1 for now
                newSegList = [] #empty for now
                maxLength = b["compressedStart"] - latestNeighborPosition - 4 #distance of the gap
                compressedNextTo = latestNeighbor["id"] #should be in front of latest neighbor
                compressedStart = latestNeighborPosition + 2
                neighborOrder = 0 #am I still using this? 
                neighborFrontList = [ b["id"] ] #this block is in front of it
                isGapBlock = True


                ### find a slot for this gap block
                #### same method as before 
                assignedSlot = -1
                j = 0
                blocksAssigned = [] 
                foundSlot = False
                for j in range(0, len(slots)):

                    if j % 2 == 0:
                        jIdx = int( math.floor( len(slots)/2.0)+math.floor(j/2.0) )
                    else:
                        jIdx = int( math.floor(len(slots)/2.0)-math.ceil(j/2.0) )                  

                    blocksAssigned = slots[jIdx]
                    bStart = compressedStart
                    bStop = compressedStart + maxLength
                    
                    overlaps = False
                    for k in range(0, len(blocksAssigned) ):

                        other = blocks[ blocksAssigned[k] ]
                        otherStart = other["compressedStart"]
                        otherStop = other["compressedStart"] + other["maxLength"]
                       
                        if bStart <= otherStart and bStop >= otherStart:
                            overlaps = True
                        if bStart >= otherStart and bStart <= otherStop:
                            overlaps = True


                    if not overlaps and not foundSlot:
    
                        assignedSlot = jIdx                   
                        foundSlot = True

                    blocksAssigned = [] 

                slot = assignedSlot
                neighborsBehindList = [] #cleared for now 
                
                #add this new gap block to the running array of gap blocks 
                gapBlocks[strNeigh] = {"id": id, 'score': score, 'numberAligned': numberAligned, 'segList': newSegList, 'maxLength': maxLength,'compressedNextTo': compressedNextTo, "compressedStart": compressedStart, "neighborOrder": neighborOrder, "isGapBlock": isGapBlock, "neighborInFrontList": neighborFrontList, "slot": slot, "neighborsBehindList": [] }
                blockCount = blockCount +1 

                consensusOrder.insert(i, id) #insert it into the consensus order
                i = i + 1  #skip   
                idx = i #move idx forward

                #add this gap block to the blocks array
                blocks.append(gapBlocks[strNeigh])    
                b["neighborsBehindList"].append(id)  ## note, this is happenning too many times for some reason
                slots[slot].append(id) #add this to the slots array 

                ### DONE creating new GAP BLOCK

            #now for either new gap blocks or retrieved ones, grab it
            gapBlock = gapBlocks[strNeigh]


            #add a segment to the gapBlock
            segId = segmentCount 
            block = gapBlock["id"]
            speciesId = sp
            start = -1 #just for gap segments
            size = gapBlock["maxLength"]
            reverseCompliment = False
            contigId = -1 #just for gap segments
            prevIdx = latestNeighborIdx
            prevSlot = latestNeighborSlot
            theClass = "static"
            if latestNeighborSlot < gapBlock["slot"]:
                theClass = "descending"
            if latestNeighborSlot > gapBlock["slot"]:
                theClass = "rising"
  
            newSeg = {"id": segId, "block": block, "contigId": contigId, "speciesId": speciesId, "start": start, "size": size, "reverseCompliment": reverseCompliment}
            gapBlock["segList"].append( {"id": newSeg["id"], "prevIdx": prevIdx, "prevSlot": prevSlot, "latestNeighborPosition": latestNeighborPosition, "latestNeighborBlock": latestNeighbor["id"], "class": theClass} )
            segments.append( newSeg )
            segmentCount = segmentCount+1
            gapBlock["numberAligned"] = gapBlock["numberAligned"] +1

            gapBlock["neighborsBehindList"].append( latestNeighbor["id"] )

            # reset the neighbors on the segment
            latestNeighborIdx = len(gapBlock["segList"])-1
            latestNeighborSlot = gapBlock["slot"]
            latestNeighborPosition = gapBlock["compressedStart"] + gapBlock["maxLength"]
            latestNeighbor = gapBlock


        #done adding segment to the gap block 

        ## not classify the segment in the current block (non-gap )
        if foundLatestNeighbor:

            if latestNeighborSlot == b["slot"]:
                b["segList"][s]["class"] = "static";

            if latestNeighborSlot < b["slot"]:
                b["segList"][s]["class"] = "descending"#"rising";

            if latestNeighborSlot > b["slot"]:
                b["segList"][s]["class"] = "rising"#"descending"
  
            b["segList"][s]["prevIdx"]  = latestNeighborIdx
            b["segList"][s]["prevSlot"] = latestNeighborSlot
            b["segList"][s]["latestNeighborPosition"] = latestNeighborPosition
            b["segList"][s]["latestNeighborBlock"] = latestNeighbor["id"]



    #update the blocks list -- i think python is by ref so prob not necessary.. 
    blocks[consensusOrder[idx]] = b
           
    i = i + 1 #move i forward 
    if i == len(consensusOrder):  #are we done? 
        done = True


#---------------------------------------------------------------------------------------------
#-------ORDER SEGMENTS TO MINIMIZE LINE CROSSINGS --------------------------------------------
#---------------------------------------------------------------------------------------------
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
        otherSegs = neighborBlock["segList"]
        for os in range(0, len(otherSegs)):
            spOs = segments[ otherSegs[os]["id"] ]["speciesId"]

            
            if spOs == sp:
                b["segList"][s]["prevIdx"] = os


        for rs in range(0, len(reorderedSegs)):
            if not found:
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


#---------------------------------------------------------------------------------------------
#-------SHIFT BLOCKS IN A SLOT---------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


#first assign blocks to an initial y coordinate- based on slot
for i in range(0, len(consensusOrder)):
    blocks[consensusOrder[i]]["yAssigned"] = blocks[consensusOrder[i]]["slot"]*len(species)*10 +10
    blocks[consensusOrder[i]]["numberAligned"] = len( blocks[consensusOrder[i]]["segList"])
    # print blocks[i]["yAssigned"]


#shift y coordinate to line up segments within slots 
#based on the most common length of 'shift'
#eg, if everyone would be straight if shifted by 10, expect for one, still shift by 10
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

    shiftCount = {}
    for j in range(0, len(yShift)):
        lookup = "%d" % yShift[j]
        if lookup in shiftCount:
            shiftCount[lookup] = shiftCount[lookup] +1 
        else :
            shiftCount[lookup] = 1
   
    mostCommonCount = 0
    mostCommonValue = 0
    for key in shiftCount:
        if shiftCount[key] > mostCommonCount:
            mostCommonCount = shiftCount[key]
            mostCommonValue = int(key)

    # print "maxShift: %d and consensusShift %d " %( maxShift, consensusShift)
    b["yAssigned"] =  b["yAssigned"] + mostCommonValue*10


#---------------------------------------------------------------------------------------------
#-------REMOVE WHITESPACE---------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

#find the slot segments
#which are contiguous sets of blocks in a slot that share some number of lines 
#these should be moved up or down as a unit 
slotSegments = []

for i in range(0, len(slots)):
    blocksAssigned = slots[i]

    slotSegments.append({"slot": i, "count": 0, "slotSegs": [] })

    previous = blocksAssigned[0]
    runningSlotSegment = {"start": blocks[previous]["compressedStart"], "stop": blocks[previous]["compressedStart"]+ blocks[previous]["maxLength"], "segs": [] }
    runningSlotSegment["segs"].append(previous)
    #runningSlotSegment.append(previous)
    for j in range(1, len(blocksAssigned)):
        current = blocksAssigned[j]

        isContiguous = False
        if blocks[current]["compressedStart"] - (blocks[previous]["compressedStart"] + blocks[previous]["maxLength"]) < 10:
            isContiguous = True 

        isASlotSegment = False
        if isContiguous:
            prevSegs = blocks[previous]["segList"]
            currSegs = blocks[current]["segList"]

            for k in range(0, len(prevSegs)):
                for l in range(0, len(currSegs)):
                    if segments[ prevSegs[k]["id"] ]["speciesId"] == segments[ currSegs[l]["id"] ]["speciesId"]:
                        isASlotSegment = True
                        break

        if isASlotSegment and isContiguous:
            runningSlotSegment["segs"].append(previous)
            runningSlotSegment["stop"] = runningSlotSegment["stop"] + blocks[previous]["maxLength"]
            # print blocks[previous]["compressedStart"]
        else:
            slotSegments[i]["slotSegs"].append(runningSlotSegment)
            slotSegments[i]["count"] = slotSegments[i]["count"] + 1
            runningSlotSegment = {"start": blocks[current]["compressedStart"], "stop": blocks[current]["compressedStart"]+ blocks[current]["maxLength"], "segs": [] }
            runningSlotSegment["segs"].append(current)

        previous = current


# for i in range(0, len(slotSegments)):
#     print "------ %d" % i
#     print slotSegments[i]

## now that you have slot segments, compute a running tally of the floor, based on the floor of other segments
## to do:


#---------------------------------------------------------------------------------------------
#-------EVALUATE------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

## evaluation of the performance of this set of slots is based on:
### for all species: sum the crossovers and wiggles 
### crossovers:  how many lines does one line cross in going from block i to block j?
### how to compute: 





#---------------------------------------------------------------------------------------------
#-------OUTPUT RESULT-------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------


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


        
