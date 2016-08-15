import json


#build a block
def createBlock(idx, score, numberAligned, segList, maxLength):
    return {'id': idx, 'score': score, 'numberAligned': numberAligned, 'segList': segList, "maxLength": maxLength}

def getSpeciesId(speciesName):
    if speciesName in speciesToIdHM:
        return speciesToIdHM[speciesName]
    else:
        id = len(speciesToIdHM) 
        species.append( {"id": id, "name": speciesName, "contigList": [] } )
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

def addSegmentToContig(segment, contigId):
    #first add segment to contig
    segId = segment["id"]
    segStart = segment["start"]
    segBlock = segment["block"]
    segLength = segment["size"] 

    insertSegmentAtIdx = 0
    segArray = contigs[contigId]["segmentList"]
    for i in range( 0, len(segArray) ):
        if segStart > segments[ segArray[i] ]["start"]:
            insertSegmentAtIdx = i
            break

    segArray.insert(insertSegmentAtIdx, segId )

    contigs[contigId]["segmentList"] = segArray

    #then produce ordered list of blocks
    orderedBlocks = contigs[contigId]["orderedBlocks"]
    orderedBlocks.insert( insertSegmentAtIdx, segBlock )
    contigs[contigId]["orderedBlocks"] = orderedBlocks

    contigs[contigId]["contigLength"] = contigs[contigId]["contigLength"] + segLength

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

#------loop through file
blockNumber = 0
segmentNumber = 0

#------ Create 0th block, which is empty and represents the start
blocks.append( createBlock(0, 0, 0, {}, 0) )


for line in f: 

    lineType = line[:1]

    if lineType == "a":

        if( blockNumber > 0 ):
            blocks.append( createBlock(blockNumber, score, numberAligned, segList, maxLength) )

        blockNumber = blockNumber+1 #increment count.  Note, there is no 0th block from data
        # print blockNumber

        #segment info to capture in subsequent lines
        segList = []; #this will store the ids of the segments in the block
        maxLength = 0 #what is the length of the longest segment in this block

        #block info
        tokens = line.split(" ") #split the line
        score = int(tokens[1].split("=")[1]) #score value for block
        numberAligned = int(tokens[3].split("=")[1]) #how many segments are in the block
        


    if lineType == "s":
        #note, there is a 0th segment
        tokens = line.split()

        speciesName = int(tokens[1].split('.')[0])
        contig = tokens[1].split('.')[1]
        speciesId = getSpeciesId(speciesName)
        contigId = getContigId(contig, speciesId)

        # addToSpecies(speciesId, speciesName, contigId, segmentNumber);
        # addToContigs(contigId, contig, speciesId, segmentNumber);

        start = int(tokens[2])
        size = int(tokens[3])
        reverseComplimentStr = tokens[4]
        revComp = False
        if reverseComplimentStr.find("+") == -1:
            revComp = True
        featureLen = int(tokens[5])

        if( revComp ):
            start = featureLen - start - size; 

        if( size > maxLength ):
            maxLength = size

        segList.append( {"id": segmentNumber} )  #for the block data structure

        segments.append( {"id": segmentNumber, "speciesId": speciesId, "contigId":contigId, "start": start, "size": size, "reverseCompliment": revComp, "block": (blockNumber)} );

        #handle species and contig info
        addSegmentToContig( segments[segmentNumber], contigId )

        segmentNumber = segmentNumber+1

#add a starting block
#and an ending block
blocks.append( createBlock(blockNumber, 0, len(species), {}, 0) )
blocks[0]["numberAligned"] = len(species)

print "NUM BLOCKS %d  LAST BLOCK  %d" %( len(blocks), blockNumber )

#add contig to species, ordered by length of contig
for i in range(0, len(contigs) ):
    speciesId = contigs[i]["speciesId"]

    insertIdx = 0
    contigList = species[speciesId]["contigList"] 
    contigLen = contigs[i]["contigLength"]
    for j in range(0, len( contigList) ):
        otherLen = contigs[ contigList[j] ]["contigLength"]
        if contigLen > otherLen:
            insertIdx = j

    species[speciesId]["contigList"].insert( insertIdx, contigs[i]["id"] )


#create ordered list of blocks for the species
for i in range(0, len(species) ):
    contigList = species[i]["contigList"]
    orderedBlocks = ["0"]  
    for j in range(0, len( contigList) ):  
        orderedBlocks = orderedBlocks + contigs[ contigList[j] ]["orderedBlocks"]
    orderedBlocks.append( "%d"%blockNumber )
    species[i]["orderedBlocks"] = orderedBlocks


allJson["blocks"] = blocks
allJson["segments"] = segments
allJson["contigs"] = contigs
allJson["species"] = species

data_string = json.dumps(allJson, sort_keys=True, indent=2)
print  data_string


# print "blocks:"
# print blocks
# print "-------"
# print "segments:"
# print segments
# print "-------"
# print "species:"
# print species
# print "-------"
# print "contigs:"
# print contigs
        
