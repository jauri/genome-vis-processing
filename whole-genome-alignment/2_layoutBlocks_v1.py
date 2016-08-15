import json
# import yaml

json_file = open("processStep1/ecoli.json")
json_data = json.load(json_file)

blocks = json_data["blocks"]
segments = json_data["segments"]
species = json_data["species"]
contigs = json_data["contigs"]

genomesAndAlignments = {}
maxLengths = [] 

startingPositions = []
blockSize = []
blockDepth = []
overlapsWith = []
globalSegmentDepth = []

slengthOfCoordinateSystem = 0

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
			# print "		j= %d" % j 
			# print "		arr.len= %d" % len(arr)
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
		# print "REMOVED %d, %d" % (e["from"], e["to"] )


# print "------- no cycle max weights"
# print maxNoCycle


# print "------- positioning segments"


#only works because 0th is the root 
blockDepth[0] = 0
for i in range(0, len(maxNoCycle) ):
	arr = maxNoCycle[i] #children
	l = maxLengths[i] #max length of this segment
	for j in range(0, len(arr)):
		e = arr[j]
		if( startingPositions[ e["to"] ] < maxLengths[ e["from"] ] + startingPositions[ e["from"] ] ):
			startingPositions[ e["to"] ] = maxLengths[ e["from"] ] + startingPositions[ e["from"] ]			
		if( blockDepth[ e["to"] ] < blockDepth[ e["from"] ] + 1 ):
			blockDepth[ e["to"] ] = blockDepth[ e["from"] ] + 1



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
# print "starting positions:"
# print startingPositions
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
	blocks[i]["overlapsWith"] = overlapsWith[i]
	if( blocks[i]["startingPosition"] + blocks[i]["maxLength"] > lengthOfCoordinateSystem ):
		lengthOfCoordinateSystem = blocks[i]["startingPosition"] + blocks[i]["maxLength"] 

# print "length of coordinate system = %d" % lengthOfCoordinateSystem

allJson = {}
allJson["blocks"] = blocks
allJson["segments"] = segments
allJson["contigs"] = contigs
allJson["species"] = species
allJson["lengthOfCoordinateSystem"] = lengthOfCoordinateSystem


data_string = json.dumps(allJson, sort_keys=True, indent=2)
print data_string
