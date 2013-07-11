# --------------------------------------------------------------------------------------------
# Kumar Thurimella
# 2012 May 24
#
# Run file:
# python processing_peaks.py <Filtered File> <Simulated DNA Text File> Probability Threshold
#
# This program processes the filtered data read in from a text file
# from the Filter Bowties pipeline. It then analyzes the data and
# finds statistically significant regions throughout the genome.
# 
# --------------------------------------------------------------------------------------------

from time import clock
from collections import defaultdict
from scipy.stats import poisson
import numpy as np
import sys, os
import itertools
import random

# Raise error if the correct number of arguments isn't correct
if (len(sys.argv) != 4):
    sys.exit("\n\nUsage: processing_peaks.py <Filtered Text File> <Simulated DNA Text File> Probability Threshold")

# Output the data in a text file using .write function    
annotationFile = os.path.splitext(sys.argv[1])[0]
out = open(annotationFile + "_processed.txt", 'w')
out.write("\t\t\tGenome\tSignificant Peak Window\tLength of Genome Covered")

# Calculate a poisson distribution (function)
def poissonTest(enrichmentList, controlList):
    return poisson.pmf(enrichmentList, controlList)

# Saturate every control (expected value) enrichment value to 1
# just to make sure the poisson mass function does not go to zero
# as it will if the expected value is 0. So whether or not the experimental
# value is significant or not, it will be read as zero and will skew our
# likeliness generator.  
def addOne(x):
    return x + 1

def filterControl(controlList, normalLength, controlLength):
    newList = []
    
    count = 0
    
    # Randomly pick fragments from the given control fragments. 
    # The number of fragments to pick is the number of experimental
    # fragments.
    while count < normalLength:
        x = random.randint(0, controlLength - 1)
        newList.append(controlList[x])        
        count = count + 1
    
    return newList

def enrichmentRegion(dictList, givenLength):
    startList = []
    stopList = []
    
    
    # Separate start and stop positions
    for i in xrange(givenLength):
        if i % 2 == 0:
            startList.append(dictList[i])
        else:
            stopList.append(dictList[i])
            
    # Sort the whole key's list         
    dictList.sort()

    # Initialize
    startTuple = []
    stopTuple = []
    
    
    # This step essentially keeps track of all the 
    # coordinates and labels them '0' if it is a start
    # coordinate and '1' if it is a stop coordinate.
    for j in xrange(len(startList)):    
        startTupleTemp = (startList[j], 0)
        stopTupleTemp = (stopList[j], 1)
        
        startTuple.append(startTupleTemp)
        stopTuple.append(stopTupleTemp)
    
    # Create a tuple using the addition operator of
    # the start and stop list i.e. merge the two 
    # lists.         
    startStopList = startTuple + stopTuple
    
    # This sorts the whole start and stop list based
    # on coordinate. Assuming that any start will come before
    # its respective stop this gives us coverage over that
    # specific coordinate. 
    startStopList.sort()
    
    # Initialize    
    counter = 0
    peakList = []

    # This keeps track of the coordinate and if it is a
    # start coordinate then add a level of coverage 
    # (counter += 1) and take away a level if we reach a 
    # stop. 
    for (i, j) in startStopList:
        if j == 0:
            counter += 1
            peakList.append((counter, i))
        else:
            counter -= 1
            peakList.append((counter, i))
    
    # Define a genome length by looking at the very last 
    # position of the given genome's fragments.       
    j = len(peakList) - 1    
    (_, b) = peakList[j]
    
    # Pre-allocate an array of zero's based on genome size  
    enrichmentList = [0] * b
    
    # Create a dictionary to essentially parse through 
    # a tuple and capture what the enrichment value was
    # at a given position on a genome
    peakDict = {}
    
    peakList.sort(key=lambda tup: tup[1])
    
    for (x,y) in peakList:
        peakDict[y] = x
        
    temp = sorted(peakDict.iterkeys())
    
    # This step combines the idea from above but specifically
    # fills in every position of the genome by looking at the 
    # gaps and reads the value from left to right.    
    for i in xrange(0, len(temp)):
        if (i+1) == len(temp):
            for j in xrange(temp[i]-1, len(temp)):
                enrichmentList[j] = peakDict[temp[i]]
        else:
            for j in xrange(temp[i]-1, temp[i+1]):
                enrichmentList[j] = peakDict[temp[i]]
    
    return enrichmentList
    

# Initialize 
startStopList = []
genomeList = []

# Read more on input/output streams on the python tutorial
# and read in specific lines from a text file. 
f = open(sys.argv[1], 'r')
g = open(sys.argv[2], 'r')

filter_peaks = f.readlines()
filter_peaks_control = g.readlines()

for x in filter_peaks:
    line = x.strip("\n").split()
    # Import all of the lines into python and manipulate the data from the 
    # filtered text file. The start and stop portions of the data are read
    # in as integer types while the genome is read in as a string.
    stop = int(line[-1])
    start = int(line[-2])
    genome = str(line[2])
    
    # Create a tuple from the start position, stop position and the
    # organism (genome) name. Append everything into a list.
    startStopTuple = (start, stop, genome)
    startStopList.append(startStopTuple)
    genomeList.append(genome)

# Initialize 
startStopListControl = []
genomeListControl = []
startStopListControlOld = []

for x in filter_peaks_control:
    line = x.strip("\n").split()
    # Import all of the lines into python and manipulate the data from the 
    # filtered text file. The start and stop portions of the data are read
    # in as integer types while the genome is read in as a string.
    stopControl = int(line[-1])
    startControl = int(line[-2])
    genomeControl = str(line[2])
    
    # Create a tuple from the start position, stop position and the
    # organism (genome) name. Append everything into a list for the control.
    startStopTupleControl = (startControl, stopControl, genomeControl)
    startStopListControlOld.append(startStopTupleControl)
    genomeListControl.append(genomeControl)

# Sort the list. 
startStopList.sort()
startStopListControlOld.sort()

# Keep track of every organism that is in the list. 
# The set function will keep all of the unique ones.
genomeList = list(set(genomeList))
genomeList.sort()

# Find all the fragments in the control files of the
# same organism that is in the experimental. 
for i in genomeList:
    for x, y, z in startStopListControlOld:
        if i == z:
            startStopListControl.append((x, y, z))

# Sort each the list of tuples and the control list of 
# tuples by the organism name.
startStopList.sort(key = lambda tup: tup[2])
startStopListControl.sort(key = lambda tup: tup[2])

# Group all the tuples by the organism name and order into
# separate lists. 
splitListControl = [list(g) for key, g in itertools.groupby(startStopListControl, lambda tup: tup[2])]
splitList = [list(cluster) for key, cluster in itertools.groupby(startStopList, lambda tup: tup[2])]

# Initialize 
lengthList = []
lengthListControl = []

# Consolidate the number of fragments in the control and
# experimental. 
for i in xrange(len(splitList)):
    lengthList.append(len(splitList[i]))
    lengthListControl.append(len(splitListControl[i]))

normalizedLists = []

# This step normalizes the control set based on the number of fragments
# in the control and experimental set. If the experimental set has more
# fragments than the control set the user is notified and is alerted to
# choose a stringent p-value. The control set is now modified.
for j in xrange(len(splitList)):
    if lengthList[j] < lengthListControl[j]:
        normalizedList = filterControl(splitListControl[j], lengthList[j], lengthListControl[j])
        normalizedLists.append(normalizedList)
    else:
        print 'Please make sure to choose a stringent p-value to filter for unnecessary reads for the organism ' + genomeList[j] + '.'
        normalizedLists.append(splitListControl[j])

startStopListControl = [item for sublist in normalizedLists for item in sublist]

# Use the dictionary data structure and organize the data by keys.
# Each key corresponds to a given organism and each key generates the
# start and stop position for that organism.
genomeDict= defaultdict(list)
for n,v,k in startStopList:
        genomeDict[k].append(n)
        genomeDict[k].append(v)



# Initialize 
genomeDictList = []
listOfLists = []

# Iterate through the keys of the text file
# which are all based on the organism name
# and from there do further computations based
# on the specific data.
for key in sorted(genomeDict.iterkeys()):
    genomeDictList = genomeDict[key]
    length = len(genomeDictList)
        
    normalEnrichmentList = enrichmentRegion(genomeDictList, length)
    
    listOfLists.append(normalEnrichmentList)
    
    genomeDictList[:] = []


# Use the dictionary data structure and organize the data by keys.
# Each key corresponds to a given organism and each key generates the
# start and stop position for that organism.
genomeDictControl= defaultdict(list)
for t,u,s in startStopListControl:
        genomeDictControl[s].append(t)
        genomeDictControl[s].append(u)


# Initialize 
genomeDictListControl = []
listOfListsControl = []

# Iterate through the keys of the control text file
# which are all based on the organism name
# and from there do further computations based 
# on the specific data. Note that functional 
# programming tool 'map' is used to saturate the
# control set. 
for keys in sorted(genomeDictControl.iterkeys()):
    genomeDictListControl = genomeDictControl[keys]
    lengthControl = len(genomeDictListControl)
        
    controlEnrichmentList = enrichmentRegion(genomeDictListControl, lengthControl)
    
    newControlEnrichmentList = map(addOne, controlEnrichmentList)
            
    listOfListsControl.append(newControlEnrichmentList)
    
    genomeDictListControl[:] = []
        
# Initialize 
length = []
fragmentControl = []
controlDistributionList = []
controlDistributionLists = []

for i in xrange(len(listOfLists)):
    length.append(len(listOfLists[i]))
    
for i in xrange(len(listOfListsControl)):
    fragmentControl.append(len(listOfListsControl[i]))

# This part accounts for any disparities in the genomes
# between the control set and the normal set and fixes 
# any discrepency. 
for i in xrange(len(length)):
    current = listOfListsControl[i]
    
    if len(current) < len(listOfLists[i]):
        listOfListsControl[i] = current + [0]*(length[i] - fragmentControl[i])
    elif len(current) > len(listOfLists[i]):
        listOfListsControl[i] = current[:length[i]]

# This step computes the probabilities based on the Poisson
# Mass Function by using the chance of occurrence based on
# the control list. The minimum p-value is printed for the
# user to choose a more stringent filter process. This
# is the most computationally taxing step as all the data
# has to be converted to numpy objects and then back to native
# python objects.
for i in xrange(len(listOfLists)):  
    numpyProbabilityList = poissonTest(listOfLists[i], listOfListsControl[i])
    whereAreNaN = np.isnan(numpyProbabilityList)
    numpyProbabilityList[whereAreNaN] = 0
    
    minval = np.min(numpyProbabilityList[np.nonzero(numpyProbabilityList)])
    print "The minimum p-value for " + genomeList[i] + " is " + str(minval) + "."
    
    probabilityList = numpyProbabilityList.tolist()
    
    genomeName = genomeList[i]
    lengthOfGenome = length[i]
    
    # This gives the windows of significant peaks based on the
    # the user input
    windowList = []
        
    for i in xrange(len(probabilityList)):
        if float(sys.argv[3]) >= probabilityList[i]:
            windowList.append(i + 1)
    
    rangeList = []
    
    # This appends all the beginning and ends of the window
    # which gives ranges instead of every individual base pair that 
    # meets the threshold. 
    if len(windowList) != 0:        
        rangeList.append(windowList[0])
        
        for i in xrange(len(windowList) - 1):
            if windowList[i] != (windowList[i+1] - 1):
                rangeList.append(windowList[i])
                rangeList.append(windowList[i+1])
                
        rangeList.append(windowList[len(windowList) - 1])
    
        
    # Write the Peaks, the Poisson probability and the size of 
    # the whole list into a text file.
    out.write("\n" + str(genomeName))     
    out.write("\n" + str(rangeList))
    out.write("\n" + str(lengthOfGenome))
    
# Check run time.
print 'The run-time is: ' + str(clock()) + ' seconds'