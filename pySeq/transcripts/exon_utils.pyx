def uniqufy(nonunique):
    """ Returns the unique elements of a list.  Does not maintain order.
    """
    s=set(nonunique)
    return(list(s))

def overlap_cases(int1, int2):
    """Returns the number of overlapping base pairs between two intervals
    """
    if int1[0] > int2[0] and int1[1] < int2[1]:
        # Case where int1 is wholly enclosed in interval 2
        return(int1[1]-int1[0])
    elif int1[0]>int2[0] and int1[0]<int2[1]:
        return(int2[1]-int1[0]) 
    else:
        return(0)

def collideIntervals(tuple i1, tuple i2): 
    """Spends most of the time here.
    """
    if i1[0] > i2[0] and i1[1] < i2[1]:
        # Case where int1 is wholly enclosed in interval 2
        return([i2])
    elif i2[0] > i1[0] and i2[1] < i1[1]:
        # Case where i2 is wholly enclosed in i1
        return([i1])
    elif i1[0] > i2[0] and i1[0] < i2[1]:
        return([(i2[0],i1[1])]) 
    elif i2[0] > i1[0] and i2[0] < i1[1]:
        return([(i1[0], i2[1])])
    else:
        return([i1, i2])

def permutation_indices(data):
    return sorted(range(len(data)), key=data.__getitem__)

def swapIntervals(intervalList, i, j):
    temp = intervalList[i]
    intervalList[i] = intervalList[j]
    intervalList[j] = temp

def collapseIntervals(intervalList):
    """ Runs collide interval across all intervals in a list, returning a list
    of only the intervals with no overlaps with any other intervals.
    This is O(n^2)
    """
    cdef int size, index = 0

    size = len(intervalList)

    while index < size - 1:
        t = collideIntervals(intervalList[index], intervalList[index+1])
        if  len(t) == 1:
            intervalList.pop(index)
            intervalList.pop(index)
            intervalList.insert(index, t[0])
            size -= 1
            index = 0
        elif len(t) == 2:
            if intervalList[index][0] > intervalList[index+1][0]:
                swapIntervals(intervalList, index, index+1)
                index = 0
            else:
                index += 1
        else:
            index += 1
    return(intervalList)



def calcOverlap(intervals):
    """
    Calculates the number of base pairs that overlap in a list of intervals

    All intervals must be unique
    :TODO Full enclosed case fails.  Probably should just use interval trees,
    but constructing trees might take longer, if matches are sparse? or would
    it be faster? Just do it.

    """

    bp = 0
    for i in intervals:
        bp += sum([overlap_cases(i, j) for j in intervals])
    return(bp)

