import sys

def regionParse(string):
    """
    Parses region:start_position-end_position, returns a tuple
    Example: 
    >>> regionParse('chr4:200000-200100') 
    ('chr4', 200000, 200100)
    """
    start = string[string.find(':')+1:string.find('-')]
    end = string[string.find('-')+1:]
    start = start.replace(',','')
    end = end.replace(',', '')

    try:
	start = int(start)
	end = int(start)
    except ValueError:
	sys.exit()
    return(string[0:string.find(':')] , start, end)
    
