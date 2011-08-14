"""
Parses Japser Matrices
"""
#import optparse
#p = optparse.OptionParser()
#p.add_option("-r", "--restrict", action="store", dest="r", help="search with only one PWM", default='')
#options, args = p.parse_args()
def parse(pwm):
    name = ''
    pwmsize = 0
    A = []
    C = []
    G = []
    T = []
    matrix = []
    matricies = {} 
    sizes = {} 
    index = {}

    for line in pwm:
	if line[0] == '>':
	    name = line.strip().split(' ')[1]
	#Parses PWMs
	elif line[0] == 'A':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    A = map(float, line)
	    pwmsize = len(A)
	elif line[0] == 'C':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    C = map(float, line)    
	elif line[0] == 'G':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    G = map(float, line)
	elif line[0] == 'T':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    T = map(float, line)
	    matrix = [A,C,G,T]
	    matricies[name] = matrix
	    sizes[name] = pwmsize
	    index[name] = len(matricies) - 1
    return (index, matricies, sizes)
