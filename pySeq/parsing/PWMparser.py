"""
Parses Japser Matrices
"""
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

def uniProbe_parse(pwm):
    """ Parses a multi-matrix file from UniProbe
    """
    import re
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
    end_of_a_matrix = False
    for line in pwm:
        if end_of_a_matrix == False:
            name = line
            print(name)
            end_of_a_matrix = True
        elif line[0] == "A":
            line = line[:-1].split('\t')[1:]
            A = map(float,line)
            pwmsize = len(A)
        elif line[0] == "C":
            line = line[:-1].split('\t')[1:]
            C = map(float,line)
        elif line[0] == "G":
            line = line[:-1].split('\t')[1:]
            G = map(float,line)
        elif line[0] == "T":
            line = line[:-1].split('\t')[1:]
            T = map(float, line)
            end_of_a_matrix = True
            matrix = [A,C,G,T]
            matricies[name] = matrix
            sizes[name] = pwmsize
            index[name] = len(matricies) - 1
    return (index, matricies, sizes)
