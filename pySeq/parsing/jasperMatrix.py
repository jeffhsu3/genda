"""
Parses jasper matrixes
"""

def main():
    import sys
    pwm = open(sys.argv[1], "r")
    first_matrix = True 
    name_matrix = {}
    matrix_list = []
    matrix_c = 0
    for line in pwm:
	if line[0] == '>':
	    if first_matrix:
		name = line.strip().split(' ')[1]
		first_matrx = False
	    else:
		matrix_list.append[A,C,G,T]		
		name_matrix[name] = matrix_c
		name = line.strip().split(' ')[1]
		matrix_c += 1
	#Parses PWMs
	elif line[0] == 'A':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    A = line
	    pwmsize = len(A)
	elif line[0] == 'C':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    C = line    
	elif line[0] == 'G':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    G = line
	elif line[0] == 'T':
	    line = line[4:-2].split(' ')
	    while '' in line:
		line.remove('')
	    T = line
    print(matrix_list)
if __name__ == '__main__':
    main()
