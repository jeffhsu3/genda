from PWMparser import *
t = open("/home/hsuj/Downloads/All_PWMs/SCI09/Gcm1_pwm_primary.txt", 'rU')
index, matrix, size = uniProbe_parse(t)
print(matrix)
print(size)
