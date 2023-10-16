
#Import modules
import numpy as np
import matplotlib.pyplot as plt
import sys

#Physical paramters
m = 2 #Mass of particle
k = 1 #Spring constant
n = 2 #Number of masses
#Set up matrix

M = np.zeros((n,n))
np.fill_diagonal(M,-2*k/m)
for i in range(1, n):
        M[i, i - 1] = k/m
        M[i - 1, i] = k/m

#M = np.matrix([[-2*k/m, k/m] , [k/m, -2*k/m]])
print(M)
#Set up matrix we need to use the QU method on
N=M
df1 =2
counter = 0

while df1>0.0001:
    N1 = N[0,0]
    N2 = N[0,1]
    N3 = N[1,0]
    N4 = N[1,1]
    #print(N1,N2,N3,N4)

    c1 = np.array([N1,N3])
    f1 = c1
    mag_f1 = np.sqrt(f1[0]**2+f1[1]**2)

    c2 = np.array([N2,N4])
    f2 = c2-(np.multiply(np.dot(c2,f1),f1))/((mag_f1)**2)
    mag_f2 = np.sqrt(f2[0]**2+f2[1]**2)

    q1 = f1/(mag_f1)
    q2 = f2/(mag_f2)
    Q1 = np.matrix([q1,q2]).reshape(2,2).transpose()
    U11 = mag_f1
    U12 = np.dot(c2,q1)
    U13 = 0
    U14 = mag_f2
    U1 = np.matrix([[U11,U12],[U13,U14]])
    #print(U1)
    #print(Q1)
    N = np.matmul(U1,Q1)
    print(N)
    previous_f1 = N1
    previous_f2 = N4
    current_f1 = N[0,0]
    current_f2 = N[1,1]
    #print(previous_f1,previous_f2)
    #print(current_f1,current_f2)
    df1 = previous_f1-current_f1
    df2 = previous_f2-current_f2
    print(df1,df2)
    counter = counter+1
    print(counter)
    #print(f"f1={f1},f2={f2},|f1|={mag_f1},|f2|={mag_f2}")
print(f"The eigenvalues are {current_f1},{current_f2}")