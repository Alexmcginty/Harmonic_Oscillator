#Import modules
import numpy as np

def Eigenvalue(n,N,constant,matrix):
    """A function to find the eigenvalues for an n dimeonsiun square matrix.
    n = number of dimensions for a square matrix.
    N = number of itterations.
    constant = A switch for if you want a system of identical particles.
    matrix = if you don't want indentical particles, you can write your matrix in here."""

    if constant==True: # If you want identical particles
        #Physical paramters
        m = float(input('Please enter the mass of the particle in kg'))
        k = float(input('Please enter the spring constant between each particle in N/m'))
        if m<0:
             m = float(input('Please enter a valid mass of the particle in kg'))
        if k<0:
             k = float(input('Please enter a valid spring constant between each particle in N/m'))
        #Set up matrix, for identical particles the matrix will have diagonals of -2k/m
        # off diagonal elements will be k/m and everything else will be 0
        A = np.zeros((n,n))
        np.fill_diagonal(A,-2*k/m)
        for i in range(1, n):
                A[i, i - 1] = k/m
                A[i - 1, i] = k/m
        counter = 0 # Keep track of ittertaions 
        
        while N>counter:
            print(f'The matrix we are solving for is {A}, and the current itteration is {counter}')
            Q = np.zeros_like(A, dtype=float) # Create matrix Q
            f = np.zeros_like(A, dtype=float) # Create matrix f
            U = np.zeros_like(A, dtype=float) # Create matrix U
     
            f[:, 0] = A[:, 0]   #This will get the first column of matrix f
            Q[:, 0] = f[:, 0] / np.linalg.norm(f[:, 0]) # This will get the first column of matrix Q

            #Now for the other columns of the matrix A I loop over each column to find the elements of matrix f
            for i in range(1,n): 
                for j in range(i): 
                    f[:,i] +=  -(np.multiply(np.dot(A[:,i] , f[:, j]),f[:, j]))/(np.linalg.norm(f[:,j]))**2
                    
                f[:,i] = f[:,i]+A[:,i]  #As A[:,i] is the column vector of the ith column of the matrix A we don't want to loop j over it  
                Q[:, i] = f[:,i] / np.linalg.norm(f[:, i]) # compute each Q vector
            
            #We now find the matrix U
            for i in range(0,n):
                for j in range(0,n):
                    #The digonal elements of the matrix U is always vector f normalised
                    if i==j:
                        U[i,i]= np.linalg.norm(f[:,i])
                    #The matrix elements above the diagonal are found here    
                    if j>i:
                        U[i,j]=np.dot(A[:,j],Q[:,i])
                        
            
            UQ = np.matmul(U,Q)  #Find the matrix UQ
            A=UQ #Redfine the matrix A to the result of the UQ for use in the next itteration.
            
            Eigen =np.zeros(n) #Define an array for the eigenvalues

            #Extract the digonal elements of the matrix UQ which are the eigenvalues of the matirx A
            for i in range(0,n):
                for j in range(0,n):
                    if i==j:
                        Eigen[i]=UQ[i,j]

            print(f'The eigenvalues are {Eigen[:]}, for itteration {counter}')            
            counter=1+counter       
        
        omega=np.sqrt(-Eigen[:]) #Convert the eigenvalues to the natural frequenices
        omega1= list(omega)
        print(f'The natural frequency of oscillation for each particle is {omega1}') #Print the natural frequencies of the last itteration
    
    #If the user doesn't want identical particles we use this procudure
    else:
        A = matrix #Read the input matrix in
        counter = 0 
        while N>counter:
            print(f'The matrix we are solving for is {A}, and the current itteration is {counter}')
            Q = np.zeros_like(A, dtype=float) # Create matrix Q
            f = np.zeros_like(A, dtype=float) # Create matrix f
            U = np.zeros_like(A, dtype=float) # Create matrix U

            #Find the first column of matrix f and Q
            f[:, 0] = A[:, 0]      
            Q[:, 0] = f[:, 0] / np.linalg.norm(f[:, 0]) 

            #Find the other columns of matrix f
            for i in range(1,n):
                for j in range(i):
                    f[:,i] +=  -(np.multiply(np.dot(A[:,i] , f[:, j]),f[:, j]))/(np.linalg.norm(f[:,j]))**2
                    
                f[:,i] = f[:,i]+A[:,i] #Add the negative summation of f to the column of A we are intrested   
    
                Q[:, i] = f[:,i] / np.linalg.norm(f[:, i]) # compute the other columns in Q
            
            #Find the matrix U
            for i in range(0,n):
                for j in range(0,n):
                    #The diagonal elements of U is always the norm of the matrix f
                    if i==j:
                        U[i,i]= np.linalg.norm(f[:,i])
                    #The upper diagonal elements of U are then found
                    if j>i:
                        U[i,j]=np.dot(A[:,j],Q[:,i])
        
            #We then compute UQ for each itteration        
            UQ = np.matmul(U,Q) 
            A=UQ #We redefine A after each itteration so we use the previous itterations QU matrix to find the next QU matrix
            
            Eigen =np.zeros(n) #Create an array for the eigenvalues
            
            #Write the diagonal elements of matrix UQ to the array Eigen
            for i in range(0,n):
                for j in range(0,n):
                    if i==j:
                        Eigen[i]=UQ[i,j]

            print(f'The eigenvalues are {Eigen[:]}, for itteration {counter}')      
            counter=1+counter       
        
        omega=np.sqrt(-Eigen[:])
        omega1= list(omega)
        print(f'The natural frequency of oscillation for each particle is {omega1}') #Print the natural frequenices of the last itteration

m1 = 8
m2 = 1
k = 8
matrix = np.array([[-2*k/m1,k/m1],[k/m2,-2*k/m2]])
Eigenvalue(2,25,False,matrix)
