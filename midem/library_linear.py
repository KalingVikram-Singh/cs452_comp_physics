#!/usr/bin/env python
# coding: utf-8

# # Import functions

# In[1]:


import math
import matplotlib.pyplot as plt
import copy
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# # Library

#  Matrix basic operations - Corrected and checked

# In[163]:


#Program to find product of 2 matrices. 
def matrixmultiply(A, B):
    C = []                                      
    for i in range(len(A)):
        C.append([])                            #In length of A-matrix, C-matrix dimension is chosen.
        for j in range(len(A[0])):
            C[i].append(0)                      #Getting 0 in all indexes of C.
            for k in range(len(A[0])):
                C[i][j] = round(A[i][k] * B[k][j]+C[i][j],1)    #For each term (row) of the C matrix, the product of indivisual terms of A and B are added and appended to i,j index of C.
    return C


#Check if a matrix is symmetric 
def is_symmetric(A):                                            
    i=0
    counter=0
    while i<len(A):
        j=0
        while j<len(A):
            if A[i][j]==A[j][i]:
                counter+=1    
            j+=1
        i+=1
    if counter==len(A)*len(A):
        return True
    else:
        return False
    
#Return the transpose of a matrix
def transpose(A):       
    i=0
    B=copy.copy(A)
    while i<len(B):
        j=i+1
        while j<len(A):
            temp=B[i][j]
            B[i][j]=B[j][i]
            B[j][i]=temp
            j+=1
        i+=1
    return B


# ### Linear Algebra (All)

# In[167]:


# A=[[4,-1,0,-1,0,0],[ -1,4,-1,0,-1,0],[0,-1,4,0,0,-1],[-1,0,0,4,-1,0],[0,-1,0,-1,4,-1],[0,0,-1,0,-1,4]]
# B=[2,1,2,2,1,2]

# L=cholesky(A)                                           #Calls Cholesky
# print_matrix(A)                                         #Prints the decomposed lower matrix L

# print()
# Y2=forwardsub_chol(L,B)                                      #Calls Forward substitution

# print("Y, after forward substitution is :",Y2)
# U=transpose(A)                                          #Take upper matrix for backward substitution
# X2=backwardsub_chol(U,Y2)
# print()
# print("The solution by Cholesky is given by X=",X2)      #Solution through Cholesky
# print()


# In[3]:


#Function for LU decomposition
def print_matrix(A):
    for i in A:
        print(i)
        
def forwardsub_LU(A_fsub,B_fsub):
    n = len(A_fsub)
    x = [0] * n
    for i in range(n):
        x[i] = B_fsub[i]
        for j in range(i):
            x[i] -= A_fsub[i][j] * x[j]
        x[i] /= A_fsub[i][i]

    return x

def backwardsub_LU(matrix, b):
    n = len(matrix)
    x = [0] * n

    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] -= matrix[i][j] * x[j]
        x[i] /= matrix[i][i]

    return x

def ludecomposition(A_lu):
    lowertriangle = [[0 for x in range(len(A_lu))]
             for y in range(len(A_lu))]
    uppertriangle = [[0 for x in range(len(A_lu))]         #We define matrices beforehand
             for y in range(len(A_lu))]

    for i in range(len(A_lu)):                 #We decompose matrix into lower and upper
        for k in range(i, len(A_lu)):
            sum1 = 0
            for j in range(i):
                sum1 += (lowertriangle[i][j] * uppertriangle[j][k])             #Using recurrence to get individual sums
            uppertriangle[i][k] = A_lu[i][k] - sum1    
 
        for k in range(i, len(A_lu)):
            if (i == k):
                lowertriangle[i][i] = 1                                         #Making diagonal terms 1 for solving equation.
            else: 
                sum1 = 0                                                        #Using recurrence to get individual sums
                for j in range(i):
                    sum1 += (lowertriangle[k][j] * uppertriangle[j][i])
 
                lowertriangle[k][i] = ((A_lu[k][i] - sum1)/uppertriangle[i][i])
    return lowertriangle,uppertriangle


# In[4]:


# # #Input@1
# E=[[1,-1,4,0,2,9],[0,5,-2,7,8,4],[1,0,5,7,3,-2],[6,-1,2,3,0,8],[-4,2,0,5,-5,3],[0,7,-1,5,4,-2]]
# F=[19,2,13,-7,-9,2]

    

# L,U=ludecomposition(E) 
# print("Lower matrix:")
# print_matrix(L) 
# print()
# print("Upper matrix:")
# print_matrix(U)

# Y1=forwardsub_LU(L, F)              # We call for forward and then backward substitution.
# X1=backwardsub_LU(U, Y1)
# print()
# print("The solution from LU decomposition is :",X1)      # Here, we have the solution
# print()

