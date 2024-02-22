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

# In[5]:


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
def issymmetric(A):                                            
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

# In[6]:


#This function decomposes A into Lower and Upper which are transpose of each other if A is symmetric
def forwardsub_chol(A_fsub,B_fsub):
    i=0
    Y=[]
    for k in range(len(A_fsub)):
        Y.append(0)
    while i<len(A):
        j=0
        temp=0
        while j<i:
            temp+=A_fsub[i][j]*Y[j]
            j+=1
        Y[i]=(B_fsub[i]-temp)/A_fsub[i][i]
        i+=1
    return Y

#Program for backward substitution
def backwardsub_chol(A_bsub,Y_bsub):
    i=len(A_bsub)-1
    X=[]
    for l in range(len(A_bsub)):
        X.append(0)
    while i>=0:
        j=i+1
        temp=0
        while j<len(A_bsub):
            temp+=A_bsub[i][j]*X_bsub[j]
            j+=1
        X[i]=(Y_bsub[i]-temp)/A_bsub[i][i]
        i-=1
    return X

def cholesky(A_cholmat):
    
    if is_symmetric(A_cholmat):                                         #Check for symmetric matrix
        i=0
        while i <len(A_cholmat):
            j=0
            temp=0
            while j<i:
                temp+=A_cholmat[j][i]*A_cholmat[j][i]
                j+=1
            A_cholmat[i][i]=(A_cholmat[i][i]-temp)**(0.5)                       #Using recurrence relations
            j=i+1
            while j<len(A_cholmat):
                k=0
                temp=0
                while k<i:
                    temp+= A_cholmat[i][k]*A_cholmat[k][j]
                    k+=1
                A_cholmat[j][i]=(A_cholmat[j][i]-temp)/A_cholmat[i][i]                  #Using recurrence relations
                A_cholmat[i][j]=A_cholmat[j][i]
                j+=1
            i+=1
        i=0
        while i <len(A_cholmat):                                        #Making other indexes 0, as we want lower matrix
            j=i+1
            while j<len(A_cholmat):
                A[i][j]=0
                j+=1
            i+=1
        print()
        print("After Cholesky decomposition we have L:")       # We still have used a single variable A for the cholesky. Although for reference, we wrote L.
        return A_cholmat


# ### Solutions to nonlinear equations (complete)

# In[2]:


# regular falsi methhod to sovle non linear equations
def bracketing_regularfalsi(a,b):        
    value1=regularfalsi_function(a)
    value2=regularfalsi_function(b)
    max_iter=150
    c=(b-a)/10
    for i in range(max_iter):
        if value1*value2>0:
            if abs(value1)>abs(value2):
                b=b+c
            else:
                a=a-c
            value1=regularfalsi_function(a)
            value2=regularfalsi_function(b)
            if i==max_iter-1:
                print("chose different a,b")
                return None
    return(a,b)



def regularfalsi_function(x):
    return math.log(x/2) - math.sin(5*x/2)
def regularfalsi(a,b,tol):    
    max_iter = 100                           # let the maximum number of iterations be this
    delta=0
    X=[]
    Y=[]
    c=a
    if regularfalsi_function(b)*regularfalsi_function(b)>=0:
        if bracketing_regularfalsi(a, b)==None:
            return None
        a,b=bracketing_regularfalsi(a, b)
        print("After bracketing:",a,b)

    for i in range(max_iter):
        X.append(c)
        Y.append(i+1)
        # print("Iterations:",i+1,"solution:",c)
        temp=c
        c= b - ((b-a)*(regularfalsi_function(b))/(regularfalsi_function(b)-regularfalsi_function(a)))
        if abs((temp-c))<tol or abs(regularfalsi_function(a))<delta:
            # print("Iterations=", i)
            X.append(c)
            Y.append(i)
            break
        elif regularfalsi_function(a)*regularfalsi_function(c)<0:
            b=c
        elif regularfalsi_function(b)*regularfalsi_function(c)<0:
            a=c
    return c,i,X,Y


# In[3]:


# Newton raphson to solve non linear equations

def newtonraphson_function(x):
    return math.log(x/2) - math.sin(5*x/2)

def newtonraphson_functiondedrivative(x):
    return 1/x - (2.5)*math.cos(5*x/2)

def newtonraphson(x,tol):
    X=[]
    Y=[]
    if newtonraphson_function(x)==0:
        return x
    max_iter = 200                             #Maximum iterations
    delta=0
#     tol=10**(-6)                            # Tolerance level=10**(-4)
    for i in range(max_iter):
        value = newtonraphson_function(x)
        derivative = newtonraphson_functiondedrivative(x)
        x_new = x - value/derivative
        if abs(x_new - x)<tol or abs(newtonraphson_function(x))<delta:
            # print("Iterations =",i+2)
            X.append(x)
            Y.append(i+1)
            break
        x = x_new
        if i==max_iter - 1:
            return "try Another guess"
        X.append(x)
        Y.append(i+1)
    return x_new,i+1,X,Y


# In[7]:


# # Input
# # print(regularfalsi(a,b,tol))
# tol = 10**(-6)
# a = 1.1
# b = 2.5
# initial_guess = 5
# tol_raphson = 10**(-6)


# solu_regfalsi,iter_regfalsi,X_regfalsi,Y_regfalsi = regularfalsi(a,b,tol)
# print(f"The root through regular falsi is: {solu_regfalsi:.6f}")

# solu_newton,iter_newton,X_newton,Y_newton = newtonraphson(initial_guess,tol_raphson)
# print(f"The root through regular falsi is: {solu_newton:.6f}")

# plt.plot(Y_newton,X_newton ,color = 'r',label="Netwon raphson")
# plt.plot(Y_regfalsi,X_regfalsi ,color = 'b',label=" Regular falsi")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.show()


# In[ ]:




