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

# In[2]:


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


# ### Numerical intergrations (revise gaussian quadrature), else complete

# In[3]:


def midpoint_function(x):
    y = x**(2) 
    return y

# Midpoint method with n divisions [correct]
def midpoint(upper,lower,n):
    h = (upper - lower)/n
    area = 0
    i = 1
    counter = 0
    while i <= n:
        counter = counter + 1
        x = h*midpoint_function(lower + (2*i-1)*h/2)
        area = area + x
        i+=1    
    return area,counter


# In[169]:


# midpoint(3,0,100)


# In[8]:


def trapezoid_function1(x):
    y = x**(2) 
    return y

# Trapezoid method with n divisions [correct]
def trapezoid1(upper,lower,n):
    h = (upper - lower)/n
    area = h*(trapezoid_function1(lower) + trapezoid_function1(upper))/2 
    i = 1
    counter = 2
    while i < n:
        counter = counter + 1
        x = h*trapezoid_function1(lower + i*h)
        area = area + x
        i+=1    
    return area

def trapezoid_function2(x):
    y = x**(3) 
    return y

# Trapezoid method with n divisions [correct]
def trapezoid2(upper,lower,n):
    h = (upper - lower)/n
    area = h*(trapezoid_function2(lower) + trapezoid_function2(upper))/2 
    i = 1
    counter = 2
    while i < n:
        counter = counter + 1
        x = h*trapezoid_function2(lower + i*h)
        area = area + x
        i+=1    
    return area


# In[13]:


# # The mass is given by int_0_2(\lamdba*dx)
# upper = 2
# lower = 0
# n = 1000
# mass = trapezoid1(upper,lower,n)
# print(f"Mass of rod is:{mass:.4f}")

# # The COM is given by int_0_2(x.\lambda).dx/mass 
# COM = trapezoid2(upper,lower,n)/mass

# print(f"The COM is at :{COM:.4f}")


# In[171]:


# trapezoid(3,0,100)


# In[172]:


#sol 2
def simpson_function(x):
    y = x**(2)
    return y

# Simpson method with n divisions 
def simpson(upper,lower,n):
    area = 0
    h = (upper - lower)/n
    xi = [lower + i * h for i in range(n + 1)]
    area = h*(simpson_function(upper) + simpson_function(lower))  
    for  i in range(1,n,2):
        area += 4*h* simpson_function(xi[i])
    
    for i in range(2, n - 1, 2):
        area += 2*h* simpson_function(xi[i])
    area = area/3
    return area


# In[174]:


# simpson(3,0,100)


# In[1]:


#SOL2

# let the n = 3
# function
def gaussian_func(x):
    return (1 + x**4)**(0.5)
def gaussian_quadrature(a,b,x_gauss,w_gauss):
#         for n =3
#     x = [-0.861136,-0.339981,0.339981,0.861136] for n =3
#     w = [0.347855,0.652145,0.652145,0.347855]
    x = copy.copy(x_gauss)
    w = copy.copy(w_gauss)
    
    integral = 0
# Integrating the function by summing at the roots of the legendre function (0 - 1 limit )
    for i in range(len(x)):
        integral = integral + w[i]*gaussian_func((b-a)/2 * (x[i]) + (a+b)/2)
        
    integral =integral*(b-a)/2

# Print the result
    
    print(f"The integral through gaussian quadrature is: {integral:.6f}")

#     return integral


# In[2]:


# #question 2
# a = 0
# b = 1
# n_s = 1000
# x_gauss = [-0.861136,-0.339981,0.339981,0.861136]           #  for n =3, hard coded the values of weights and roots for n = 3, legendre polynomials
# w_gauss = [0.347855,0.652145,0.652145,0.347855]
# gaussian_quadrature(a,b,x_gauss,w_gauss)


# In[ ]:




