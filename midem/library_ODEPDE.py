#!/usr/bin/env python
# coding: utf-8

# # Import functions

# In[2]:


import math
import matplotlib.pyplot as plt
import copy
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve


# # Library

# Matrix basic operations - Corrected and checked

# In[3]:


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


# # ODE

# BVP

# In[4]:


# implicit euler


# Gauss jordan for solving linear equations and inverse
# Partial pivoting
def pivot_gauss(C_pivs,D_pivs):
    C_piv,D_piv = copy.deepcopy(C_pivs),copy.deepcopy(D_pivs)
    for i in range(len(C_piv)-1):
        if C_piv[i][i] == 0:
            for j in range(i+1,len(C_piv)):
                if C_piv[j][i] > C_piv[i][i]:
                    C_piv[i],C[j] = C_piv[j],C[i]
                    D_piv[i],D_piv[j] = D_piv[j],D_piv[i]
    return C_piv,D_piv

def gauss_jordan(C_gjs,D_gjs):
    C_gj,D_gj = copy.deepcopy(C_gjs),copy.deepcopy(D_gjs)
    C_gj,D_gj=pivot_gauss(C_gj,D_gj)             #Doing partial pivoting
    for i in range(len(C_gj)):
        piv = C_gj[i][i]         #Making diagonal terms unity
        for j in range(i,len(C_gj)):
            C_gj[i][j] = C_gj[i][j]/piv
        D_gj[i] = D_gj[i]/piv   
        for k in range(len(C_gj)):              #Column to zero except diagonal
            if (k == i) or (C_gj[k][i] == 0):
                continue
            else:
                factor = C_gj[k][i]
                for l in range(i,len(C_gj)):
                    C_gj[k][l] = C_gj[k][l] - factor*C_gj[i][l]
                D_gj[k] = D_gj[k] - factor*D_gj[i]
    return copy.deepcopy(D_gj)

# Parameters
alpha = 0.01  # Thermal diffusivity
L = 1.0       # Length of the rod
T = 1.0       # Total time
Nx = 50       # Number of spatial points
Nt = 100      # Number of time steps

# Initial condition
def initial_condition(x):
    return np.sin(np.pi*x)

# Boundary conditions
def boundary_conditions(t):
    return 0.0  # Fixed temperature at both ends
def semi_implicit(alpha,L,T,Nx,Nt):
    dx = L / (Nx - 1)
    dt = T / Nt

    # Set up grid
    x = np.linspace(0, L, Nx)
    u = initial_condition(x)

    # Implicit method (Backward Euler)
    A = np.eye(Nx) - alpha * dt / dx**2 * np.diag(np.ones(Nx - 1), -1) + alpha * dt / dx**2 * np.diag(np.ones(Nx - 1), 1)
    A[0, 0] = A[-1, -1] = 1.0  # Boundary conditions

    for n in range(1, Nt + 1):
        b = u.copy()
        b[0] = b[-1] = boundary_conditions(n * dt)
        u = gauss_jordan(A, b)

    # Plot the results
    plt.plot(x, u, label='Final solution')
    plt.xlabel('x')
    plt.ylabel('Temperature')
    plt.legend()
    plt.show()


# In[5]:


# # Parameters
# alpha = 0.01  # Thermal diffusivity
# L = 1.0       # Length of the rod
# T = 1.0       # Total time
# Nx = 50       # Number of spatial points
# Nt = 100      # Number of time steps
# semi_implicit(alpha,L,T,Nx,Nt)


# In[6]:


#Solving ODE BVP using shooting with RK4

def dydx_sh_RK(x,y,z):
    return z
def dzdx_sh_RK(x,y,z):
    return -0.01*(20 - y)

def RK4_sol2(ul,y,z,N,par):  # y,z is initial value guess, ul is the end value of L, N is no of steps and 
    yi=y                                #par is for co0nvenience.
    zi=z 
    xi=0
    X=[xi]
    Y=[yi]
    Z=[zi]
    h=(ul - 0)/N
    while xi<=ul:

        k1y = h * dydx_sh_RK(xi , yi , zi)
        k1z = h * dzdx_sh_RK(xi , yi , zi)

        k2y = h* dydx_sh_RK(xi + h/2 , yi + k1y/2 , zi + k1z/2)
        k2z = h* dzdx_sh_RK(xi + h/2 , yi + k1y/2 , zi + k1z/2)

        k3y = h* dydx_sh_RK(xi + h/2 , yi + k2y/2 , zi + k2z/2)
        k3z = h* dzdx_sh_RK(xi + h/2 , yi + k2y/2 , zi + k2z/2)

        k4y = h* dydx_sh_RK(xi + h , yi + k3y , zi + k3z)
        k4z = h* dzdx_sh_RK(xi + h , yi + k3y , zi + k3z)

        xi+= h
        yi+= (k1y+2*k2y+2*k3y+k4y)/6
        zi+= (k1z+2*k2z+2*k3z+k4z)/6

        X.append(xi)
        Y.append(yi)
        Z.append(zi)

    if par==0:
        return X,Y,Z
    elif par==1:
        return zi,yi
    
def shooter(L,Y0,guess1,guess2,N,tempend):
    z1,y1=RK4_sol2(L,Y0,guess1,N,1)
    z2,y2=RK4_sol2(L,Y0,guess2,N,1)
    count=0
    while count < 50 and abs(y1-tempend)>=0.001 and abs(y2-tempend)>=0.001:
        count+= 1
        z_new = z2 + ((z1 - z2)*(tempend - y2))/(y1 - y2)         #interpolating the values to get the spot
        z_new2, y_new = RK4_sol2(L,Y0,z_new,N,1)
#         print(y_new,count)
        if abs(y_new - tempend) < 0.001:
            X,Y,Z= RK4_sol2(L,Y0,z_new,N,0)
            break
        else:
            if y_new < tempend:
                z2 = z_new
                y2 = y_new
            else:
                z1 = z_new
                y1 = y_new

    plt.plot(X,Y ,color = 'r')
    temp=0
    for i in range (0, len(Y)):
#         print(Y[i],X[i])
        if abs(Y[i] - 100) < 0.1:
            temp = i
            break
    plt.axvline(x = X[temp])
    plt.axhline(y = 100)
    plt.xlabel('L')
    plt.ylabel('T(C)')
    plt.show()
    print(f"The temperature is 100 celcius at L = {X[temp]:.6f}")   


# In[7]:


# up_lim = 10
# value_up_lim  =  200
# value_low_lim  = 40
# guess1 = -121
# guess2 = -122
# N = 10000
# shooter(up_lim,value_low_lim,guess1,guess2,N,value_up_lim )


# PDE

# Heat Solution with crack nicolson (with*beta uxx = ut) and given boundary conditions

# In[12]:


# Gauss jordan for solving linear equations and inverse
# Partial pivoting
def pivot_gauss(C_pivs,D_pivs):
    C_piv,D_piv = copy.deepcopy(C_pivs),copy.deepcopy(D_pivs)
    for i in range(len(C_piv)-1):
        if C_piv[i][i] == 0:
            for j in range(i+1,len(C_piv)):
                if C_piv[j][i] > C_piv[i][i]:
                    C_piv[i],C[j] = C_piv[j],C[i]
                    D_piv[i],D_piv[j] = D_piv[j],D_piv[i]
    return C_piv,D_piv

def gauss_jordan(C_gjs,D_gjs):
    C_gj,D_gj = copy.deepcopy(C_gjs),copy.deepcopy(D_gjs)
    C_gj,D_gj=pivot_gauss(C_gj,D_gj)             #Doing partial pivoting
    for i in range(len(C_gj)):
        piv = C_gj[i][i]         #Making diagonal terms unity
        for j in range(i,len(C_gj)):
            C_gj[i][j] = C_gj[i][j]/piv
        D_gj[i] = D_gj[i]/piv   
        for k in range(len(C_gj)):              #Column to zero except diagonal
            if (k == i) or (C_gj[k][i] == 0):
                continue
            else:
                factor = C_gj[k][i]
                for l in range(i,len(C_gj)):
                    C_gj[k][l] = C_gj[k][l] - factor*C_gj[i][l]
                D_gj[k] = D_gj[k] - factor*D_gj[i]
    return copy.deepcopy(D_gj)

def cranck_boundary(x,t,v,N):
    v[:, 0] = 0
    v[int(N/2),0] =300
    v[0, :] = 0
    v[-1, :] = 0
    return copy.copy(v)

def crank_nicolson(N,Nt,lx,tx,alpha,beta):
    h = lx / N
    k = tx / Nt 

# Initialize grid
    x = np.linspace(0, lx, N+1)
    t = np.linspace(0, tx, Nt+1)
    u = np.zeros((N+1, Nt+1))
#     print(u[1])
# Set initial conditions
    u = cranck_boundary(x,t,u,N)
    

# Crank-Nicolson method
    for j in range(Nt):
        A = np.diag(1 + beta*alpha*np.ones(N-1)) + np.diag(-beta*alpha*np.ones(N-2)/2, k=-1) + np.diag(-beta*alpha*np.ones(N-2)/2, k=1)
        B = np.diag(1 - beta*alpha*np.ones(N-1)) + np.diag(beta*alpha*np.ones(N-2)/2, k=-1) + np.diag(beta*alpha*np.ones(N-2)/2, k=1)
        u[1:N, j+1] =gauss_jordan(A, np.dot(B, u[1:N, j]))      # This gives the solution of the equation using 
                                                                    #  inverse of A thruogh gauss_jordan decompostion of:  
#     print(u)                                                              # Au(j+1) = Bu(j)
#     print(A)
# Display the solution in a table
    print("Solution in a Table:")
    for i in range(N+1):
        for j in range(Nt+1):
            print(f"{u[i, j]:.4f}", end="\t")


# Plot the solution
    X, T = np.meshgrid(x, t)
    plt.contourf(X, T, u.T, levels=20, cmap='coolwarm')
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('Temperature Distribution (Crank-Nicolson)')
    plt.show()


# In[14]:


# # Sol 3
# Nx,Nt,Lx,Tx,alpha,beta = 10,1000,2,2,0.05,1
# crank_nicolson(Nx,Nt,Lx,Tx,alpha,beta)


# In[ ]:




