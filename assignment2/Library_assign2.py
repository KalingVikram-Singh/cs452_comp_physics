#!/usr/bin/env python
# coding: utf-8

# In[4]:


import math
import matplotlib.pyplot as plt
import copy
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# In[5]:


# Question 1 (Cholesky)

#This function decomposes A into Lower and Upper which are transpose of each other and checks if A is symmetric.
#Check if a matrix is symmetric 
def is_symmetric(A_sm):                                            
    i=0
    counter=0
    while i<len(A_sm):
        j=0
        while j<len(A_sm):
            if A_sm[i][j]==A_sm[j][i]:
                counter+=1    
            j+=1
        i+=1
    if counter==len(A_sm)*len(A_sm):
        return True
    else:
        return False
    
#Return the transpose of a matrix
def transpose(A_tr):       
    i=0
    B_tr=copy.copy(A_tr)
    while i<len(B_tr):
        j=i+1
        while j<len(A_tr):
            temp=B_tr[i][j]
            B_tr[i][j]=B_tr[j][i]
            B_tr[j][i]=temp
            j+=1
        i+=1
    return B_tr


def print_matrix(A_z):
    for i in A_z:
        print(i)
        
def forwardsub_chol(A_fsub,B_fsub):
    n = len(A_fsub)
    x = [0] * n
    for i in range(n):
        x[i] = B_fsub[i]
        for j in range(i):
            x[i] -= A_fsub[i][j] * x[j]
        x[i] /= A_fsub[i][i]

    return x

def backwardsub_chol(matrix, b):
    n = len(matrix)
    x = [0] * n

    for i in range(n - 1, -1, -1):
        x[i] = b[i]
        for j in range(i + 1, n):
            x[i] -= matrix[i][j] * x[j]
        x[i] /= matrix[i][i]

    return x

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
                A_cholmat[i][j]=0
                j+=1
            i+=1
        print()
        print("After Cholesky decomposition we have A:")     
        return A_cholmat
    else:
        return "A is not symmetric"
    
def truncate_float(f, n):
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d + '0' * n)[:n]])


# In[8]:


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

# print(f"The solution by Cholesky is given by X :",X2)      #Solution through Cholesky
# print(f"The solution by Cholesky is given by X (truncated upto 6th decimal place):")      #Solution through Cholesky
# for i in X2:
#     truncated_value = truncate_float(i, 6)
#     print(truncated_value)  


# In[10]:


# Question 1 (Gauss Siedel)

#function for gauss-seidel
def gauss_seidel(A_gs,B_gs,X):  
    n=len(A_gs)
    tol=10**(-6)
    X_old=[]
    diff=1
    for i in range(n):
        X_old.append(0)
    for k in range(100):                
        for i in range(n):                  
            sum1=B_gs[i]                   # Assigning b[i] before hand so as to reduce steps where we subtract b[i] from other number.    
            sum2=0
            j=0
            k=i+1
            while j<i:
                sum1-=A_gs[i][j]*X[j]
                j+=1
            while k<n:
                sum2-=A_gs[i][k]*X_old[k]                  #Using iterations
                k+=1
            X_old=copy.copy(X)
            X[i]=(sum1+sum2)/A_gs[i][i]
        diff=abs(X_old[i]-X[i])             #Check of tolerance
        if diff<tol:
            return X
    return "Non-ending"


# In[23]:


# C=[[4,-1,0,-1,0,0],[ -1,4,-1,0,-1,0],[0,-1,4,0,0,-1],[-1,0,0,4,-1,0],[0,-1,0,-1,4,-1],[0,0,-1,0,-1,4]]
# D=[2,1,2,2,1,2]
# X0=[0,0,0,0,0,-1]   #Guess

# sol = gauss_seidel(C,D,X0)
# print(f"The solution of the equations is given by:{sol}")
# print("The solution of the equations is given by(truncated upto 6th decimal place):")
# for i in sol:
#     truncate_value =truncate_float(i, 6)
#     print(truncate_value)


# In[14]:


#Function for LU decomposition
def print_matrix(A_pm):
    for i in A_pm:
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


# In[24]:


# # #Input@1
# E=[[11,3,0,1,2],[0,4,2,0,1],[3,2,7,1,0],[4,0,4,10,1],[2,5,1,3,13]]
# F=[51,20,15,15,92]

    

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
# print(f"The solution by LU decomposition is given by X (truncated upto 6th decimal place):")      #Solution through Cholesky
# for i in X1:
#     truncated_value = truncate_float(i, 6)
#     print(truncated_value)  


# In[16]:


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


# In[25]:


# # Input
# G=[[11,3,0,1,2],[0,4,2,0,1],[3,2,7,1,0],[4,0,4,10,1],[2,5,1,3,13]]
# H=[51,20,15,15,92]
# L = [[0 for i in range(len(G[i]))] for i in range(len(G))]
# A_inp = copy.deepcopy(G)
# B_inp = H.copy()
# # 
# # To get solution:
# sol_g =gauss_jordan(A_inp,B_inp)
# print("The solution through gauss jordan is (truncated to 6th decimal place): ")
# for i in sol_g:
#     trunc_value = truncate_float(i,6)
#     print(trunc_value)


# In[19]:


# Question 3
    
def truncate_float(f, n):
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d + '0' * n)[:n]])


def matrix_multiply(X,Y):
    # Initialize the result matrix with zeros
    result = [[0] * len(Y[0]) for _ in range(len(X))]

    # Iterate through rows of X
    for i in range(len(X)):
        # Iterate through columns of Y
        for j in range(len(Y[0])):
            # Iterate through rows of Y
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result

def dot_product(v1, v2):                       #Gives dot product of two vectors
    return sum(x*y for x, y in zip(v1, v2))

def matrix_vector_mul(A, v):                   #Multiplies matrix and vectors
    return [dot_product(row, v) for row in A]

def vector_subs(v1, v2):                       #Substract 2 vector
    return [x - y for x, y in zip(v1, v2)]
 
def vector_addi(v1, v2):                        #Add two vector
    return [x + y for x, y in zip(v1, v2)]

def scalar_mul(s, v):                           #Multiply vector with a scalar
    return [s * x for x in v]

def conjugate_gradient(A, b, x0, tol, max_iter):
    r = vector_subs(b, matrix_vector_mul(A, x0))
    d = copy.deepcopy(r)
    rold = dot_product(r, r)

    for i in range(max_iter):
        Ad = matrix_vector_mul(A, d)
        alpha = rold / dot_product(d, Ad)
        x0 = vector_addi(x0, scalar_mul(alpha, d))
        r = vector_subs(r, scalar_mul(alpha, Ad))
        rnew = dot_product(r, r)
        if (rnew) < tol**2:
            break
        d = vector_addi(r, scalar_mul(rnew/rold, d))
        rold = rnew
    return x0


def inverse_via_conjugate_gradient(A, tol=1e-4, max_iter=1000):
    n = len(A)
    identity_matrix = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    inverse_matrix = []

    for i in range(n):
        # Each column of the identity matrix is used as 'b'
        b = [row[i] for row in identity_matrix]
        x0 = [-1] * n  # Initial guess
        x = conjugate_gradient(A, b, x0, tol, max_iter)
        inverse_matrix.append(x)

    # Transpose the result to get the inverse matrix
    inverse_matrix = list(map(list, zip(*inverse_matrix)))
    return inverse_matrix


# In[26]:


# # Question 3
# #Need to make code to invert the matrix

# I = [[2, -3, 0, 0, 0, 0],
#      [-1, 4, -1, 0, -1, 0],
#      [0, -1, 4, 0, 0, -1],
#      [0, 0, 0, 2, -3, 0],
#      [0, -1, 0, -1, 4, -1],
#      [0, 0, -1, 0, -1, 4]]

# J = [-5/3, 2/3, 3, -4/3, -1/3, 5/3]
# x0 = [1,0,0,0,0,0]
# tol_cg = 10**(-4)
# iter_cg = 1000
# x = conjugate_gradient(I,J, x0,tol_cg,iter_cg)
# print("Solution: ", x)
# print()
# print("The solution(truncated to 4th decimal):")

# for i in x:
#     trunc_value = truncate_float(i,4)
#     print(trunc_value)

# print()
# x0 = [0,0,0,0,0,1]   
# X_CG_I =  [0,1,0,0,0,0]

# A_inv = inverse_via_conjugate_gradient(I)
# print(f"The inverse of A is:")
# for i in A_inv:
#     print(i)

# print()
# print("To check if its the inverse, multiply A with A_inv to get an identity function.")
# print()
# print(f"The matrix porduct gives:")
# Prod =matrix_multiply(I,A_inv)
# for i in Prod:
#     print(i)
# print()
# print("Thus, we get the inverse with a upto 4th decimal place correction. The matrix A *A_inv gives us identity with 10^(-4) errors" )


# In[21]:


# Question 4
    
def truncate_float(f, n):
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d + '0' * n)[:n]])


def matrix_multiply(X,Y):
    # Initialize the result matrix with zeros
    result = [[0] * len(Y[0]) for _ in range(len(X))]

    # Iterate through rows of X
    for i in range(len(X)):
        # Iterate through columns of Y
        for j in range(len(Y[0])):
            # Iterate through rows of Y
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result

def dot_product(v1, v2):                       #Gives dot product of two vectors
    return sum(x*y for x, y in zip(v1, v2))

def matrix_vector_mul(A, v):                   #Multiplies matrix and vectors
    return [dot_product(row, v) for row in A]

def vector_subs(v1, v2):                       #Substract 2 vector
    return [x - y for x, y in zip(v1, v2)]
 
def vector_addi(v1, v2):                        #Add two vector
    return [x + y for x, y in zip(v1, v2)]

def scalar_mul(s, v):                           #Multiply vector with a scalar
    return [s * x for x in v]

def conjugate_gradient_onfly(b, N,x0, tol, max_iter,m):
    r = [0.0 for i in range(N)]
    
    for i in range(N):
        Ax0 = 0
        for j in range(N):
            Ax0 += A(i,j,m)*x0[j]  
        r[i] = b[i] - Ax0    

        
    d = copy.deepcopy(r)
    rold = dot_product(r, r)
    residue = []
    iterations = []
    Ad = [0]*N

        
    
    for i in range(max_iter):
        for j in range(N):
            product_2 = 0
            for k in range(N):
                product_2 += A(j,k,m)*d[k]   
            Ad[j] = product_2
            
#         print(dot_product(d,Ad))
        alpha = rold / (dot_product(d, Ad))
        x0 = vector_addi(x0, scalar_mul(alpha, d))
        r = vector_subs(r, scalar_mul(alpha, Ad))
        rnew = dot_product(r, r)
        residue.append(rnew)
        iterations.append(i)
        if (rnew) < tol**2:
            break
        temp = d.copy()
        d = vector_addi(r, scalar_mul(rnew/rold, temp))
        rold = rnew
    return x0,residue,iterations

def A(x,y,m):
    return (delta((x+1)%50,y%50) + delta((x-1) %50,y%50) - 2*delta(x%50,y%50))/2 + (m**(2))*delta(x%50,y%50)

def delta(x,y):
    if x==y:
        return 1 
    else:
        return 0

def inverse_via_conjugate_gradient_2(N, m,tol=1e-6, max_iter=1000):
    n = N
    identity_matrix = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    inverse_matrix = []

    for i in range(n):
        # Each column of the identity matrix is used as 'b'
        b = [row[i] for row in identity_matrix]
        x0 = [5] * n  # Initial guess
        x,residue,iterations = conjugate_gradient_onfly( b, N , x0, tol, max_iter,m)
        inverse_matrix.append(x)

    # Transpose the result to get the inverse matrix
    inverse_matrix = list(map(list, zip(*inverse_matrix)))
    return inverse_matrix,residue,iterations


# In[27]:


# # Question 3
# #Need to make code to invert the matrix
# N = 50
# m = 0.2
# A_inv,residue,iterations = inverse_via_conjugate_gradient(N,m)
# print(f"The inverse of A is:")
# print()
# for i in A_inv:
#     print(i)
#     print()

# #  Plot the solution
# plt.title("Residue vs Iterations")
# plt.xlabel("Iterations")
# plt.ylabel("Residue")
# plt.plot(iterations, residue, color="red", marker="o", label="residue points")
# plt.legend()
# plt.show()

