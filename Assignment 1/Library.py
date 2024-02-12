#!/usr/bin/env python
# coding: utf-8

# # Import functions

# In[1]:


import math
import matplotlib.pyplot as plt
import copy
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


# # Library started

# Miscellaneous functions

# In[2]:


# Sum of an AP: 
def AP(a,cd,n):                                  
    term=a  
    sum2=term
    for i in range(n-1):                                                         
        term=term+cd
        sum2=sum2+term
    return sum2


#Sum of a HP    
def HP(h,cd,n):
    sum1=1/h
    term=1/h
    for i in range(n-1):
        term=1/term+cd
        term=1/term
        sum1+=term
    return sum1


#Sum of GP  
def GP(g,r,n):
    i=1
    sum1=0
    while i<=n:
        sum1=sum1+g
        g=g*r
        i=i+1
    term=g/r
    return sum1


#Sum of 1st n odd numbers   
def sumodd(n):
    j=1
    sum1=0
    i=1
    while j<=n:
        sum1=sum1+i
        i=i+2
        j=j+1
    return sum1


#Factorial of N:    
def factorial(n):
    if n<0:
        return None
    i=1
    fact=n
    if n==0:
        print("The factorial is : 1")
    else:
        while i<n:
            fact=fact*i
            i=i+1
        return fact

    
##################################################### COMPLEX OPERATIONS ###################################################

#Class for operations of complex numbers 
class mycomlpex:

    def __init__(self,r1,i1,r2,i2):
        self.a=r1
        self.b=i1
        self.c=r2
        self.d=i2

    def print_number(self,):
        print("z1=",self.a,"+ i",self.b," and z2=",self.c,"+",self.d,"i")

    def mod(self):
        z1=(self.a**2+self.b**2)**(0.5) 
        z2=(self.c**2+self.d**2)**(0.5)
        print("The modulus are: z1:",z1,"and z2:",z2)

    def addition(self):
        r=(self.a+self.c)
        i=(self.b+self.d)
        if i>=0:
            print("Addition:",r,"+",i,"i")
        else:
            print("Addition:",r,i,"i")

    def substraction(self): 
        r=(self.a-self.c)
        i=(self.b-self.d)
        if i>=0:
            print("Substraction:",r,"+",i,"i")
        else:
            print("Substraction:",r,i,"i")

    def product(self):
        r=(self.a*self.c - self.b*self.d)
        i=(self.b*self.c + self.a*self.d)
        if i>=0:
            print("Product:",r,"+",i,"i")
        else:
            print("Product :",r,i,"i")


#  Matrix basic operations - Corrected and checked

# In[3]:


#Program to swap individual rows    [Correct]
def swap_row(A,i,j):
    temp=A[i]
    A[i]=A[j]
    A[j]=temp
    return A


#Program to find sum of 2 matrices.
def matrixsummation(A,B):
    C = []
    for i in range(len(A)):                     #In length of A-matrix, C-matrix dimension is chosen.
        C.append([])
        for j in range(len(A[0])):
            C[i].append(A[i][j] + B[i][j])      #For each term (row) of the C matrix, the sum of indivisual terms of A and B are added and appended to the jth row.
    return C


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


#Dot product - This has to be 2 column vectors multiplied together.
def dotproduct(C,D):
    i=0
    sum1=0
    while i<len(C):                 #Iteration with i<[length of the matrix]
        sum1+=C[i]*D[i] 
        i=i+1
    return sum1
 
#Prints matrix 
def printmatrix(A):
    for i in A:
        print(i)


# Numerical intergrations

# In[4]:


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


# In[5]:


def trapezoid_function(x):
    y = x**(2) 
    return y

# Trapezoid method with n divisions [correct]
def trapezoid(upper,lower,n):
    h = (upper - lower)/n
    area = h*(midpoint_function(lower) + midpoint_function(upper))/2 
    i = 1
    counter = 2
    while i < n:
        counter = counter + 1
        x = h*midpoint_function(lower + i*h)
        area = area + x
        i+=1    
    return area,counter


# In[6]:


def simpson_function(x):
    y = math.sqrt(1+x**(4))
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


# In[7]:


#SOL2

# let the n =4
# function
def gaussian_func(x):
    return (1 + x**4)**(0.5)
def gaussian_quadrature(a,b):
# Define the two-point Gaussian quadrature weights and abscissas
    x = [-0.861136,-0.339981,0.339981,0.861136]
    w = [0.347855,0.652145,0.652145,0.347855]
    integral = 0
# Perform the integral transformation
    for i in range(len(x)):
        integral = integral + w[i]*gaussian_func((b-a)/2 * (x[i]) + (a+b)/2)
        
    integral =integral*(b-a)/2 
# Print the result
    print(f"The integral through gaussian quadrature is: {integral:.6f}")


# Solutions to nonlinear equations

# In[8]:


def fixedpoint_function(x):
    return math.exp(-x)

def fixedpoint_method(guess):
    x = 0
    y = guess
    counter = 0
    while abs(y - x)>0.00001:
        counter = counter + 1
        x = y
        y = fixedpoint_function(x)
    return y,counter
        


#     ODE

# In[9]:


def ODEfunction(x,y):                                 # Defining the differential equation
    return (5*x**(2) - y)/(math.exp(x+y))

def RK4_uncoupled(x,y0,h,x1):                         # Runge-Kutta-4 for uncoupled equations and one boundary condition, x0,x1 are boundaries and y0 is the boundary condition
    X=[]
    Y=[]
    i = x
    y=y0
    while i <= x1:
        Y.append(y)
        k1 = h*ODEfunction(i,y)
        k2 = h*ODEfunction(i+h/2,y+k1/2)
        k3 = h*ODEfunction(i+h/2,y+k2/2)
        k4 = h*ODEfunction(i+h,y+k3)
        X.append(i)
        i = i+h
        y = y+(k1 + 2*k2 + 2*k3 + k4)/6
    return X,Y

def printfunction(X,Y):
    plt.plot(X, Y)
    plt.show()


# Poisson's equation

# In[10]:


def poisson_func(x,y,i,j):
    return x[i] * np.exp(y[j])
    

def poisson(lx,ly,nx,ny,lx0,ly0):

    # Creating the grid
    x = np.linspace(0, 2, nx)
    y = np.linspace(0, 1, ny)
    hx = (lx-0)/(nx-1)
    hy = (ly - 0)/(ny-1)
    # Initialize the solution matrix
    u = np.zeros((nx, ny))

    # Set boundary conditions
    u[:, 0] = x
    u[:, -1] = x * np.exp(1)
    u[0, :] = 0
    u[-1, :] = 2 * np.exp(y)

    # Solve the Poisson's equation using finite difference
    for it in range(1000):
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                u[i, j] = 0.5*((u[i+1, j] + u[i-1, j])*(hy**(2)) + (u[i, j+1] + u[i, j-1])*(hx**(2)) - (hx**(2) * hy**(2))* poisson_func(x,y,i,j))/(hx**(2)+ hy**(2))         #Considering hx = hy

    # Display the solution in a table
    print("Solution in a Table:")
    for row in u:
        for value in row:
            print(f"{value:.4f}", end="\t")
        print()

    # Display the solution in a 3-D plot
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, u, cmap='viridis')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Solution')
    ax.set_title("Solution to Poisson's Equation")
    plt.show()
    
    return x,y,u


# In[11]:


from scipy.linalg import solve

def cranck(L,T,nx,nt,alpha):
# Discretize the spatial and temporal domains
    x = np.linspace(0, L, nx+1)
    t = np.linspace(0, T, nt+1)
# Set up the grid
    hx = L / nx
    ht = T / nt

# Initialize the solution matrix
    u = np.zeros((nx+1, nt+1))
    # Set initial condition
    u[:, 0] = 4 * x - x**2 / 2
    u[0, :] = 0
    u[7, :] = 0

    # Constructging the tridiagonal matrix
    A = np.diag(1 + 2*alpha*np.ones(nx-1)) + np.diag(-alpha*np.ones(nx-2), k=-1) + np.diag(-alpha*np.ones(nx-2), k=1)

    print("A = ",A)

    # Time-stepping loop
    for n in range(0, nt):
        b = np.copy(u[1:-1, n])
        b[0] += alpha * (u[0, n] + u[0, n+1])
        b[-1] += alpha * (u[-1, n] + u[-1, n+1])                        #Using gauss sieldel from numpy library to invert the matrix
        u[1:-1, n+1] = solve(A, b)

    # Display the solution in a table
    print("Solution in a Table:")
    for i in range(nx+1):
        for j in range(nt+1):
            print(f"{u[i, j]:.4f}", end="\t")
        print()

# Display the solution in a contour plot
    X, T = np.meshgrid(x, t)
    fig, ax = plt.subplots()
    contour = ax.contourf(X, T, u, cmap='viridis')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Time')
    ax.set_title('Solution to Heat Equation using Crank-Nicolson')
    plt.colorbar(contour)
    plt.show()

