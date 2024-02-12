#!/usr/bin/env python
# coding: utf-8

# In[1]:


#question 1
guess = 12

#question 2
a = 0
b = 1
n_s = 1000
x_gauss = [-0.861136,-0.339981,0.339981,0.861136]           #  for n =3, hard coded the values of weights and roots for n = 3, legendre polynomials
w_gauss = [0.347855,0.652145,0.652145,0.347855]

#question 3
h=[0.5,0.2,0.05,0.01]
x_ = 10

# question 4
Nx,Nt,Lx,Tx,alpha = 100,2000,8,1,0.03125 # alpha = 0.03125, LU decomposition is used to solve the equation/inverse the matrix 

# question 5
lx,ly,nx,ny,lx0,ly0,iterations = 2,1,6,6,0,0,1000   # lx,ly,lx0,ly0 give the boundaries and nx,ny give the grid sizes
                                                    # iterations is necessary to get accurate values by iteratively improving
                                                    # the solution to get the accurate enough result, else it will give very off answers
                                                    # in the only 1 iteration

