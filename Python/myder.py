# coding: utf-8

# In[74]:


import math as math
import numpy as np
import time
import matplotlib.pyplot as plt
from array import array
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import mpld3
mpld3.enable_notebook()

# define 1st order derivatives
def forwardiff(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(1,N):
        slop = (f(a+k*h)-f(a+(k-1)*h))/h
        g[k]=slop
    return g
    
def backwardiff(f,a,b,N):
    h = (b-a)/N
    g= np.zeros(N)
    for k in range(0,N-1):
        slop = (f(a+(k+1)*h)-f(a+k*h))/h
        g[k]=slop
    return g

def centraldiff(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(0,N):
        slop = (f(a+(k+1/2)*h)-f(a+(k-1/2)*h))/h
        g[k]=slop
    return g


def fordiff2(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(1,N):
        slop = (f(a+k*h)-f(a+(k-1)*h))/h
        g[k]=slop
    return g

# define 2nd order derivatives by central differentiation
def centraldiff2(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(0,N):
        slop = (f(a+(k+1)*h)-f(a+(k-1)*h))/(2*h)
        g[k]=slop
    return g

# define 3rd order derivatives by central differentiation
def centraldiff3(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(0,N):
        slop = (27*f(a+(k+1/2)*h)-27*f(a+(k-1/2)*h)+f(a+(k-3/2)*h)-f(a+(k+3/2)*h))/(24*h)
        g[k]=slop
    return g

# define 4th order derivative by central differentiation
def centraldiff4(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(0,N):
        slop = (8*f(a+(k+1)*h)-8*f(a+(k-1)*h)+f(a+(k-2)*h)-f(a+(k+2)*h))/(12*h)
        g[k]=slop
    return g

# define 5th order derivative by central differentiation
def centraldiff5(f,a,b,N):
    h = (b-a)/N
    g = np.zeros(N)
    for k in range(0,N):
        slop = ((75/64)*f(a+(k+1/2)*h)-(75/64)*f(a+(k-1/2)*h)-(3/640)*f(a+(k-5/2)*h)+(3/640)*f(a+(k+5/2)*h)+(25/384)*f(a+(k-3/2)*h)-(25/384)*f(a+(k+3/2)*h))/h
        g[k]=slop
    return g

# define higher order derivatives by forward differentiation method
def forhigh(f,a,b,N):
    h = (b-a)/N
    n = f.shape[0]
    g = np.zeros(n)
    for i in range(1,n-1):
        slop = (f[i]-f[i-1])/h
        g[i-1] = slop
    g[n-1] = 0    
    return g

# define higher order derivatives by backward differentiation method
def backhigh(f,a,b,N):
    h = (b-a)/N
    n = f.shape[0]
    g = np.zeros(n)
    g[0] = 0
    for i in range(2,n):
        slop = (f[i]-f[i-1])/h
        g[i-1] = slop
    return g
        


