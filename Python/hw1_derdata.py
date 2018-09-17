
# coding: utf-8

# In[17]:


import math as math
import numpy as np
import matplotlib.pyplot as plt
from array import array
import matplotlib.pyplot as plt
import scipy.io as spi
get_ipython().run_line_magic('matplotlib', 'inline')
import mpld3
mpld3.enable_notebook()

derdata = spi.loadmat('/Users/anindita/Desktop/graduate studies/NEU/academics/Fall 18/Computational physics/Matlab codes/derdata.mat')
x=derdata['X'];
y=derdata['Y'];
dy = np.diff(y,1,0)/np.diff(x,1,0)
n = dy.shape[0]
ddy = np.diff(np.diff(y,1,0),1,0)/np.diff(np.diff(x,1,0),1,0)
ddy_alter = np.diff(y,2,0)/np.diff(x,2,0)
n1 = ddy.shape[0]
n2 = ddy_alter.shape[0]

# plot Y, Y'(x) and Y''(X)
plt.plot(x,y,'g+-', label='Y')
plt.plot(x[0:n],dy, 'y^-', label='Y\'(X)')
plt.plot(x[0:n1],ddy, 'b*-', label='Y\"(X) by double application of diff function')
plt.plot(x[0:n2],ddy_alter, 'r--', label='Y\"(X) by original diff function')
plt.xlabel('X')
plt.legend()
plt.show()

# compare 2nd derivatives by two methods
plt.plot(x[0:n1],ddy, 'b*-', label='Y\"(X) by double application of diff function')
plt.plot(x[0:n2],ddy_alter, 'r--', label='Y\"(X) by original diff function')
plt.xlabel('X')
plt.title('Comparison of 2nd derivatives')
plt.legend()
plt.show()

