#!/usr/bin/env python
# coding: utf-8

# # Cálculo de derivadas numéricas y errores
# ### Autor: Luis M. de la Cruz Salas
# ### Rev: 11/02/2019

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')


# * Aproximación de la primera derivada usando diferencias finitas hacía adelante (Forward):
# 
# $\displaystyle
# \frac{\partial u(x)}{\partial x} \approx \lim\limits_{h\to 0} \frac{u(x+h) - u(x)}{h}
# $

# In[3]:


def forwardFD(u,x,h):
    """ Esquema de diferencias finitas hacia adelante
    u : es la función a evaluar
    x : es la posición
    h : es el tamaño de malla.
    """
    return (u(x+h)-u(x))/h


# * Aproximación de la primera derivada usando diferencias finitas hacía atrás (Backward):
# 
# $\displaystyle
# \frac{\partial u(x)}{\partial x} \approx \lim\limits_{h\to 0} \frac{u(x) - u(x-h)}{h}
# $

# In[4]:


def backwardFD(u,x,h):
    """ Esquema de diferencias finitas hacia atrás
    u : es la función a evaluar
    x : es la posición
    h : es el tamaño de malla.
    """
    return (u(x)-u(x-h))/h


# * Aproximación de la primera derivada usando diferencias finitas hacía centradas (Centered):
# 
# $\displaystyle
# \frac{\partial u(x)}{\partial x} \approx \lim\limits_{h\to 0} \frac{u(x+h) - u(x-h)}{2h}
# $

# In[5]:


def centeredFD(u,x,h):
    """ Esquema de diferencias finitas centradas
    u : es la función a evaluar
    x : es la posición
    h : es el tamaño de malla.
    """
    return (u(x+h)-u(x-h))/(2*h)


# In[6]:


def u(x):
    """ Esta función es la que deseo evaluar
    """
    return np.sin(x)


# In[7]:


def uprime(x):
    """ Esta es la derivada de la función f(x) """
    return np.cos(x)


# In[8]:


def error(ue, un):
    return np.fabs(un - ue)


# In[11]:


h = np.array([0.1,0.05,0.01,0.005,0.001])
h


# In[13]:


np.ones(5)


# In[14]:


ef = error(uprime(np.ones(5)),forwardFD(u,1,h))
eb = error(uprime(np.ones(5)),backwardFD(u,1,h))
ec = error(uprime(np.ones(5)),centeredFD(u,1,h))
print(ef,eb,ec)


# In[15]:


import pandas as pd

dframe = pd.DataFrame(np.array([h,ef,eb,ec]))
dframe


# In[16]:


dframe = pd.DataFrame(np.array([h,ef,eb,ec]).transpose())
dframe


# In[17]:


dframe = pd.DataFrame(np.array([h,ef,eb,ec]).transpose(), columns=['$h$','$D_+$', '$D_-$','$D_0$'])
dframe
print(dframe)

# In[18]:


plt.plot(h,ef,'o-',label='Forward')
plt.plot(h,eb,'s-',label='Backward')
plt.plot(h,ec,'v-',label='Centered')
plt.loglog()
plt.legend()
plt.xlabel('log($h$)')
plt.ylabel('log(error)')
plt.show()


# In[20]:


N = 10
A = 0
B = 2 * np.pi
x = np.linspace(A,B,N)
h = (B-A) / (N-1)

plt.subplot(211)
plt.plot(x,u(x),'o-')
plt.grid()

plt.subplot(212)
plt.plot(x,uprime(x),'s-')
plt.plot(x,forwardFD(u,x,h), label='Forward')
plt.plot(x,centeredFD(u,x,h), label='Centered')
plt.grid()
plt.legend()

plt.subplots_adjust(hspace=0.5)
plt.show()


# In[ ]:




