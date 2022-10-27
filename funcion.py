from cmath import pi
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import signal
import streamlit as st 
from scipy.integrate import quad
from math import *


#def sc(r,k):
 #   h=k*np.cos(i*r)
  #  return h

#def ss(r,k):
 #   h=k*np.sin(i*r)
  #  return h 

def senabs(li,ls,a,f,p):
    t=np.arange(li,ls,p)
    y=np.abs(a*np.sin(f*t))
    n="np.abs(a*np.sin(f*t))"
    return t,y,n

#x,y,fc,fs= senabs(-np.pi,np.pi,1,1,0.001)

t,y= senabs(-np.pi,np.pi,1,1,0.001)
k=np.abs(np.sin(t))
def fc(t):  
    h=k*np.cos(i*t)
    return h
def fs(t):
    h=k*np.sin(i*t)
    return h 

#fc=lambda x:square(x)*cos(i*x)  

#fs=lambda x:square(x)*sin(i*x)

#x=np.arange(-np.pi,np.pi,0.001)

#y=np.abs(np.sin(x))


n=5

An=[] 

Bn=[]

sum=0

for i in range(n):

    an=quad(fc,-np.pi,np.pi)[0]*(1.0/np.pi)

    An.append(an)

for i in range(n):

    bn=quad(fs,-np.pi,np.pi)[0]*(1.0/np.pi)

    Bn.append(bn) 

for i in range(n):

    if i==0.0:
        sum=sum+An[i]/2

    else:

        sum=sum+(An[i]*np.cos(i*t)+Bn[i]*np.sin(i*t))



plt.plot(t,sum,'g')

plt.plot(t,y,'r--')

plt.title(("serie de furier con",n,"armonicos"))

plt.show()