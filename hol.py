from cmath import pi
from ctypes import sizeof
from re import M
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import streamlit as st 
import scipy as sp
import sympy as sym
from scipy import signal
from scipy import integrate 
from scipy.integrate import quad
import scipy.fftpack as fourier
from math import *
import scipy.io.wavfile as waves

def expo(T,a,p):
	t=np.arange(0,T,p)
	y=a*(np.exp(t))
	return t,y

def senabs(T,a,p):
    f=(2*np.pi)/T
    t=np.arange(0,T,p)
    y=np.abs(a*np.sin(f*t))
    return t,y

def senabsp(T,a,p):
    f=(2*np.pi)/T
    t=np.arange(-T,T,p)
    y=np.abs(a*np.sin(f*t))
    return t,y

def trian(T,a,p):
    f=(2*np.pi)/T
    t=np.arange(0,T,p)
    y = a*signal.sawtooth(f * t)
    return t,y

def trianp(T,a,p):
    f=(2*np.pi)/T
    t=np.arange(-T,T+0.01,p)
    y = a*signal.sawtooth(f * t)
    return t,y

def ramp(T,a,p):
    t = np.arange(0,T,p)
    y = t*0
    pend=a/(T-T*(2/3))
    y[t <=T-T*(2/3)]=(t[t<=T-T*(2/3)])*pend
    y[(t>T-T*(2/3)) & (t<T-(T/3))]=a
    y[t>=T-(T/3)]=(t[t>=T-T*(1/3)]-T)*(-pend)
    return t,y

def rec(T,a,p):
    f=(2*np.pi)/T
    t = np.arange(0,T,p)
    y = a*signal.square(f*t)
    return t,y

def recp(T,a,p):
    f=(2*np.pi)/T
    t = np.arange(-T-p,T,p)
    y = a*signal.square(f*t)
    return t,y


def serief(n,t,y,T,p):

    A = np.zeros((n))
    B = np.zeros((n))
    c = np.zeros(n)
    t5=np.zeros(n)
    fase= np.zeros(n)
    m=np.size(t)
    a0=0
    for i in range(m):
        a0=a0+(1/T)*y[i]*p

    for i in range(n):
        for j in range (m):
            A[i]=A[i] + ((2/T)*y[j]*np.cos(i*t[j]*f))*p   
            B[i]=B[i] + ((2/T)*y[j]*np.sin(i*t[j]*f))*p

        c[i] = ((A[i]**2)+(B[i]**2))**0.5
        t5[i]=i+1
        if A[i]==0:
            fase[i]= (np.pi)/2
        else: 
            fase[i] = -np.arctan(B[i]/A[i])

    t1=np.arange(-T,T,p)
    xf=0*t1-a0

    for i in range(n):
        xf= xf+A[i]*np.cos(i*f*t1)+B[i]*np.sin(i*f*t1)

    return t1,xf,c,fase,t5

def señalc(w0,w1,w2,fs,A0,A1,A2):
    Ts=1/fs  
    t=np.arange(0,1,Ts)
    y=A0*np.sin(w0*t) + A1*np.cos(w1*t) + A2*np.sin(w2*t)
    return t,y

def señald(w0,w1,w2,fs,A0,A1,A2,n):
    t=np.arange(0,n-1,1)
    y=A0*np.sin(w0*t/fs) + A1*np.cos(w1*t/fs) + A2*np.sin(w2*t/fs)
    return t,y

def tf(y,fs):
    ft=fourier.fft(y)
    oft=sp.fft.fftshift(ft) 
    hzft=abs(oft)
    phift=np.angle(oft)
    t1=fs*np.arange(-0.5-(1/len(hzft)),0.5-(1/len(hzft)),(1/len(hzft)))
    t2=fs*np.arange(-0.5-(1/len(phift)),0.5-(1/len(phift)),(1/len(phift)))
    return t1,t2,hzft,phift

fig, ax = plt.subplots()
def graph3(r,y):
	fig, ax = plt.subplots() #graph para cont (los invocas en el boton)
	ax.plot(r,y)
	ax.legend()
	the_plot3.pyplot(fig)
def graph(r,y,r2,y2):
        fig, ax = plt.subplots()
        ax.plot(r, y, label="Signal") #graph para serie
        ax.plot(r2, y2, label="System")
        ax.legend()
        the_plot.pyplot(fig)
def graph2(r,y):
	fig, ax = plt.subplots()
	ax.plot(r, y, label="Convolution") #graph para cont
	ax.legend()
	the_plot2.pyplot(fig)
def graph2disc(r,y):
	fig, ax = plt.subplots()
	ax.scatter(r, y, label="Convolution") #graph para disc
	ax.legend()
	the_plot2.pyplot(fig)
def graphdisc(r,y):
	fig, ax = plt.subplots()
	ax.scatter(r, y, label="Convolution") #graph para disc
	ax.legend()
	the_plot1.pyplot(fig)
def graph3disc(r,y):
	fig, ax = plt.subplots()
	ax.scatter(r, y, label="Convolution") #graph para disc
	ax.legend()
	the_plot3.pyplot(fig)
def graph1(r,y):
	fig, ax = plt.subplots()
	ax.plot(r, y, label="Convolution") #graph para cont
	ax.legend()
	the_plot1.pyplot(fig)
def option_type():
	option_type = st.selectbox(
	'Type',
	('En que lo puedo ayudar?','Serie', 'Transformada',))
	return option_type

ones='1.'
one='1'
two='2'
three= '3'
p1=0.04
option_type=option_type()
def graph(r,y,r2,y2):
	fig, ax = plt.subplots()
	ax.plot(r, y, label="Signal")
	ax.plot(r2, y2, label="System")
	ax.legend()
	the_plot.pyplot(fig)
def type(one):
	option1 = st.selectbox(
        'Funcion'+one,
        ('Exponencial', 'Sinusoidal rectificada', 'Triangular', 'Rectangular periódica', 'Rampa trapezoidal')) #lo que cambis aqui toca cambiarlo en los if de abajo
	return option1
def discorcont(one):
	option2 = st.selectbox(
        'Funcion'+one,
        ('Discreta','Continua')) #lo que cambies aqui toca cambiarlo en los if de abajo
	return option2

def contrec(which):

    st.header("Funcion"+': '+which)
    li = st.slider(' Periodo'+' '+which,1,10, 5)
    a = st.slider('Amplitud'+' '+which,1,20, 10)
    n = st.slider('Número de armónicos'+' '+which,1,100, 10)
    return li,a,n

def contrec2(which):

    st.header("Funcion"+': '+which)
    f = st.slider(' frecuencia'+' '+which,1,1000,500)
    a = st.slider('Amplitud'+' '+which,1,20, 10)
    return f,a

col1, col2 = st.columns(2)
with col1:	
    p = 0.001
    if option_type== 'Serie':
        a =type(one)
        if a=='Sinusoidal rectificada':  #a estos if me referia
            T,a,n = contrec(one)
            
            t,y = senabs(T,a,p)
            r,k = senabsp(T,a,p)
            t1,ysf,c,fase,t5= serief(n,t,y,T,p)

        if a== 'Triangular':
            li,a,n =contrec(one)
            t,y = trian(T,a,p)
            r,k = trianp(T,a,p)
            t1,ysf,c,fase,t5= serief(n,t,y,T,p)
        if a=='Rectangular periódica':
            li,a,n = contrec(one)
            t,y = rec(T,a,p)
            r,k = recp(T,a,p)
            t1,ysf,c,fase,t5= serief(n,t,y,T,p)

        if a=='Exponencial':
            li,a,n =contrec(one)
            t,y = expo(T,a,p)
            r,k = expo(T,a,p)
            t1,ysf,c,fase,t5= serief(n,t,y,T,p)

        if a=='Rampa trapezoidal':
            li,a,n=contrec(one)
            t,y = rec(T,a,p)
            r,k = recp(T,a,p)
            t1,ysf,c,fase,t5= serief(n,t,y,T,p)
        
    if option_type== 'Transformada':
        c = discorcont(one)
        if c=='Discreta':
            w0,A0 = contrec2(one)
            w1,A1 = contrec2(two)
            w2,A2 = contrec2(three)
            st.header('------------------------------')
            nm1 = st.slider('# de muestras', 0, 50, 20)
            fs = st.slider(' Frecuencia de muestreo', 1, 1000, 500)
            t2,y2=señald(w0,w1,w2,fs,A0,A1,A2,nm1)
            t3,t4,Hz,Ph=tf(y2,fs)
            ny = fs/2
                
        else:
            w0,A0 = contrec2(one)
            w1,A1 = contrec2(two)
            w2,A2 = contrec2(three)
            st.header('------------------------------')
            fs = st.slider(' Frecuencia de muestreo', 1, 1000, 500)
            t2,y2=señalc(w0,w1,w2,fs,A0,A1,A2)
            t3,t4,Hz,Ph=tf(y2,fs)
            ny = fs/2
with col2:
    if st.button('Graficar '): 
        if option_type== 'Serie':  
            try:
                graph(r,k,t1,ysf)
                graph3(t5,c)
                graph1(t5,fase)
            except:
                pass
        if option_type== 'Transformada':
                try:
                    if (w0 & w1 & w2)<= ny:
                        st.header("VOLANDOSE NYQUIST MI REY?")
                        
                    else:

                        if c == 'Discreta':
                            
                            graph3disc(t2,y2)
                            graph2disc(t3,Hz)
                            graphdisc(t4,Ph)

                        else:
                            graph2(t2,y2)
                            graph1(t3,Hz)
                            graph3(t4,Ph)

                except:

                    pass