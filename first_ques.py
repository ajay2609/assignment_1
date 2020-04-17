
# coding: utf-8

# In[188]:


import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import sys 
get_ipython().magic('matplotlib inline')
from matplotlib.pylab import figure


# In[198]:


import scipy.integrate as integrate 
from scipy.integrate import quad


def viogt(g,w,N,w_r):
        a=(g*w_r*10**(-3))/(b*8*np.pi)
        u=((1/w-1/w_r)*c)/(b/w_r)
        I=integrate.quad(lambda y:np.exp(-y**2)/(a**2 + (y-u)**2),-np.inf,np.inf)[0]
#        y=np.linspace(-1.0,1.0,1000000)
#        func=np.array(np.exp(-y**2)/(a**2 + (y-u)**2))  
#        I=integrate.trapz(func,y)
        return a/(b*np.sqrt(np.pi)) * 2.0*I * w_r *N * f12*10**(-15)  

c=3*10**5               #in km/s
b=12                    # in km/s
f12 = 0.42


# In[199]:


viogt = np.vectorize(viogt)
w_range =np.linspace(1214.0, 1217.0,1000)
for i in range(13,20):
    plt.plot(w_range,np.exp(-viogt(6.265*10**(-2),w_range,10**i,1215.67)),label='LogN=%i'%i)
plt.legend()    
plt.show()


# In[200]:


import re
for i in range(10,50,10):
    b=i
    plt.plot(w_range,[np.exp(-viogt(6.265*10**(-2),wav,10**13,1215.67)) for wav in w_range],label='b=%i'%i)

plt.title(r'voigt profiles for Ly-alpha Transition with N=e+13/$cm^2$')
plt.ylabel('normalized flux')
plt.xlabel('wavelength')    
plt.legend()
plt.savefig('viogt_b.png')
plt.show()    


# In[230]:


K=[]
K9=[]
logN=[]
p=9
while p<23:
    b=12
    f12=0.42
    w_range =np.linspace(1200.0, 1230.0,1000)
    K.append( [(1.0-np.exp(-viogt(6.265*10**(-2),l,10**p,1215.67))) for l in w_range] )
    b=30
    K9.append( [(1.0-np.exp(-viogt(6.265*10**(-2),l,10**p,1215.67))) for l in w_range] )
    logN.append(p)
    p+=0.5


# In[231]:


intgr=[]
intgr9=[]
for i in range(0,len(logN)):
    H=np.array(K[i])
    J=np.array(K9[i])
    intgr.append(np.log10(integrate.trapz(H,dx=0.002)))
    intgr9.append(np.log10(integrate.trapz(J,dx=0.002)))


# In[232]:


plt.figure(figsize=(5,5))
plt.plot(np.array(logN[0:-3]+np.log10(1215.67*0.42)) ,intgr[0:-3]-np.log10(1215.67),label=r"COG for b=12")
plt.plot(np.array(logN[0:-3]+np.log10(1215.67*0.42)) ,intgr9[0:-3]-np.log10(1215.67),label=r"COG for b=30")
plt.ylabel(r'log[W$\lambda$]')
plt.xlabel(r'Log[N $\lambda$ f]')
plt.legend()
plt.show()


# In[204]:


plt.plot(np.array(logN[0:-3]+np.log10(1215.67*0.42)) ,intgr1[0:-3]-np.log10(1215.67),label=r"COG for b=12")
plt.ylabel(r'log[W$\lambda$]')
plt.xlabel(r'Log[N $\lambda$ f]')
plt.title(r"Curve of Growth for Lymen-$\alpha$ line with b=12km/s and analytical COG ")
m=1
c=1-np.log10(1.132)-21.82
x=np.linspace(11,17,100)
y=m*x+c
plt.plot(x,y,'m:',label="linear region",)
m=0.5
c=-18-np.log10(1.88)+3.97378
x=np.linspace(21,23,100)
y=m*x+c
plt.plot(x,y,'-.',label="damped region")
x=np.linspace(17,22,100)
y=0.9*np.log10(x-(16))-4.3
plt.plot(x,y,'--',label="saturated region")
plt.legend()
plt.show()


# In[225]:


K1=[]
K2=[]
K3=[]
K4=[]
K=[]
logN=[]
b=15
p=13
while p<20:
    w_range =np.linspace(1606.5, 1610.5,1000)
    g=2.74*10**(-2)
    f12=0.0577
    K1.append( [(1.0-np.exp(-viogt(g,l,10**p,1608.45))) for l in w_range] )
    g=1.13*10**(-1)
    f12=0.133
    w_range =np.linspace(1524.7, 1528.7,1000)
    K2.append( [(1.0-np.exp(-viogt(g,l,10**p,1526.7))) for l in w_range] )
    g=2.88*10**(-2)
    f12=0.1278
    w_range =np.linspace(1332.5, 1336.5,1000)
    K3.append( [(1.0-np.exp(-viogt(g,l,10**p,1334.53))) for l in w_range] )
    g=5.65*10**(-2)
    f12=0.048
    w_range =np.linspace(1300.2, 1304.2,1000)
    K4.append( [(1.0-np.exp(-viogt(g,l,10**p,1302.16))) for l in w_range] )
    g=2.265*10**(-2)
    f12=0.42
    w_range =np.linspace(1213.7, 1217.7,1000)
    K.append( [(1.0-np.exp(-viogt(g,l,10**p,1215.67))) for l in w_range] )
    logN.append(p)
    p+=0.1


# In[226]:


intgr1=[]
intgr2=[]
intgr3=[]
intgr4=[]
intgr5=[]
for i in range(0,len(logN)):
    H=np.array(K1[i])
    J=np.array(K2[i])
    M=np.array(K3[i])
    q=np.array(K4[i])
    z=np.array(K[i])
    intgr1.append(np.log10(integrate.trapz(H,dx=0.002)))
    intgr2.append(np.log10(integrate.trapz(J,dx=0.002)))
    intgr3.append(np.log10(integrate.trapz(M,dx=0.002)))
    intgr4.append(np.log10(integrate.trapz(q,dx=0.002)))
    intgr5.append(np.log10(integrate.trapz(z,dx=0.002)))


# In[233]:


plt.figure(figsize=(10,10))
plt.plot(np.array(logN[0:]+np.log10(1608.45*0.0577)) ,intgr1[0:]-np.log10(1608.45),label=r"Fe II")
plt.plot(np.array(logN[0:]+np.log10(1526.7*0.133)) ,intgr2[0:]-np.log10(1526.7),label=r"Si II")
plt.plot(np.array(logN[0:]+np.log10(1334.53*0.1278)) ,intgr3[0:]-np.log10(1334.53),label=r"C II")
plt.plot(np.array(logN[0:]+np.log10(0.048*1304.2)) ,intgr4[0:]-np.log10(1304.2),label=r"O I")
plt.plot(np.array(logN[0:]+np.log10(1215.67*0.42)) ,intgr5[0:]-np.log10(1215.67),label=r"Lyman$\alpha$")
plt.title('Curve  of growth for different ions')
plt.ylabel(r'log[W/$\lambda$]')
plt.xlabel(r'Log[N $\lambda$ f]')
plt.legend()
plt.show()

