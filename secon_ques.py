#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import math
import scipy.special as s
import sys 
get_ipython().magic(u'matplotlib inline')
import mpld3
mpld3.enable_notebook()
import scipy.integrate as integrate 
from scipy.integrate import quad
c=3*10**5 


# In[2]:


import scipy.integrate as integrate 
from scipy.integrate import quad


def viogt(b,w,N):
        a=(3.03*10**-3)/b
        u=((1/w-1/1215.67)*c)/(b/1215.67)
        I=integrate.quad(lambda y:np.exp(-y**2)/(a**2 + (y-u)**2),-np.inf,np.inf)[0]
#        y=np.linspace(-1.0,1.0,1000000)
#        func=np.array(np.exp(-y**2)/(a**2 + (y-u)**2))  
#        I=integrate.trapz(func,y)
        return a/(b*np.sqrt(np.pi)) * 2.0*I * 1215.67 *N * f12*10**(-16)  

c=3*10**5               #in km/s
b=12                    # in km/s
f12 = 0.42





import os 
from astropy.table import Table
import matplotlib.pylab as plt

fptr= Table.read('hlsp_igm_hst_cos_1es1553_g130m-g160m_v3_spec.dat', format='ascii')
print(fptr.colnames)



plt.plot(fptr['col1'],fptr['col4'])
plt.ylim([1.58e-14, 1.61*10**(-14)])
plt.xlim([1160,1280])
plt.show()




# In[5]:


for i in range(1500,3700):
    fptr['col4'][i]=1.6e-14


# In[6]:


for i in range(2627,2727):
    fptr['col2'][i]=0.0


# In[7]:



plt.plot(fptr['col1'],fptr['col2'])
plt.show()


# In[8]:


global f_norm,wav
f_norm= np.divide(fptr['col2'],fptr['col4'])


# In[9]:


wav= fptr['col1']


# In[10]:


plt.plot(wav,f_norm)
plt.show()


# In[14]:


def eq_width(index):
    """
        calculates equivalent widths for detected lines
    """
    counter_1=index
    counter_2=index
    V=[]
    wav_v=[]
    Y=[]
    wav_y=[]
    while True:
        V.append(1.0-f_norm[counter_1])
        wav_v.append(wav[counter_1])
        if f_norm[counter_1]>0.95:
            sud=np.array(V)
            f_intgr= integrate.trapz(sud[0:-1],x=np.array(wav_v[0:-1]))
            break                
        counter_1+=1 
        
    while True:
        Y.append(1.0-f_norm[counter_2])
        wav_y.append(wav[counter_2])
        if f_norm[counter_2]>0.95:
            yud=np.array(sorted(Y))
            s_intgr = integrate.trapz(yud[0:-1],x=np.array(sorted(wav_y[0:-1])))
            break
        counter_2-=1
#    print(sud[0:-1],yud[0:-1])    
    return f_intgr + s_intgr   
                                      


# In[15]:


print eq_width(579)


# In[16]:


def wave2index(l):
    i=max(min(int((l-1135.43)/0.03),20117),0)
    while abs(wav[i]-l)>0.02:
        if(wav[i]<l):
            i=i+1
        else:
            i=i-1
    return i


# In[17]:


wave2index(1223.5)


# In[18]:


counter_line = 0 
Flag=0
indi=152
for y in f_norm[152:-1]:
    if Flag==0:
        indi+=1
        if y <= 0.5:
            Flag=1
            width = eq_width(indi)
            counter_line += 1 
            print width ,wav[indi]
            
    else :
        indi+=1
        if y>0.95:
            Flag=0    


# In[29]:


plt.plot(wav[2409:2944],f_norm[2409:2944],label=r'lyman $\alpha$ spectrum data')
viogt = np.vectorize(viogt)
w_range =np.linspace(1207.5, 1223.5,1000)
plt.plot(w_range,np.exp(-viogt(30,w_range,10**21.45)),'r',label='viogt profile fit with b=30 LogN=21.45')
plt.xlabel('wavelength( A)')
plt.ylabel('normalized flux')
plt.legend()    
plt.show()


# In[19]:


gptr= Table.read('absortion_lines.dat', format='ascii')
print(gptr.colnames)


# In[20]:


W_lines=[]
for w in gptr['col2']:
    W_lines.append(eq_width(wave2index(w)))


# In[21]:


print W_lines


# In[22]:


import csv
with open('linewidth.csv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(gptr['col1'],gptr['col2'],W_lines))


# In[23]:


aptr= Table.read('line_data.dat', format='ascii')
print(aptr.colnames)


# In[24]:


import scipy.integrate as integrate 
from scipy.integrate import quad

def Viogt(g,w,N,w_r,F):
        a=(g*w_r)/(b*8*np.pi)
        u=((1/w-1/w_r)*c)/(b/w_r)
        I=integrate.quad(lambda y:np.exp(-y**2)/(a**2 + (y-u)**2),-np.inf,np.inf)[0]
#        y=np.linspace(-1.0,1.0,1000000)
#        func=np.array(np.exp(-y**2)/(a**2 + (y-u)**2))  
#        I=integrate.trapz(func,y)

        return a/(b*np.sqrt(np.pi)) * 2.0*I * w_r *N * F*10**(-16) 


# In[36]:



def plotwave(l,N):
    counter=0
    for w in np.array(aptr['col2']):    
        if l==int(w):
            break
        counter+=1     
    l=w
    F=aptr['col3'][counter]
    Y=aptr['col4'][counter]
    w_range =np.linspace(l-2,l+2,1000)
    plt.plot(w_range,[np.exp(-Viogt(Y*10**(-10),e,10**N,l,F)) for e in w_range])
    Varray=[np.exp(-Viogt(Y*10**(-8),e,10**N,l,F)) for e in w_range]
    I1=wave2index(l-2)
    I2=wave2index(l+2)
    wl=np.add(0.12,wav[I1:I2])
    plt.plot(wl,f_norm[I1:I2])    
    plt.show()
    return np.min(f_norm[I1:I2]),np.min(Varray),F,Y
    


# In[39]:


b=31
plotwave(1144,15.5)


# In[33]:


'''
b=30
r=1144
for i in range(154,160):
    a,c,f_,y_ = plotwave(r,i/10.0)
    if abs(a-c)<0.1:
        break
w_range =np.linspace(r-2,r+2,1000)
plt.plot(w_range,np.exp(-Viogt(y_*10**(-8),w_range,10**(i/10.0),r,f_)))
I1=wave2index(r-2)
I2=wave2index(r+2)
wl=np.add(0.12,wav[I1:I2])
plt.plot(wl,f_norm[I1:I2])    
plt.show()



# In[31]:


b=30
r=1144
i=155
w_range =np.linspace(r-2,r+2,1000)
plt.plot(w_range,np.exp(-Viogt(y*10**(-10),w_range,10**(i/10.0),r,0.083)))
I1=wave2index(r-2)
I2=wave2index(r+2)
wl=np.add(0.12,wav[I1:I2])
plt.plot(wl,f_norm[I1:I2])    
plt.show()


# In[35]:


print(f_)


# In[ ]:



'''
