# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 15:34:10 2021

@author: xmu_xu
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#data=np.concatenate((np.random.normal(0.2,.05,2500),np.random.normal(0.7,.2,5000)))
data = np.loadtxt(open("new1_data.csv","rb"),delimiter=",",skiprows=0)
data = np.array(data)

#print(data)

y,x,_=plt.hist(data, 100, density=True, alpha=.4, facecolor='green', label='Density of data')

x=(x[1:]+x[:-1])/2 # for len(x)==len(y)

def gauss(x,mu,sigma,A):

    #return A*np.exp(-(x-mu)**2/2/np.pi/sigma**2)
    return A*np.exp(-(x-mu)**2/2/sigma**2)/np.sqrt(2*np.pi*sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):

    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

expected=(0.05,.01,10,0.7,.01,8) 

params,cov=curve_fit(bimodal, x, y, expected, method='trf', maxfev=50000) # // method lm trf or dogbox

sigma=np.sqrt(np.diag(cov))
#plt.plot(x,bimodal(x,*params),'r--',lw=2.5,label='Guassin Fitting')
#plt.legend(loc='best')
#plt.show()

print(params)

xf_x = []
xf_y = []
xf_y1 = []

xf_ave = params[0]
xf_sigma = params[1]
xf_A = params[2]

xf_ave1 = params[3]
xf_sigma1 = params[4]
xf_A1 = params[5]

xf_p = 0.0
xf_p1 = 0.0

for i in range(100):
    xf_b = x[i]
    xf_x.append(xf_b)
    xf_c = xf_A*np.exp(-1*(xf_b-xf_ave)**2/(2*xf_sigma**2))/np.sqrt(2*np.pi*xf_sigma**2)
    xf_c1 = xf_A1*np.exp(-1*(xf_b-xf_ave1)**2/(2*xf_sigma1**2))/np.sqrt(2*np.pi*xf_sigma1**2)
    xf_y.append(xf_c)
    xf_y1.append(xf_c1)
    
    xf_p += 0.01*xf_c
    xf_p1 += 0.01*xf_c1

xf_pp = xf_p/(xf_p+xf_p1)
xf_pp1 = xf_p1/(xf_p+xf_p1)

plt.plot(xf_x, xf_y,'r--', linewidth=2.5, label='pA=%6.3f' % xf_pp)
plt.plot(xf_x, xf_y1,'b-', linewidth=2.5, label='pN=%6.3f' % xf_pp1)
plt.savefig('fitting.svg', bbox_inches='tight', dpi=300)
plt.legend(loc='best')
plt.show()