#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 23:37:28 2020

@author: caitlin
"""
from mpl_toolkits import mplot3d

import numpy as np
import matplotlib.pyplot as plt
import random as ran


# beta=(1/(k_B*T))
# t= k_B*T

#%%

#J_range=np.arange(0,0.1,0.1)  #interaction strength
J_range=np.array([0.5])         #case when J=0

for i in range(len(J_range)):
    J=J_range[i]

    def free_energy_1D(t,H,J):
        return -t*(J/t + np.log(np.cosh(H/t) + np.sqrt(np.sinh(H/t)**2 + np.exp(-4.*J/t))))

    t = np.linspace(0.1, 1, 30)
    h = np.linspace(-1, 1, 30)

    T, H = np.meshgrid(t, h)
    F = free_energy_1D(T, H, J)
    
    fig = plt.figure(1)
    ax = plt.axes(projection='3d')
    ax.contour3D(T, H, F, 50, cmap='binary')
    ax.plot_wireframe(T, H, F, color='black')
    ax.set_title('wireframe')
    ax.view_init(10, 140)
    ax.set_xlabel('t')
    ax.set_ylabel('H')
    ax.set_zlabel('F')


    ax.plot_surface(T, H, F, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    ax.set_title('Free energy');
    ax.view_init(10, 155)


#%%2D projection Free e vs mag field for different t

t = np.linspace(0.01, 1, 100)
h = np.linspace(-1, 1, 1000)

#free_energy_1D(t,H,J)
J=0

plt.figure(2)
plt.plot(h,free_energy_1D(0.001, h, J),label="t=0")
plt.plot(h,free_energy_1D(0.1, h, J),label="t=0.1")
plt.plot(h,free_energy_1D(0.5, h, J),label="t=0.5")
plt.plot(h,free_energy_1D(1., h, J),label="t=1.")
plt.xlabel("h")
plt.ylim(-1.2,0.0def susceptibility(t ,H, J):
    beta=(1/(t))
    chi = beta*(np.cosh(beta*H)*np.exp(-4*beta*J))/(((np.sinh(beta*H)**2)+np.exp(-4*beta*J))**1.5)
    return chi
5)
plt.ylabel("Free E")
plt.title("J="+str(J))

plt.legend()



#%% Vary J

t = np.linspace(0.01, 1, 100)
h = np.linspace(-1, 1, 100)

#free_energy_1D(t,H,J)
plt.figure(2)
plt.plot(h,free_energy_1D(t, h, 0),label="j=0.")
plt.plot(h,free_energy_1D(t, h, 0.1),label="j=0.1")
plt.plot(h,free_energy_1D(t, h, 0.3),label="j=0.3")
plt.plot(h,free_energy_1D(t, h, 0.5),label="j=0.5")
plt.plot(h,free_energy_1D(t, h, 1),label="j=1")
plt.xlabel("h")
plt.ylabel("Free E")
plt.legend()

#%%   MAGNETISATION

J=0.1

def magnetization(t, H, J):
    beta=(1/(t))
    m=np.sinh(beta*H)/(np.sqrt((np.sinh(beta*H))**2+np.exp(-4*beta*J)))
    return m

t = np.linspace(0.01, 1, 30)
h = np.linspace(-1, 1, 30)

T, H = np.meshgrid(t, h)
M = magnetization (T, H, J)

fig = plt.figure(3)
ax = plt.axes(projection='3d')
ax.contour3D(T, H, M, 50, cmap='binary')
ax.plot_wireframe(T, H, M, color='black')
ax.set_title('wireframe')
ax.view_init(10, 140)
ax.set_xlabel('t')
ax.set_ylabel('H')
ax.set_zlabel('M')

ax.plot_surface(T, H, M, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Magnetisation, J='+str(J));
ax.view_init(10, 230)

#t = np.linspace(0.01, 1, 100)
h = np.linspace(-1, 1, 100)

#free_energy_1D(t,H,J)
plt.figure(2)
plt.plot(h,magnetization(0.01, h, J),label="T=0.01")
plt.plot(h,magnetization(0.1, h, J),label="T=0.1")
plt.plot(h,magnetization(0.5, h, J),label="T=0.5")
plt.plot(h,magnetization(1, h, J),label="T=1")
plt.xlabel("h")
plt.ylabel("Magentisation")
plt.title("Magnetisation, J="+str(J))
plt.legend()


#%%   SUSCEPTIBILITY

J=0.1


t = np.linspace(0.01, 1, 30)
h = np.linspace(-1, 1, 300)


def susceptibility(t ,H, J):
    beta=(1/(t))
    chi = beta*(np.cosh(beta*H)*np.exp(-4*beta*J))/(((np.sinh(beta*H)**2)+np.exp(-4*beta*J))**1.5)
    return chi


T, H = np.meshgrid(t, h)
chi = susceptibility (T, H, J)

fig = plt.figure(4)
ax = plt.axes(projection='3d')
ax.contour3D(T, H, chi, 50, cmap='binary')
ax.plot_wireframe(T, H, chi, color='black')
ax.set_title('wireframe')
ax.view_init(10, 140)
ax.set_xlabel('t')
ax.set_ylabel('H')
ax.set_zlabel(r'$\chi$')

ax.plot_surface(T, H, chi, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Susceptibility, J='+str(J));
ax.view_init(10, 130)

plt.figure(2)
plt.plot(h,susceptibility(0.01, h, J),label="T=0.01")
plt.plot(h,susceptibility(0.1, h, J),label="T=0.1")
plt.plot(h,susceptibility(0.5, h, J),label="T=0.5")
plt.plot(h,susceptibility(1, h, J),label="T=1")
plt.ylim(0,5)


plt.xlabel("h")
plt.ylabel("Susceptibility")
plt.title("Susceptibility, J="+str(J))
plt.legend()


#%%Limiting cases

  #interaction strength




t = np.linspace(0.1, 10, 300)
#J_range=np.arange(0.01,1,0.1)

J=0.1

def susceptibility_red(t , J):
    beta=(1/(t))
    chi = beta*(np.exp(-4*beta*J))/((np.exp(-4*beta*J))**1.5)
    return chi

T, TOJ = np.meshgrid(t, t/J)

chi = susceptibility_red(T, J)

fig = plt.figure(5)
ax = plt.axes(projection='3d')
ax.contour3D(T, TOJ, T*chi, 50, cmap='binary')
ax.plot_wireframe(T, TOJ, T*chi, color='black')
ax.set_title('wireframe')
ax.view_init(10, 140)
ax.set_xlabel('t')
ax.set_ylim(0,10)
ax.set_ylabel('t/j')
ax.set_zlabel('t*chi')


ax.plot_surface(T, TOJ, T*chi, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('Reduced Susceptibility');
ax.view_init(10, 300)

    #%%
plt.figure(2)

plt.plot(t/0.1,t*susceptibility_red(t, 0.1),label="J=0.1")
plt.plot(t/0.5,t*susceptibility_red(t, 0.5),label="J=0.5")

plt.plot(t/1,t*susceptibility_red(t, 1),label="J=1")
plt.ylim(0,10)
plt.xlim(0,10)



plt.xlabel("t/J")
plt.ylabel(r"t$\chi(T,0)$")
plt.title("Reduced Susceptibility")
plt.legend()
   
#%% Correlation length



def cor_length(J,t):
    return -1/(np.log(np.tanh(J/t)))  

t = np.linspace(0.01, 20, 300)

plt.xlabel("t")
plt.ylabel(r"t$\xi(t,0)$")
#plt.plot(t, cor_length(0.01, t), label="J=0.01")
plt.plot(t, cor_length(0.1, t), label="J=0.1")
plt.plot(t, cor_length(1, t), label="J=1")
plt.plot(t, cor_length(5, t), label="J=5")

plt.ylim(0,15)
plt.legend()
plt.title("Correlation Length" )
