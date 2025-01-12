import plotly.graph_objects as go
import numpy as np
import scipy
import os
import scipy
import random
from scipy import linalg
import matplotlib.pyplot as plt

G = 1
ro_0 = 1
a = 1
#R = 0.001

def ro(r):
    #return ro_0 * np.exp(-r/r_0)
    c = -1
    A = 1
    #return 1 / (c * r**2 + A * r ** 2 * np.log(r))
    #return 216 / (A**3 * r **6 - 6 * c * A**2 * r**4 + 12 * c**2 * A * r**2 - 8 * c**3)
    return ro_0/(r/a)**1 / (1+r/a)**3

def Hernquist_U(r):
  return -4 * np.pi * G * ro_0 * a ** 2 / (2 * (1 + r / a))

def Tsh(r, Ee, Le):
  val = 2 * Ee - 2 * Hernquist_U(r) - Le ** 2 / r **2
  #val[val < 0] = np.inf
  if val <= 0:
    return 0
  return 1 / np.sqrt(val)

'''
def f(E, L):
  if 2 * E - 2 * Hernquist_U(R) - L ** 2 / R **2 < abs(Hernquist_U(R)) / 10 * 0:
#return 0
  inte, err = scipy.integrate.quad(lambda x: Tsh(x, E, L), 0, np.inf)
'''
E0 = 3
L0 = 1
i1 = 0
i2 = 150
def target(args):
    E0 = args[0]
    L0 = args[1]
    norm_ver_l = scipy.integrate.quad(lambda x: x**2 * np.exp(-x**2/L0**2), 0, np.inf)[0]
    ro_ar = []
    for R in np.arange(0.01,10,0.01)[i1:i2]:
        #print(R *10)
        E_min = Hernquist_U(R) * 1.1
        L_max = (-Hernquist_U(R) * 2 * R ** 2) ** 0.5 * 1.1

        de, dl = abs(-E_min) / 50, abs(L_max / 50)
        E_r = np.arange(E_min, 0, de)
        L_r = np.arange(0, L_max, dl)
        file_name = 'B_R=' + str(R) + '.npy'
        arr = np.load(file_name)
        ver_e = np.array([np.exp(E_r / E0) / E0 * de])
        #ver_e = (E_r - E0) / 3*2 / E0**2

        #print('inf')
        #ver_e[ver_e > -2/3/E0] = 0
        #ver_e[ver_e < 0] = 0
        #ver_e = np.array([ver_e])
        #ver_l = np.array([L_r ** 2 * np.exp(-L_r ** 2 / L0 ** 2) / norm_ver_l * dl])
        #ver_l = np.array([np.exp(-L_r / L0) / L0 * dl])
        ver_l = np.array([2 * np.pi / (np.pi * L0) ** 1.5 * L_r**0.5 * np.exp(-L_r / L0)])
        ver = ver_l.T @ ver_e
        ro_ar.append((arr * ver).sum() /4 / np.pi / R**2)
    #ret = ((abs((np.array(ro_ar) - ro(np.arange(0.01,10,0.01)[i1:i2]))) / np.array(ro_ar)).sum())/(i2-i1)
    #ret = ((abs((np.array(ro_ar) - ro(np.arange(0.01,10,0.01)[i1:i2])))) * (np.arange(0.01,10,0.01)[i1:i2])**2).sum() * 0.01
    ret = 0
    for ik in [10,20,30,40,50,70,100,200,300,400,500,600,700]:
        delta_M = (((np.array(ro_ar)[i1:ik] - ro(np.arange(0.01, 10, 0.01)[i1:ik]))) * (np.arange(0.01, 10, 0.01)[i1:ik]) ** 2).sum() * 0.01
        M = (((ro(np.arange(0.01, 10, 0.01)[i1:ik]))) * (np.arange(0.01, 10, 0.01)[i1:ik]) ** 2).sum() * 0.01
        ret += abs(delta_M / M)
    if np.isnan(ret):
        return 10000000000
    return ret

res2 = scipy.optimize.minimize(target, np.array([random.randint(0,4000)/100,random.randint(0,4000)/100]), method= 'COBYLA', bounds=[(0,np.inf), (0,np.inf)])


def ploter(args):
    E0 = args[0]
    L0 = args[1]
    #norm_ver_l = scipy.integrate.quad(lambda x: x**2 * np.exp(-x**2/L0**2), 0, np.inf)[0]
    ro_ar = []
    for R in np.arange(0.01,10,0.01)[0:750]:
        #print(R *10)
        E_min = Hernquist_U(R) * 1.1
        L_max = (-Hernquist_U(R) * 2 * R ** 2) ** 0.5 * 1.1


        de, dl = abs(-E_min / 50), abs(L_max / 50)
        E_r = np.arange(E_min, 0, de)
        L_r = np.arange(0, L_max, dl)
        file_name = 'B_R=' + str(R) + '.npy'
        arr = np.load(file_name)
        ver_e = np.array([np.exp(E_r / E0) / E0 * de])
        #print(E_r)
        #print()
        #print(L_r)
        #ver_e = (E_r - E0) / 3 * 2 / E0 ** 2
        #ver_e[ver_e > -2 / 3 / E0] = 0
        #ver_e[ver_e < 0] = 0
        #ver_e = np.array([ver_e])
        #ver_l = np.array([L_r ** 2 * np.exp(-L_r ** 2 / L0 ** 2) / norm_ver_l * dl])
        #ver_l = np.array([np.exp(-L_r / L0) / L0 * dl])
        ver_l = np.array([2 * np.pi / (np.pi * L0) ** 1.5 * L_r ** 0.5 * np.exp(-L_r / L0)])
        ver = ver_l.T @ ver_e
        ro_ar.append((arr * ver).sum() / 4 / np.pi / R**2)

    #plt.plot(np.arange(0.01,10,0.01)[0:750], (ro_ar-ro(np.arange(0.01, 10, 0.01)[0:750])) / ro(np.arange(0.01, 10, 0.01)[0:750]), color = 'r', label = 'delt')
    plt.plot(np.arange(0.01, 10, 0.01)[0:750], ro(np.arange(0.01, 10, 0.01)[0:750]), color='b', label='start')
    plt.plot(np.arange(0.01, 10, 0.01)[0:750],ro_ar, color='g',label='our')
    for i2 in [10,50,100,200,300,400,700]:
        delta_M = (((np.array(ro_ar)[i1:i2] - ro(np.arange(0.01, 10, 0.01)[i1:i2]))) * (np.arange(0.01, 10, 0.01)[i1:i2]) ** 2).sum() * 0.01
        M = (((ro(np.arange(0.01, 10, 0.01)[i1:i2]))) * (np.arange(0.01, 10, 0.01)[i1:i2]) ** 2).sum() * 0.01
        print('на', i2 / 100,'a: dm / m = ', delta_M / M)


print(res2)
ploter(res2.x)
#ploter([1,1])
plt.legend()
plt.show()