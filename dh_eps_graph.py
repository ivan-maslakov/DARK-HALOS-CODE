import numpy as np
import scipy
import matplotlib.pyplot as plt
E_min = -30
ro_0 = 10
r_0 = 1
G = 6.67*10**-11
G = 1
m = 1.78e-24
m = 1
k = 1.380649 * 10 ** (-23)
k = 1
T_out = 10**12 * 11604
T_out = 1
Cm = 1
#T_out = 5
#E_min = 10 * k * T_out
#ro_0 = 4
a = 1

def T_sht(E, M, r):
    vspp = 2/m * (E-U(r)) - M**2/m/r**2
    if vspp < 0:
        return 0
    return 1 / (vspp) ** 0.5


def ro(r):
    #return ro_0 * np.exp(-r/r_0)
    c = -1
    A = 1
    #return 1 / (c * r**2 + A * r ** 2 * np.log(r))
    #return 216 / (A**3 * r **6 - 6 * c * A**2 * r**4 + 12 * c**2 * A * r**2 - 8 * c**3)
    #return ro_0/(r/a)**1 / (1+r/a)**3
    return Cm * r **(-2.867) * T_out


def U(r):
    #return -4 * np.pi * G * ro_0 * r_0 ** 3 / r * (np.exp(-r / r_0) * ((r/r_0 + 1) ** 2 + 1) + 2)
    #return -4 * np.pi * G * (scipy.integrate.quad(lambda x: x ** 2 * ro(x), 0, r)[0] / r + scipy.integrate.quad(lambda x: x * ro(x), r, np.inf)[0])
    #a = 1
    #return -4 * np.pi * G * ro_0 * a ** 2 / (2 * (1 + r / a))
    return -4 * np.pi * G * Cm * T_out * r **(-0.867) * 0.61778
'''
def T_dif(r):
    return scipy.integrate.quad(lambda x: ,np.inf)[0]
'''

def eps0(E, T = T_out):
    if E < 0:
        return 0
    return 2 * np.pi / (np.pi * k * T) ** (2/3) * E**0.5 * np.exp(-E / k / T)

def sig0(E, T):
    return scipy.integrate.quad(lambda x: eps0(x, T), 0, E)[0]

def eps(E, T_out):
    if E > 0:
        return 0
    return scipy.integrate.quad(lambda x: sig0(E - U(x), T_out * (a/x)**0.867) * 4 * np.pi * x**2 * ro(x), a/10, np.inf)[0]

def F0(v,T):
    return scipy.integrate.quad(lambda x: 4 * np.pi * (m/2/np.pi/k/T)**1.5 * x**2 * np.exp(-m*x**2 / 2/k/T), 0, v)[0]

def theta(x):
    if x > 0:
        return 1
    return 0

def integrant_a(L, alf, T):
    return scipy.integrate.quad(lambda x: F0(L / m / x / np.sin(alf), T * (a/x)**0.867) * 4 * np.pi * x**2 * ro(x / m) * theta(-U(x) - L/2/x**2), 0, a*5)[0]

def mu(L, T_out):
    alfs = np.arange(0,np.pi/2, np.pi/40)
    dalf = alfs[1] - alfs[0]
    ans = 0
    for al in alfs:
        ans += integrant_a(L, al, T_out) * dalf
    #return 2 * scipy.integrate.quad(lambda al: integrant_a(L, al, T_out), 0, np.pi/2)[0]
    return ans

r_min, r_max = 0.0, 100
N_r = 500
r_range = np.arange(r_min, r_max, (r_max - r_min) / N_r)
dr = (r_max - r_min) / N_r

fig2, axes = plt.subplots(nrows=3, ncols=2, figsize=(9, 9))

ax0 = (axes.flatten())[1]
ax0_u = []
for ax0_r in r_range:
    ax0_u.append(U(ax0_r))
ax0.scatter(r_range, ax0_u, s = 5, color = 'b')
print('start')
print(ax0_u)
print('finish', r_range[1] - r_range[0])
ax0.set_title('U(r), r_0 = 1, ro_0 = 1')


ax1 = (axes.flatten())[0]
ax1.scatter(r_range, ro(r_range), s = 5, color = 'b')
ax1.set_title('ro(r), r_0 = 1, ro_0 = 1')
ax2 = (axes.flatten())[2]
ax3 = (axes.flatten())[3]
ax4 = (axes.flatten())[4]

for T_out in [1, 5, 15]:
    El = k * T_out * 50
    Ll = (El)**0.5 * 1
    E_range = np.arange(-El,0, El / 100)
    L_range = np.arange(0, Ll, Ll / 10)
    #L_range = []
    ax2_epss = []
    for ax2_e in -E_range:
        ax2_epss.append(eps0(ax2_e, T_out))
    ax2.plot(-E_range, ax2_epss, label = T_out)
    plt.grid()
    #ax2.set_title('eps0(T), r_0 = 1, ro_0 = 1')



    ax3_epss = []
    for ax3_e in E_range:
        print(ax3_e)
        ax3_epss.append(eps(ax3_e, T_out))
    ax3_epss = np.array(ax3_epss)
    ax3_epss = ax3_epss / ax3_epss.sum() / abs(E_range[1] - E_range[0])
    ax3.plot(E_range, ax3_epss, label = T_out)
    print(ax3_epss)
    plt.legend()
    plt.grid()
    #ax3.set_title('eps(E), r_0 = 1, ro_0 = 1')

    ax4_mus = []
    for ax4_l in L_range:
        print(ax4_l / Ll)
        ax4_mus.append(mu(ax4_l, T_out))
    ax4_mus = np.array(ax4_mus)
    ax4_mus = ax4_mus / ax4_mus.sum() / abs(L_range[1] - L_range[0])
    ax4.plot(L_range, ax4_mus, label=T_out)
    plt.legend()
    plt.grid()
    '''
    '''
plt.legend()
plt.show()