import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
from scipy.optimize import curve_fit
import pylab as pltl
####################################
show = True

write = True
####################################
dia_wiel = 0.66
massa_fietser = 115
fout_massa_fietser = 5
ro_lucht = 1.2
doorsnede_fieser = 0.57
fout_doorsnede_fietser = 0.25
helling = 2.5*np.pi/180
fout_helling = 0.25*np.pi/180
g = 9.81
#####################
# v(t) niet naukeurig omdat ricobepaling
# CR is heeeel gevoelig!!!


c1 = np.loadtxt('data/zonderv0.txt')
b1 = np.loadtxt('data/1bar.txt')
b2 = np.loadtxt('data/2bar.txt')
b3 = np.loadtxt('data/3bar.txt')
t1 = c1*10**(-3)
t2 = b1*10**(-3)
t3 = b2*10**(-3)
t4 = b3*10**(-3)
x1 = np.array(range(len(t1)))*np.pi*0.66
x2 = np.array(range(len(t2)))*np.pi*0.66
x3 = np.array(range(len(t3)))*np.pi*0.66
x4 = np.array(range(len(t4)))*np.pi*0.66

def model_x(x, v_e, tau):
    """logaritmische functie"""
    return v_e*tau*np.log((1+np.exp(2*x/tau))/((1+1)*np.exp(x/tau)))
popt1_x, pcov1_x = curve_fit(model_x, t1, x1)
popt2_x, pcov2_x = curve_fit(model_x, t2, x2)
popt3_x, pcov3_x = curve_fit(model_x, t3, x3)
popt4_x, pcov4_x = curve_fit(model_x, t4, x4)
fout1_x = np.sqrt(np.diag(pcov1_x))
fout2_x = np.sqrt(np.diag(pcov2_x))
fout3_x = np.sqrt(np.diag(pcov3_x))
fout4_x = np.sqrt(np.diag(pcov4_x))
labels = ("v_e[x]", "tau[x]")
units = ("", "")
if write == True:
    print("v_e, tau uit x(t):")
    for l, p, s, u in zip(labels, popt1_x, fout1_x, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    for l, p, s, u in zip(labels, popt2_x, fout2_x, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    for l, p, s, u in zip(labels, popt3_x, fout3_x, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    for l, p, s, u in zip(labels, popt4_x, fout4_x, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")

fig, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(5, 3))
ax.set_title('x(t) datapunten')
ax.plot(t1, x1, '.', markersize = 1, label = 'referentie')
ax.plot(t2, x2, '.', markersize = 1, label = '1 bar')
ax.plot(t3, x3, '.', markersize = 1, label = '2 bar')
ax.plot(t4, x4, '.', markersize = 1, label = '3 bar')
ax.set_xlabel('tijd [s]')
ax.set_ylabel('afstand [m]')
ax.legend()
if show == True:
    plt.show()

v1 = np.gradient(x1, t1)
v2 = np.gradient(x2, t2)
v3 = np.gradient(x3, t3)
v4 = np.gradient(x4, t4)
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(5, 3))
ax.set_title('v(t) datapunten')
ax.plot(t1, v1, '.', markersize=1, label='referentie')
ax.plot(t2, v2, '.', markersize=1, label='1 bar')
ax.plot(t3, v3, '.', markersize=1, label='2 bar')
ax.plot(t4, v4, '.', markersize=1, label='3 bar')
ax.set_xlabel('tijd [s]')
ax.set_ylabel('snelheid [m/s]')
ax.legend()
if show == True:
    plt.show()

def model_v(x, v_e, tau):
    """snelheid"""
    return v_e*((np.exp(2*x/tau)-1))/(np.exp(2*x/tau)+1)
popt1_v, pcov1_v = curve_fit(model_v, t1, v1)
popt2_v, pcov2_v = curve_fit(model_v, t2, v2)
popt3_v, pcov3_v = curve_fit(model_v, t3, v3)
popt4_v, pcov4_v = curve_fit(model_v, t4, v4)
fout_v = np.sqrt(np.diag(pcov1_v))
labels = ("v_e[v]", "tau[v]")
units = ("", "")
if write == True:
    print("v_e, tau uit v(t):")
    for l, p, s, u in zip(labels, popt1_v, fout_v, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    for l, p, s, u in zip(labels, popt2_v, fout_v, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    for l, p, s, u in zip(labels, popt3_v, fout_v, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    for l, p, s, u in zip(labels, popt4_v, fout_v, units):
        print(f"{l}: {p:.3f} ± {s:.3f} {u}")

#CD bepalen
def CD(popt):
    return (2*massa_fietser)/(ro_lucht*doorsnede_fieser*popt[0]*popt[1])
def fout_CD(popt, fouten):
    fout_tau = fouten[1]
    fout_ve = fouten[0]
    return CD(popt)/massa_fietser*np.sqrt(fout_massa_fietser**2 + massa_fietser**2*(fout_doorsnede_fietser/doorsnede_fieser)**2 + massa_fietser**2*(fout_tau**2)+massa_fietser**2*(fout_ve**2))
CD1_x = str(CD(popt1_x))+ '±' +str(fout_CD(popt1_x, fout1_x))
CD2_x = str(CD(popt2_x))+ '±' +str(fout_CD(popt2_x, fout2_x))
CD3_x = str(CD(popt3_x))+ '±' +str(fout_CD(popt3_x, fout3_x))
CD4_x = str(CD(popt4_x))+ '±' +str(fout_CD(popt4_x, fout4_x))
CD1_v = str(CD(popt1_v))+ '±' +str(fout_CD(popt1_v, fout_v))
CD2_v = str(CD(popt2_v))+ '±' +str(fout_CD(popt2_v, fout_v))
CD3_v = str(CD(popt3_v))+ '±' +str(fout_CD(popt3_v, fout_v))
CD4_v = str(CD(popt4_v))+ '±' +str(fout_CD(popt4_v, fout_v))
if write == True:
    print("CD:")
    print(CD1_x, CD2_x, CD3_x, CD4_x)
    print(CD1_v, CD2_v, CD3_v, CD4_v)

#CR bepalen
def CR(popt):
    return np.tan(helling) - popt[0]/(popt[1]*g*np.cos(helling))
def fout_CR(popt, fouten):
    return np.sqrt((popt[1]*g-np.sin(helling)*popt[0])**2*(fout_helling)**2 + fouten[0]**2 + popt[0]**2*(fouten[1]/popt[1])**2)/(popt[1]*g*np.cos(helling))
CR1_x = str(CR(popt1_x))+ '±' +str(fout_CR(popt1_x, fout1_x))
CR2_x = str(CR(popt2_x))+ '±' +str(fout_CR(popt2_x, fout2_x))
CR3_x = str(CR(popt3_x))+ '±' +str(fout_CR(popt3_x, fout3_x))
CR4_x = str(CR(popt4_x))+ '±' +str(fout_CR(popt4_x, fout4_x))
CR1_v = str(CR(popt1_v))+ '±' +str(fout_CR(popt1_v, fout_v))
CR2_v = str(CR(popt2_v))+ '±' +str(fout_CR(popt2_v, fout_v))
CR3_v = str(CR(popt3_v))+ '±' +str(fout_CR(popt3_v, fout_v))
CR4_v = str(CR(popt4_v))+ '±' +str(fout_CR(popt4_v, fout_v))
if write == True:
    print("CR:")
    print(CR1_x, CR2_x, CR3_x, CR4_x)
    print(CR1_v, CR2_v, CR3_v, CR4_v)
def plot_x_v():
    # data to plot
    n_groups = 4
    CR_x = (CR(popt1_x), CR(popt2_x), CR(popt3_x), CR(popt4_x))
    CR_v = (CR(popt1_v), CR(popt2_v), CR(popt3_v), CR(popt4_v))
    # create plot
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.35
    opacity = 0.8
    rects1 = plt.bar(index, CR_x, bar_width,
    alpha=opacity,
    color='b',
    label='$C_R$ uit $x$')
    rects2 = plt.bar(index + bar_width, CR_v, bar_width,
    alpha=opacity,
    color='g',
    label='$C_R$ uit $v$')
    plt.xlabel('Meting')
    plt.ylabel('Waarde $C_R$')
    plt.title('Verschil van $C_R$ uit x(t) en v(t)')
    plt.xticks(index + bar_width, ('referentie', '1 bar', '2 bar', '3 bar'))
    plt.legend()
    plt.tight_layout()
    plt.show()
if show == True:
    plot_x_v()

def plot_D_R():
    labels = ["1 bar", "2 bar", "3 bar", "referentie"]
    CD_x = (CD(popt2_x), CD(popt3_x), CD(popt4_x), CD(popt1_x))
    CR_x = (CR(popt2_v), CR(popt3_v), CR(popt4_v), CR(popt1_v))
    v_e_x = (popt2_x[0], popt3_x[0], popt4_x[0], popt1_x[0])
    fout_op_CD = (fout_CD(popt1_x, fout1_x), fout_CD(popt2_x, fout2_x), fout_CD(popt3_x, fout3_x), fout_CD(popt4_x, fout4_x))
    fout_op_CR = (fout_CR(popt1_x, fout1_x), fout_CR(popt2_x, fout2_x), fout_CR(popt3_x, fout3_x), fout_CR(popt4_x, fout4_x))
    fout_op_ve = (fout1_x[0], fout2_x[0], fout3_x[0], fout4_x[0])
    fig, axes = pltl.subplots(3, 1)
    ax1, ax2, ax3 = axes
    ax1.set_title('$C_D$, $C_R$ en $v_e$ per meting')
    ax1.barh(labels, CD_x,alpha=0.8, align="center",color='b')
    ax1.errorbar(CD_x, labels, markersize=3, xerr=fout_op_CD, fmt="o", color="dodgerblue")
    ax1.set_yticks(labels)
    ax1.set_xlabel('waarde van $C_D$')
    ax2.barh(labels, CR_x,alpha=0.8, align="center",color='b')
    ax2.errorbar(CR_x, labels, markersize=3, xerr=fout_op_CR, fmt="o", color="dodgerblue")
    ax2.set_yticks(labels)
    ax2.set_xlabel('waarde van $C_R$')
    ax3.barh(labels, v_e_x, alpha=0.8, align="center", color='g')
    ax3.errorbar(v_e_x, labels, markersize=1, xerr=fout_op_ve, fmt="o", color="lime")
    ax3.set_yticks(labels)
    ax3.set_xlabel('waarde van $v_e$ $[m/s]$')
    pltl.show()
if show == True:
    plot_D_R()

a1 = np.gradient(v1, t1)
a2 = np.gradient(v2, t2)
a3 = np.gradient(v3, t3)
a4 = np.gradient(v4, t4)
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(5, 3))
ax.set_title('a(t) datapunten')
ax.plot(t1, a1, '.', markersize=1, label='referentie')
ax.plot(t2, a2, '.', markersize=1, label='1 bar')
ax.plot(t3, a3, '.', markersize=1, label='2 bar')
ax.plot(t4, a4, '.', markersize=1, label='3 bar')
ax.set_xlabel('tijd [s]')
ax.set_ylabel('versnelling $[m/s^2]$')
ax.legend()
if show == True:
    plt.show()

W1_p = massa_fietser*scipy.integrate.cumulative_trapezoid(a1, x=x1)
W2_p = massa_fietser*scipy.integrate.cumulative_trapezoid(a2, x=x2)
W3_p = massa_fietser*scipy.integrate.cumulative_trapezoid(a3, x=x3)
W4_p = massa_fietser*scipy.integrate.cumulative_trapezoid(a4, x=x4)
W1 = max(W1_p)
W2 = max(W2_p)
W3 = max(W3_p)
W4 = max(W4_p)
t1_W = np.delete(t1, 0)
t2_W = np.delete(t2, 0)
t3_W = np.delete(t3, 0)
t4_W = np.delete(t4, 0)
for i in range(len(W1_p)):
    if W1_p[i] == W1:
        t_maxW1 = t1_W[i]
for i in range(len(W2_p)):
    if W2_p[i] == W2:
        t_maxW2 = t2_W[i]
for i in range(len(W3_p)):
    if W3_p[i] == W3:
        t_maxW3 = t3_W[i]
for i in range(len(W4_p)):
    if W4_p[i] == W4:
        t_maxW4 = t4_W[i]
if write == True:
    print("arbeid max:")
    print(W1, W2, W3, W4)
    print(t_maxW1, t_maxW2, t_maxW3, t_maxW4)
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(5, 3))
ax.set_title('W(t) datapunten')
ax.plot(t1_W, W1_p, '.', markersize=1, label='referentie')
ax.plot(t2_W, W2_p, '.', markersize=1, label='1 bar')
ax.plot(t3_W, W3_p, '.', markersize=1, label='2 bar')
ax.plot(t4_W, W4_p, '.', markersize=1, label='3 bar')
ax.set_xlabel('tijd [s]')
ax.set_ylabel('arbeid $[J]$')
ax.legend()
if show == True:
    plt.show()

for i in range(len(t1)):
    if t1[i] == t_maxW1:
        distance1 = x1[i]
        velocity1 = v1[i]
for i in range(len(t2)):
    if t2[i] == t_maxW2:
        distance2 = x2[i]
        velocity2 = v2[i]
for i in range(len(t3)):
    if t3[i] == t_maxW3:
        distance3 = x3[i]
        velocity3 = v3[i]
for i in range(len(t4)):
    if t4[i] == t_maxW4:
        distance4 = x4[i]
        velocity4 = v4[i]
h1 = distance1*np.sin(helling)
h2 = distance2*np.sin(helling)
h3 = distance3*np.sin(helling)
h4 = distance4*np.sin(helling)
E_pot1 = h1*massa_fietser*g
E_pot2 = h2*massa_fietser*g
E_pot3 = h3*massa_fietser*g
E_pot4 = h4*massa_fietser*g
E_kin1 = 0.5*massa_fietser*velocity1**2
E_kin2 = 0.5*massa_fietser*velocity2**2
E_kin3 = 0.5*massa_fietser*velocity3**2
E_kin4 = 0.5*massa_fietser*velocity4**2
fout_h1 = np.sqrt(distance1**2*np.cos(helling)**2*(fout_helling)**2+np.sin(helling)**2*0.01)
fout_h2 = np.sqrt(distance2**2*np.cos(helling)**2*(fout_helling)**2+np.sin(helling)**2*0.01)
fout_h3 = np.sqrt(distance3**2*np.cos(helling)**2*(fout_helling)**2+np.sin(helling)**2*0.01)
fout_h4 = np.sqrt(distance4**2*np.cos(helling)**2*(fout_helling)**2+np.sin(helling)**2*0.01)
fout_E_pot1 = g*np.sqrt(h1**2*fout_massa_fietser**2+massa_fietser**2*fout_h1**2)
fout_E_pot2 = g*np.sqrt(h2**2*fout_massa_fietser**2+massa_fietser**2*fout_h2**2)
fout_E_pot3 = g*np.sqrt(h3**2*fout_massa_fietser**2+massa_fietser**2*fout_h3**2)
fout_E_pot4 = g*np.sqrt(h4**2*fout_massa_fietser**2+massa_fietser**2*fout_h4**2)
fout_E_kin1 = np.sqrt((1/2*velocity1**2)**2*fout_massa_fietser**2+velocity1**2*fout1_x[0]**2*massa_fietser**2)
fout_E_kin2 = np.sqrt((1/2*velocity2**2)**2*fout_massa_fietser**2+velocity2**2*fout2_x[0]**2*massa_fietser**2)
fout_E_kin3 = np.sqrt((1/2*velocity3**2)**2*fout_massa_fietser**2+velocity3**2*fout3_x[0]**2*massa_fietser**2)
fout_E_kin4 = np.sqrt((1/2*velocity4**2)**2*fout_massa_fietser**2+velocity4**2*fout4_x[0]**2*massa_fietser**2)

if write == True:
    print("hoogte:", h1)
    print("fout h",fout_h1)
    print("fout Epot",fout_E_pot1)
    print("Epot",E_pot1)
    print(fout_E_kin1, E_kin1)



if write == True:
    print("Epot, Ekin:")
    print(E_pot1, E_pot2, E_pot3, E_pot4)
    print(E_kin1, E_kin2, E_kin3, E_kin4)

E_verloren1 = E_pot1 - E_kin1 - W1
E_verloren2 = E_pot2 - E_kin2 - W2
E_verloren3 = E_pot3 - E_kin3 - W3
E_verloren4 = E_pot4 - E_kin4 - W4
rendement1 = W1/E_pot1
rendement2 = W2/E_pot2
rendement3 = W3/E_pot3
rendement4 = W4/E_pot4
fout_op_rendement1 = rendement1*np.sqrt((fout_E_pot1/E_pot1)**2 + (fout_E_pot1/E_pot1)**2)
fout_op_rendement2 = rendement2*np.sqrt((fout_E_pot2/E_pot2)**2 + (fout_E_pot2/E_pot2)**2)
fout_op_rendement3 = rendement3*np.sqrt((fout_E_pot3/E_pot3)**2 + (fout_E_pot3/E_pot3)**2)
fout_op_rendement4 = rendement4*np.sqrt((fout_E_pot4/E_pot4)**2 + (fout_E_pot4/E_pot4)**2)



if write == True:
    print("verloren energie:")
    print(E_verloren1, E_verloren2, E_verloren3, E_verloren4)
    print("rendement:")
    print(rendement1,rendement2,rendement3,rendement4)
    print("fout rendement:")
    print(fout_op_rendement1, fout_op_rendement2, fout_op_rendement3, fout_op_rendement4)

labels = ["referentie", "3 bar", "2 bar", "1 bar"]
fig, axes = pltl.subplots(1, 1)
rendement = rendement1,rendement4,rendement3,rendement2
fout_op_rendement = fout_op_rendement1, fout_op_rendement2, fout_op_rendement3, fout_op_rendement4
axes.set_title('Rendement fiets per meting')
axes.bar(labels, rendement, alpha=1, align="center",color='tab:purple')
axes.set_xticks(labels)
axes.set_ylabel('rendement [%]')
axes.errorbar(labels, rendement, markersize=0, yerr=fout_op_rendement, fmt="o", color="r")
if show == True:
    pltl.show()