import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
dia_wiel = 0.66
massa_fietser = 115
fout_massa_fietser = 5
ro_lucht = 1.2
doorsnede_fieser = 0.57
fout_doorsnede_fietser = 0.25
helling = 2.5
fout_helling = 1
g = 9.81
def gettimes(file):
    T_list = []
    T_list_mili = np.loadtxt("data/"+file+".txt", dtype=int)
    for T in T_list_mili:
        T_list.append(T/1000)
    return T_list

def displacement(T_list):
    X_list = []
    for i in range(len(T_list)):
        X = 0.66*np.pi*i
        X_list.append(X)
    return X_list

def velocity(T_list):
    V_list = []
    X_list = displacement(T_list)
    for i in range(len(T_list)):
        if T_list[i] != 0:
            V = X_list[i]/T_list[i]
            V_list.append(V)
        else:
            V_list.append(0)
    return V_list
def plotit(T_list, arg_list, title, as_x, as_y):
    T = T_list
    y = arg_list
    plt.title(title)
    plt.xlabel(as_x)
    plt.ylabel(as_y)
    plt.plot(T, y, "ob")
    plt.show()

def model(t, v_e, tau):
    """v(t) model"""
    return v_e*(np.exp(2*t/tau)-1)/(np.exp(2*t/tau)+1)
def fit(X, Y, dY, show):
    _, ax = plt.subplots(nrows=1, ncols=1, dpi=120, figsize=(5, 3))
    # We genereren een linspace om het model te plotten
    x = np.linspace(np.min(X), np.max(X), 100)

    # We schalen ook de stroom voor een beter leesbare x-as
    ax.errorbar(X, Y, yerr=dY, label="data",

                # De errorbars wat mooier maken :)
                marker="o", markersize=4, fmt=" ", color="black", ecolor="black", capsize=2, capthick=0.6,
                linewidth=0.6)

    popt, pcov = curve_fit(model, X, Y, sigma=dY, absolute_sigma=True)

    fout = np.sqrt(np.diag(pcov))
    labels = ("v_e", "tau")
    units = ("m/s", "s")

    # label,parameter, standaarddeviatie, unit
    # zip: arrays samenvoegen
    #for l, p, s, u in zip(labels, popt, fout, units):
        #print(f"{l}: {p:.3f} ± {s:.3f} {u}")
    # Schatting van de fitparametes
    ax.plot(x, model(x, popt[0], popt[1]), 'b--', label="model")
    #ax.plot(x, model(20, x, 20), 'b--', label="model")

    ax.set_ylabel("$V_{k}$ [V]")
    ax.set_xlabel("$I$ [mA]")
    ax.legend()
    plt.tight_layout();
    if show == True:
        plt.show()
    return popt[0], popt[1]


def C_D(m, ro, A, v_e, tau):
    return (2*m)/(ro*A*v_e*tau)

def C_R(v_e, tau, theta):
    return -(v_e/(tau*g)+np.sin(theta))/np.cos(theta)



##################################################################################################

##############################
file = "3bar"
##############################

T_list = gettimes(file)
X_list = displacement(gettimes(file))
V_list = velocity(gettimes(file))
dX = [0.1 for X in X_list]
dV = [0.011 for V in V_list]
#plotit(T_list, X_list, "title", "x-as", "y-as")
v_e = fit(T_list, V_list, dV, False)[0]
tau = fit(T_list, V_list, dV, False)[1]


