import glob
import astropy.io.fits as pf
import numpy as np
import scipy
import matplotlib.pyplot as plt
import import_functies
###############################################
save_data = False
save_data_fout = False
###############################################
show_data = False
show_soortelijke_weerstand = True
###############################################
lengte = 336.64
fout_lengte = 5
diameter = 0.0008
fout_diameter = 0.00002
###############################################

doorsnede = np.pi*(diameter/2)**2
fout_doorsnede = doorsnede*fout_diameter/diameter*np.sqrt(2)


def laad_data(file, kolommen=2):
    lijsten = []
    for i in range(kolommen):
        lijsten.append(np.loadtxt(file, usecols=i, unpack=True))
    return lijsten

def plot_data(x, y, fout_X=[None], fout_Y=[None], titel='title', x_label='x label', y_label='y label',Functie_naam='datapunten', X_naam='X naam', Y_naam='Y naam'):
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
    if fout_Y[0] != None:
        ax.errorbar(x, y, yerr=fout_Y, label="Fout op "+Y_naam, marker="o", markersize=0, fmt=" ",color="b", ecolor="black", capsize=2, capthick=2, linewidth=1)
    if fout_X[0] != None:
        ax.errorbar(x, y, xerr=fout_X, label="Fout op "+X_naam, marker="o", markersize=0, fmt=" ",color="b", ecolor="black", capsize=2, capthick=2, linewidth=1)
    ax.plot(x, y, '.', markersize=5, label=Functie_naam, color="black")
    ax.set_title(titel)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend()
    plt.show()



def fit_data(x, y, fout_X=[None], fout_Y=[None], titel='title', x_label='x label', y_label='y label',Functie_naam='datapunten', X_naam='X naam', Y_naam='Y naam'):

    def lin_model(x, rho_c, alpha):
        return rho_c * (1+alpha*(x-273.15))
    min_values = (0, 0)
    max_values = (1, 1)
    bounds = (min_values, max_values)
    fit_waarden_lin, fit_fout_lin_ = scipy.optimize.curve_fit(lin_model, x, y, bounds=bounds, maxfev=10000)
    fit_fout_lin = np.sqrt(np.diag(fit_fout_lin_))
    print('fit waarden lin: ' + str(fit_waarden_lin)+'+/-'+str(fit_fout_lin))
    def fitfuctie_lin(x):
        return fit_waarden_lin[0] * (1+fit_waarden_lin[1]*(x-273.15))

    def exp_model(x, rho_0, alpha):
        return rho_0 * np.exp(x * alpha)
    min_values = (0, 0)
    max_values = (1, 1)
    bounds = (min_values, max_values)
    fit_waarden_exp, fit_fout_exp_ = scipy.optimize.curve_fit(exp_model, x, y, bounds=bounds, maxfev=10000)
    fit_fout_exp = np.sqrt(np.diag(fit_fout_exp_))
    print('fit waarden exp: ' + str(fit_waarden_exp)+'+/-'+str(fit_fout_exp))
    def fitfuctie_exp(x):
        return fit_waarden_exp[0] * np.exp(x * fit_waarden_exp[1])

    temperaturen = np.linspace(x[0]-10, x[-1]+10, 1000)
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
    ax.plot(temperaturen, fitfuctie_lin(temperaturen), '-', markersize=1, label='lineaire fit', color="r")
    ax.plot(temperaturen, fitfuctie_exp(temperaturen), '-', markersize=1, label='exponentiële fit', color="g")
    if fout_Y[0] != None:
        ax.errorbar(x, y, yerr=fout_Y, label="Fout op "+Y_naam, marker="o", markersize=0, fmt=" ",color="r", ecolor="black", capsize=2, capthick=2, linewidth=1)
    if fout_X[0] != None:
        ax.errorbar(x, y, xerr=fout_X, label="Fout op "+X_naam, marker="o", markersize=0, fmt=" ",color="r", ecolor="black", capsize=2, capthick=2, linewidth=1)
    ax.plot(x, y, '.', markersize=5, label=Functie_naam, color="black")
    ax.set_title(titel)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.legend()
    plt.show()
    return fit_waarden_lin, fit_fout_lin


if __name__ == "__main__":
    if save_data:
        import_functies.insert_data()
    if save_data_fout:
        import_functies.insert_data_fout()

    temperatuur, weerstand = laad_data("data/data.txt")
    temperatuur_fout, weerstand_fout = laad_data("data/data_fout.txt")
    if show_data:
        plot_data(temperatuur, weerstand, temperatuur_fout, weerstand_fout, 'Meetresultaten', 'Temperatuur $[K]$', 'Weerstand $[\Omega ]$', 'Datapunten','temperatuur', 'weerstand')

    if show_soortelijke_weerstand:
        specifieke_weerstand = 10**(7)*weerstand*doorsnede/lengte
        specifieke_weerstand_fout = np.sqrt((weerstand_fout/weerstand)**2+(fout_lengte/lengte)**2+(fout_doorsnede/doorsnede)**2)*specifieke_weerstand
        print(specifieke_weerstand)
        print(specifieke_weerstand_fout)
        plot_data(temperatuur, specifieke_weerstand, temperatuur_fout, specifieke_weerstand_fout, 'Specifieke weerstand', 'Temperatuur $[K]$',r"$\rho $ x$ \, 10^{-7}$ $[\Omega m]$", 'Datapunten', 'temperatuur', 'weerstand')
        fit_data(temperatuur, specifieke_weerstand, temperatuur_fout, specifieke_weerstand_fout, 'Specifieke weerstand', 'Temperatuur $[K]$',r"$\rho $ x$ \, 10^{-7}$ $[\Omega m]$",'Datapunten','temperatuur', 'weerstand')


'''
temperatuur(weerstand):                                                        geeft de temperatuur voor een gegeven weerstand
insert_data():                                                                 maak datalijst met temperaturen en weerstanden
insert_data_fout():                                                            idem met fouten
laad_data("filepath", kolommen)                                                Geeft een lijst van arrays terug met daarin de data die in de kolommen staat
plot_data(X, Y, fout_X, fout_Y, titel, x_label, y_label, X_naam, Y_naam):      maakt een grafiek met de volgende (optionele) data
'''

'''
    def exp_model(x, rho_0, alpha):
        return rho_0 * np.exp(x * alpha)
    min_values = (0, 0)
    max_values = (0.02, 0.005)
    bounds = (min_values, max_values)
    fit_waarden_exp, fit_fout_exp = scipy.optimize.curve_fit(exp_model, x, y, bounds=bounds, maxfev=10000)
    print(fit_waarden_exp)
    def fitfuctie_exp(x):
        return fit_waarden_exp[0] * np.exp(x * fit_waarden_exp[1])
'''
