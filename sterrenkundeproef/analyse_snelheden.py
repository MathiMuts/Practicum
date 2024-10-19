import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
##############################################################################################
"""
In dit programma wordende snelheden geanaliseerd om daaruit de periode,
snelheid en massa(functie) te berekenen.
Deze worden geschreven in 'calculated data/%s/analyse' % sternaam.
De grootst mogelijke onzekerheid op v is op papier bepaald en is dan len(snelheden)
aantal keer in onzekerheid_snelheden.txt geplaatst
"""
##############################################################################################
show = True                     #toon grafieken
##############################################################################################
sternaam = 'BD+39_4926'         # 'BD+39_4926'   'BD-11_162'     'IRAS06165+3158'    'SS_Lep'
c = 299792458
G = 6.67259 * (10**(-11))
massa_zon = 1.989*(10**30)
massa_witte_dwerg = massa_zon*0.6
##############################################################################################

#De data inladen
dagen_file = 'calculated data/%s/dagen.txt' % sternaam
onzekerheid_snelheid_file = 'calculated data/%s/onzekerheid_snelheid.txt' % sternaam
snelheden_file = 'calculated data/%s/snelheden.txt' % sternaam
analyse = 'calculated data/%s/analyse.txt' % sternaam
dagen_original = np.loadtxt(dagen_file, usecols=0, unpack=True)         #array with days
dagen = dagen_original-dagen_original[0]
snelheden = np.loadtxt(snelheden_file, usecols=0, unpack=True)          #array with speeds
onzekerheid_snelheden = np.loadtxt(onzekerheid_snelheid_file, usecols=0, unpack=True)

#Deze funcite fit een sinus op snelheid i.f.v. dagen
# om daaruit de frequentie, amplitude, beginfase en verschuiving terug te geven.
def fit_it(min_amp, max_amp, min_hoek_v, max_hoek_v, min_beginfase, max_beginfase, min_versch, max_versch):                              #geef dataset en krijg sinus met bounds: min values, max values voor fitten
    def model(x, a, b, c, d):
        """Harmonisch model"""
        return a * np.sin(b*x + c) + d

    min_values = (min_amp, min_hoek_v, min_beginfase, min_versch)
    max_values = (max_amp, max_hoek_v, max_beginfase, max_versch)
    bounds = (min_values, max_values)

    popt, pcov = curve_fit(model, dagen, snelheden, bounds=bounds, maxfev=10000, sigma=onzekerheid_snelheden)

    fout = np.sqrt(np.diag(pcov))
    labels = ("amplitude", "hoeksnelheid", "beginfase", "verschuiving")
    units = (" ", " ", " ", " ")
    snelheden_km = snelheden/1000
    onzekerheid_snelheden_km = onzekerheid_snelheden/1000
    d = np.linspace(dagen[0], dagen[-1], 10000)
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
    ax.set_title("Periodische beweging van "+ sternaam)
    ax.plot(dagen, snelheden_km, '.', markersize = 1, label = 'Snelheid')
    ax.errorbar(dagen, snelheden_km, yerr=onzekerheid_snelheden_km, label="data", marker="o", markersize=0, fmt=" ", color="r", ecolor="Black", capsize=2, capthick=0.6,linewidth=0.6)
    ax.plot(d, (popt[0]*np.sin(popt[1]*d+popt[2])+popt[3])/1000, '.', markersize = 1, label = 'Snelheid')
    ax.set_xlabel('Dag')
    ax.set_ylabel('Snelheid [km/s]')
    plt.show()
    return labels, popt, fout, units

#Dit onderdeel geeft de geschatte parameters mee om een zinvolle sinus-fit terug te krijgen per ster
if sternaam == 'BD+39_4926':
    labels, waarde, fout, units = fit_it(10000, 20000,0.005, 0.01, 0, 7, -50000, -20000)
elif sternaam == 'BD-11_162':
    labels, waarde, fout, units = fit_it(8000, 15000,0.005, 0.02, 0, 7, 0, 10000)
elif sternaam == 'IRAS06165+3158':
    labels, waarde, fout, units = fit_it(15000, 25000,0.001, 0.05, 0, 7, -20000, -10000)
elif sternaam == 'SS_Lep':
    labels, waarde, fout, units = fit_it(15000, 25000,0.001, 0.05, 0, 7, 15000, 25000)

#geeft betekenisvolle namen aan de variabelen
amplitude, hoeksnelheid, beginfase, verschuiving = waarde
fout_amplitude, fout_hoeksnelheid, fout_beginfase, fout_verschuiving = fout
hoeksnelheid = hoeksnelheid
fout_hoeksnelheid = fout_hoeksnelheid

#Periode en massafunctie berekenen
periode = (2*np.pi/hoeksnelheid)*24*60*60
fout_periode = periode*fout_hoeksnelheid/hoeksnelheid
massafunctie = (periode*amplitude**3)/(2*np.pi*G)
fout_massafunctie = np.sqrt(((amplitude**3)/(2*np.pi*G))**2*fout_periode**2+(3*periode*amplitude**2/(2*np.pi*G))**2*fout_amplitude**2)

#berekende data printen
for l, p, s, u in zip(labels, waarde, fout, units):
    print(f"{l}: {p:.3f} ± {s:.3f} {u}")
print("periode: "+str(periode)+" ± "+str(fout_periode))
print("massafunctie: " + str(massafunctie/massa_zon)+" ± "+str(fout_massafunctie/massa_zon))

#De massa van de onzichtbare component plotten i.f.v. de hellingshoek
massas = []
fout_massas = []
ies = [i*np.pi/180 for i in range(10, 91)]
ies_degr = [i for i in range(10, 91)]
for i in ies:
    d = massafunctie/(np.sin(i)**3)
    coefficients = [1, -d, -2*d*massa_witte_dwerg, -d*massa_witte_dwerg**2]
    roots = np.roots(coefficients)
    roots_real = roots[np.isreal(roots)].real
    massas.append(roots_real[0])
    teller = 3 * roots_real[0] ** 2 * (np.sin(i)) ** 3 * (massa_witte_dwerg + roots_real[0]) ** 2 - roots_real[0] ** 3 * (np.sin(i)) ** 3 * 2 * (
                massa_witte_dwerg + roots_real[0])
    noemer = (massa_witte_dwerg + roots_real[0]) ** 4
    fout_massa = np.sqrt((fout_massafunctie)**2 * ((teller)/(noemer))**(-2))
    fout_massas.append(fout_massa)
massas = [i/massa_zon for i in massas]
fouten_massa = [i/massa_zon for i in fout_massas]
if show == True:
    fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
    ax.set_title('Massa van ' + sternaam+' i.f.v. de inclinatiehoek')
    ax.plot(ies_degr, massas, '.', markersize = 2.5, label = 'massa')
    ax.errorbar(ies_degr, massas, yerr=fouten_massa, label="data", marker="o", markersize=0, fmt=" ", color="r", ecolor="Black", capsize=1, capthick=0.6,linewidth=0.6)
    ax.set_xlabel('inclinatiehoek')
    ax.set_ylabel('zonnemassas')
    plt.show()

#Alle berekende data in analyse.txt schrijven
with open(analyse, 'w') as f:
    for l, p, s, u in zip(labels, waarde, fout, units):
        f.write(f"{l}: {p:.3f} ± {s:.3f} {u}" + "\n")
    f.write("periode: "+str(periode)+" ± "+str(fout_periode) + "\n")
    f.write("massafunctie: " + str(massafunctie/massa_zon)+" ± "+str(fout_massafunctie/massa_zon))



