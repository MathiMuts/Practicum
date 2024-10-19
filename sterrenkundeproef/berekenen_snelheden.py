import glob
import astropy.io.fits as pf
import numpy as np
import scipy
import matplotlib.pyplot as plt
##############################################################################################
"""
In dit programma wordende snelheden en bijhorende dagen berekend voor elke dataset.
Deze worden geschreven in 'calculated data/%s/snelheden' % sternaam.
De berekening van v gebeurt a.d.h.v. redshift (blueshift)
van het gemeten spectrum met de verwachtte data uit een labo.
"""
##############################################################################################
show = True                #toon grafieken
bereken = False             #True = voer de berekening van v uit (+/-2u)
##############################################################################################
sternaam = 'BD+39_4926'     # 'BD+39_4926'   'BD-11_162'     'IRAS06165+3158'    'SS_Lep'
c = 299792458
aantal_verfijningen_berekening_v = 10   #aantal verfijningen van het optimalisatieprobleem
verfijn_op_interval = (-100000, 100000) #het interval waarop we v zoeken
##############################################################################################

#De data inladen
datalist = glob.glob("datasets/%s/*.fits" % sternaam)
datalist.sort()
inputfile = "datasets/%s/LijnLijst.txt" % sternaam
lines, weights = np.loadtxt(inputfile, usecols=(0, 1), unpack=True)     #spectraallijnen, gewicht

#Over alle datasets herhalen
for data_index in datalist[:1]:
    spec = pf.getdata(data_index)                                       #intensiteit
    header = pf.getheader(data_index)                                   #golflengte, datum in julieaanse dagen, $unitx$


    # De golflengtes en header vergaren
    def golflengte_grid(spec, header):

        # Lees meta-data
        ref_pix = int(header['CRPIX1'])-1   # index van de referentiepixel
        ref_val = float(header['CRVAL1'])   # ln(golflengte) van de referentiepixel
        ref_del = float(header['CDELT1'])   # breedte van de pixel in eenheid van CRVAL1
        JD = header['BJD']  # datum van de waarneming in Juliaanse dagen
        unitx = header['CTYPE1']

        numberpoints = spec.shape[0]

        # Maak een golflengtegrid aan
        wavelengthbegin = ref_val - ref_pix*ref_del
        wavelengthend = wavelengthbegin + (numberpoints-1)*ref_del
        wavelengths = np.linspace(wavelengthbegin,wavelengthend,numberpoints)
        wavelengths = np.exp(wavelengths)

        return wavelengths, JD, unitx
    wavelengths, JD, unitx = golflengte_grid(spec, header)

    #plot: intensiteit(wavelengths)
    if show == True:
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
        ax.set_title('Intensiteit per golflengte')
        ax.plot(wavelengths, spec, '.', markersize = 1, label = 'intensiteit')
        ax.set_xlabel('Golflengte $[Å]$')
        ax.set_ylabel('Intensiteit $[J/m^2]$')
        ax.legend()
        plt.tight_layout()
        plt.show()
        #plt.savefig('Spectrum.png')

    #plot: weights(lines)
    if show == True:
        fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
        ax.set_title('weights(lines)')
        ax.plot(lines, weights, '.', markersize = 1, label = 'intensiteit')
        ax.set_xlabel('wavelengths')
        ax.set_ylabel('intensiteit')
        ax.legend()
        plt.tight_layout()
        plt.show()
        #plt.savefig('Spectrum.png')

    #interpolate intensiteit(wavelengths)
    g = scipy.interpolate.interp1d(wavelengths, spec)

    #deze functie bepaalt v voor elke dataset
    def bereken_v():

        #Geeft de waargenomen intensiteit van een golflengte
        def intensiteit(Lambda_shifted):
            return g(Lambda_shifted)

        #Geeft de som (rekening houdend met gewichen) van alle intensiteiten voor alle golflengtes in lijnlijst
        def som_doppler(v):
            return sum(weights[i]*intensiteit(lines[i] * (v / c + 1)) for i in range(len(lines)))

        #Optimaliseert de som naar een zo klein mogelijke waarde
        def optimaliseet_het(start, stop, step):
            mogelijke_v = np.linspace(start, stop, step)
            som_intensiteit = [som_doppler(v) for v in mogelijke_v]
            v_min = mogelijke_v[som_intensiteit.index(min(som_intensiteit))]
            if show == True:
                fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
                ax.set_title('snelheid: ' + str(start) + ", " + str(stop))
                ax.plot(mogelijke_v, som_intensiteit, '.', markersize=1, label='snelheid')
                ax.set_xlabel('mogelijke_v')
                ax.set_ylabel('som_intensiteiten')
                ax.legend()
                plt.show()
            return v_min

        #bepaalt het nieuwe interval om v op te zoeken
        def grenzen(oude_grenzen):
            interval = np.absolute(oude_grenzen[1]-oude_grenzen[0])
            v_temp = optimaliseet_het(oude_grenzen[0], oude_grenzen[1], 100)
            return (v_temp-interval/10, v_temp+interval/10)

        #voert de functies uit en geeft v terug
        oude_grenzen = verfijn_op_interval
        for i in range(aantal_verfijningen_berekening_v):
            print(i)
            oude_grenzen = grenzen(oude_grenzen)
        v = optimaliseet_het(oude_grenzen[0], oude_grenzen[1], 100)
        return v

    #Dit deel voert de functies uit en schrijft de gevonden gegevens naar de .txt bestanden
    write_file_snelheden = 'calculated data/%s/snelheden.txt' % sternaam
    write_file_dagen = 'calculated data/%s/dagen.txt' % sternaam
    if bereken == True:
        with open(write_file_snelheden, 'w') as f:
            f.write(str(bereken_v()) + "\n")
        with open(write_file_dagen, 'w') as g:
            g.write(str(header['BJD'])+"\t" + str(header['DATE-OBS'])+"\n")
    print("index: " + str(datalist.index(data_index)))