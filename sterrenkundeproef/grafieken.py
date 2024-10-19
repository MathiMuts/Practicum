import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

y = [40, 37, 28, 19, 14, 12, 13, 12, 10, 12, 11, 13, 12, 14, 19, 28, 37, 40]
x = [0, " ", "  ", "   ", "$v-\Delta v$", "    ", "     ", "               ", "      ", "       ", "        ", "         ", "          ", "$v+\Delta v$", "           ", "            ", "             ", "              "]
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
ax.plot(x, y, '.', markersize = 5, label = '', color='r')
ax.add_patch(Rectangle((4, 10), 9, 4.5,edgecolor = 'green',fill=False,lw=1))
ax.add_patch(Rectangle((4, 0), 0, 10,edgecolor = 'green',fill=False,lw=1,linestyle='dotted'))
ax.add_patch(Rectangle((13, 0), 0, 10,edgecolor = 'green',fill=False,lw=1,linestyle='dotted'))
ax.add_patch(Rectangle((13.25, 10), 0, 4.5,edgecolor = '#0A865E',fill=False,lw=1))
ax.text(13.5, 11.5, '$\Delta som\, I$', fontsize=10, color='#0A865E')
ax.set_xlim(0, 17)
ax.set_ylim(0, 50)
ax.set_xlabel('Snelheid $[km/s]$')
ax.set_ylabel('Intensiteit $[J/m^2]$')
ax.set_title('Lokale fout op v')
plt.show()

y = [55, 47, 38, 25, 14, 12, 13, 12, 10, 12, 11, 13, 12, 14, 25, 38, 47, 55]
x = [0, " ", "  ", "   ", "$\lambda '-\Delta \lambda '$", "    ", "     ", "               ", "      ", "       ", "        ", "         ", "          ", "$\lambda '+\Delta \lambda '$", "           ", "            ", "             ", "              "]
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=150, figsize=(5, 3))
ax.plot(x, y, '.', markersize = 5, label = '')
ax.add_patch(Rectangle((4, 10), 9, 4.5,edgecolor = 'red',fill=False,lw=1))
ax.add_patch(Rectangle((4, 0), 0, 10,edgecolor = 'red',fill=False,lw=1,linestyle='dotted'))
ax.add_patch(Rectangle((13, 0), 0, 10,edgecolor = 'red',fill=False,lw=1,linestyle='dotted'))
ax.add_patch(Rectangle((13.25, 10), 0, 4.5,edgecolor = 'orange',fill=False,lw=1))
ax.text(13.5, 11.5, "$\Delta I_{\lambda '}$", fontsize=10, color="orange")
ax.set_xlim(0, 17)
ax.set_ylim(0, 60)
ax.set_xlabel('Golflengte $[nm]$')
ax.set_ylabel('Intensiteit $[J/m^2]$')
ax.set_title('Lokale fout op I')
plt.show()