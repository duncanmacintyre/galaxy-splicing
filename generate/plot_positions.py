import numpy as np
import pandas
import re
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from whereswitch import whereswitch

# load the data
df = pandas.read_csv('galaxy_masses.csv')

# angular size distance to protocluster in kpc
# http://www.astro.ucla.edu/%7Ewright/CosmoCalc.html
# D_A = 1419.1*1e3
# D_A = 1421.3*1e3 # H0=67.8, Ωm = 0.309, ΩΛ = 0.691

# conversion factor from arcsec to kpc
ratio_kpc_arcsec = 7.

# k is a logical array of those galaxies not in SPIREc
k = df['label'].map(lambda s: re.match(r'SPIREc.', s) is None)

# red_RGB = np.array([1, 0, 0], ndmin=2)
# blue_RGB = np.array([0, 0, 1], ndmin=2)
# max_v = df['v'].abs().max()
# def cmap(v):
#     print(np.where(
#             v > 0,
#             v * red_RGB / np.abs(v).max(),
#             -v * blue_RGB / np.abs(v).max()
#         ))
#     return([0, 0, 0])

def plotposweights(df, circles=(90, 210), centre='C20', centre_radius=210, extra_cm_radius=750, colourbar=False, circle_text=True):

    if circles is None:
        circles = []

    # use SkyCoord class to parse the strings with ra and dec
    c = SkyCoord(ra=df['RA'], dec=df['DEC'],
                 unit=(u.hourangle, u.deg))

    # compute total masses of galaxies
    M = df['Mgas'] + df['M*'] + df['Mvir']

    # get galaxy coordinates in arcsec
    ra, dec = c.ra.arcsec, c.dec.arcsec

    # determine central point
    # this will be (0, 0) and the centre of the circles

    # use midpoint between C1 and C6
    # centre_ra = ra[(df['label'] == 'C1') | (df['label'] == 'C6')].mean()
    # centre_dec = dec[(df['label'] == 'C1') | (df['label'] == 'C6')].mean()

    # find ra, dec coordinates of galaxy to use as temporary centre
    centre_ra = ra[(df['label'] == centre)]
    centre_dec = dec[(df['label'] == centre)]

    # compute coordinates in kiloparsecs
    x = ratio_kpc_arcsec * (ra - centre_ra) * np.cos(c.dec.radian) # position north-south
    # the cosine accounts for the lines of declination getting shorter the larger the declination
    y = ratio_kpc_arcsec * (dec - centre_dec) # position east-west

    # determine which galaxies are within centre_radius of the galaxy called centre
    l = np.sqrt(x**2 + y**2) < centre_radius

    # compute centre of mass of region within this radius using weighted average
    cmass_x, cmass_y = (x[l] * M[l]).sum() / M[l].sum(), (y[l] * M[l]).sum() / M[l].sum()
    # change coordrinates so that centre of mass is at (0, 0)
    x, y = x - cmass_x, y - cmass_y

    fig, ax = plt.subplots(constrained_layout=True)

    # plot the circles
    angle = np.linspace(0, 2*np.pi, num=1000)
    for r in circles:
        ax.plot(r*np.cos(angle), r*np.sin(angle), color='grey', linestyle=':')
        if circle_text:
            ax.annotate(str(r) + ' kpc', np.array((1,1))*np.sqrt(0.5)*(r+5), c='grey', rotation=45,
                        horizontalalignment='center', verticalalignment='center')

    # plot an X at the centre of mass of the central area
    ax.scatter(0, 0, c='brown', marker='x')

    # plot another X if extra_cm_radius given
    if extra_cm_radius is not None:
        l2 = np.sqrt(x**2 + y**2) < extra_cm_radius
        cmass2_x, cmass2_y = (x[l2] * M[l2]).sum() / M[l2].sum(), (y[l2] * M[l2]).sum() / M[l2].sum()
        ax.scatter(cmass2_x, cmass2_y, c='darkgrey', marker='x')

    # plot the galaxies
    # the marker area is proportional to galaxy mass
    s = ax.scatter(x, y,
               c=df['v'],
               s=M,
               marker='.')

    # add colour bar to show velocities, if set
    if colourbar:
        fig.colorbar(s, ax=ax)

    # plot labels for galaxies
    for ix, iy, itxt in zip(x, y, df['label']):
        ax.annotate(str(itxt), (ix, iy), xytext=(ix+0.2, iy+0.2)) 

    ax.set_aspect('equal')

    # x-axis should go from positive to negative (to match Hill+ 2020)
    lim = ax.get_xlim()
    ax.set_xlim(lim[1], lim[0])

    # distance from centre of mass
    # r = np.sqrt(x**2 + y**2)
    # fig2, [ax2a, ax2b] = plt.subplots(2, constrained_layout=True)
    # ax2a.scatter(r[df['Mgas_method'] == 'CO43'],
    #              M[df['Mgas_method'] == 'CO43'],
    #              label='CO43',
    #              marker = 'o')
    # ax2a.scatter(r[df['Mgas_method'] == 'Cii'],
    #              M[df['Mgas_method'] == 'Cii'],
    #              label='Cii',
    #              marker = '*')
    # ax2a.scatter(r[df['Mgas_method'] == 'SFR'],
    #              M[df['Mgas_method'] == 'SFR'],
    #              label='SFR',
    #              marker = '+')
    # ax2b.scatter(r[df['Mgas_method'] == 'CO43'],
    #              df['v'][df['Mgas_method'] == 'CO43'],
    #              label='CO43',
    #              marker = 'o')
    # ax2b.scatter(r[df['Mgas_method'] == 'Cii'],
    #              df['v'][df['Mgas_method'] == 'Cii'],
    #              label='Cii',
    #              marker = '*')
    # ax2b.scatter(r[df['Mgas_method'] == 'SFR'],
    #              df['v'][df['Mgas_method'] == 'SFR'],
    #              label='SFR',
    #              marker = '+')
    # ax2a.set_xlabel('Distance from centre [kpc]')
    # ax2b.set_xlabel('Distance from centre [kpc]')
    # ax2a.set_ylabel('Mass [1e10 M☉]')
    # ax2b.set_ylabel('Radial velocity [km/s]')
    # ax2a.legend()
    # ax2b.legend()



plotposweights(df, circles=(90, 95), centre_radius=95)

plt.show()
