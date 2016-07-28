import matplotlib.pyplot as pl
import astropy.constants as co
import astropy.units as un
import pandas as pd
import numpy as np

import tardis

tar = tardis.run_tardis("simple_example.yml")
#tar = tardis.run_tardis("tardis_example.yml")

freq = tar.runner.spectrum.frequency
lam  = tar.runner.spectrum_virtual.wavelength
mc_L_nu = tar.runner.spectrum_virtual.luminosity_density_nu
mc_L_lam = tar.runner.spectrum_virtual.luminosity_density_lambda

in_L_nu = tar.runner.L_nu
in_L_lam = co.c.to("Angstrom/s")/lam**2 * in_L_nu

if True:
    print "Integrator Luminosity {:.6e}".format(np.trapz(in_L_nu*un.erg,freq.cgs))
    print "MC         Luminosity {:.6e}".format(np.trapz(mc_L_nu.cgs,freq.cgs) )
    pl.plot(lam,in_L_lam,label="Source")
    pl.plot(lam,mc_L_lam,label="MC")
    pl.xlabel(u"A")
    pl.ylabel("erg/A/s")
    pl.title("Comparison")
    pl.ylim(0,4e38)
    pl.legend(loc="best")
    pl.show()

def plot_srcfun(mdl):
    sel = mdl.runner.wave[:,0].argsort()
    [pl.axvline(x=X,ymin = .85,ymax=.95,color='r') for X in tar.atom_data.lines.wavelength]
    pl.stem(1e8*mdl.runner.wave[sel],mdl.runner.att_S_ul)
    pl.title("Source function")
    pl.gca().set_yscale('log')
    pl.gca().set_xscale('log')
    pl.ylim(0,1e-3)    
    pl.xlabel("A")
    pl.show()


def complot():
    p10lam = np.load("10plam.npy")
    p20lam = np.load("20plam.npy")
    p50lam = np.load("50plam.npy")
    p100lam = np.load("100plam.npy")
    p200lam = np.load("200plam.npy")

    lam     = np.load("lam.npy")
    mc_L_lam = np.load("MClam.npy")

    pl.plot(lam,p10lam,label="10 p")
    pl.plot(lam,p20lam,label="20 p")
    pl.plot(lam,p50lam,label="50 p")
    pl.plot(lam,p100lam,label="100 p")
    pl.plot(lam,p100lam,label="200 p")
   
    pl.plot(lam,mc_L_lam,linewidth=1.2,label="MC")
 
    pl.xlabel(u"A")
    pl.ylabel("erg/A/s")
    pl.title("Comparison")

    pl.ylim(0,3.7e38)
    pl.legend(loc="best")
    pl.show()


def summary():
    lam      = np.load("lam.npy")
    mc_L_lam = np.load("MClam.npy")  

    p10lam   = np.load("10plam.npy")
    p20lam   = np.load("20plam.npy")
    p50lam   = np.load("50plam.npy")
    p100lam  = np.load("100plam.npy")
    p200lam  = np.load("200plam.npy")

    srcs   = [p10lam, p20lam, p50lam, p100lam, p200lam] 
    labels = ["10 p","20 p","50 p","100 p","200 p"]

    for i,src in enumerate(srcs):
        pl.plot(lam,src,label=labels[i])
        pl.plot(lam,mc_L_lam,linewidth=1.2,label="MC")
        pl.ylim(0,6e38)
        pl.xlabel(u"A")
        pl.ylabel("erg/A/s")
        pl.title(labels[i])
        pl.legend(loc="best")
        pl.show()


def lam_vary_summary():
    lam05  = np.load("lam.npy")
    lam15  = np.load("1500lam.npy")
    lam55  = np.load("5500lam.npy")
    lam100 = np.load("10000lam.npy")
    int05  = np.load("20plam.npy")
    int15  = np.load("1500lam20p.npy")
    int55  = np.load("5500lam20p.npy")
    int100 = np.load("10000lam20p.npy")

    lams = [lam05, lam15, lam55, lam100]
    ints = [int05, int15, int55, int100]
    labels = ["500 lambda", "1.5k lambda","5.5k lambda","100k lambda"]

    for i,(lam,src) in enumerate(zip(lams,ints)):
        pl.plot(lam,src,label=labels[i])
    pl.ylim(0,4e38)
    pl.xlabel(u"A")
    pl.ylabel("erg/A/s")
    pl.title("Varying number of wavelengths")
    pl.legend(loc="best")
    pl.show()
