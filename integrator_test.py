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

upper_level_index = tar.atom_data.lines.set_index(['atomic_number', 'ion_number', 'level_number_upper']).index.copy()
e_dot_lu          = pd.DataFrame(tar.runner.Edotlu, index=upper_level_index)
e_dot_u           = e_dot_lu.groupby(level=[0, 1, 2]).sum()
e_dot_u.index.names = ['atomic_number', 'ion_number', 'source_level_number'] # To make the q_ul e_dot_u product work, could be cleaner
transitions       = tar.atom_data.macro_atom_data[tar.atom_data.macro_atom_data.transition_type == -1].copy()
transitions_index = transitions.set_index(['atomic_number', 'ion_number', 'source_level_number']).index.copy()
tmp  = tar.plasma.transition_probabilities[(tar.atom_data.macro_atom_data.transition_type == -1).values]
q_ul = tmp.set_index(transitions_index)
wave = tar.atom_data.lines.wavelength_cm[transitions.transition_line_id].values.reshape(-1,1)


if True:
    print "Integrator Luminosity {:.6e}".format(np.trapz(in_L_nu*un.erg,freq.cgs))
    print "MC         Luminosity {:.6e}".format(np.trapz(mc_L_nu.cgs,freq.cgs) )
    pl.plot(lam,in_L_lam,label="Source")
    pl.plot(lam,mc_L_lam,label="MC")
    pl.ylim(0,6e38)
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

    pl.ylim(0,6e38)
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
