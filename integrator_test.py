import scipy.interpolate as intr
import matplotlib.pyplot as pl
import astropy.constants as co
import astropy.units as un
import pandas as pd
import numpy as np

import tardis

def ratiocompare():
    pl.plot(lam,in_L_lam/mc_L_lam)
    pl.xlabel(u"A")
    pl.title("Ratio")
    pl.show()

tar = tardis.run_tardis("simple_example.yml")

freq = tar.runner.spectrum.frequency
lam  = tar.runner.spectrum_virtual.wavelength
mc_L_nu = tar.runner.spectrum_virtual.luminosity_density_nu
mc_L_lam = tar.runner.spectrum_virtual.luminosity_density_lambda

in_L_nu = tar.runner.L_nu
in_L_lam = co.c.to("Angstrom/s")/lam**2 * in_L_nu
BB = tardis.util.intensity_black_body(tar.runner.spectrum_virtual.frequency,6416)
Fnu = BB*np.pi
Lbb_nu = 4*np.pi*tar.runner.r_inner_cgs.min()**2 * Fnu
Lbb_lam = Lbb_nu * co.c.to("Angstrom/s")/lam**2 


if True:
    print "Integrator Luminosity {:.6e}".format(np.trapz(in_L_nu*un.erg,freq.cgs))
    print "MC         Luminosity {:.6e}".format(np.trapz(mc_L_nu.cgs,freq.cgs) )
    print "BB         Luminosity {:.6e}".format(np.trapz(Lbb_nu*un.erg,freq.cgs))
    ratiocompare()

def saveit(name):
    saved = {"lam":lam.value,"mc_L_lam":mc_L_lam.value,"in_L_lam":in_L_lam.value}
    np.savez(name,**saved)

def saverateio(name):
    an = np.loadtxt("analytic_solution.dat")
    np.save(name,in_L_lam.value/an[:,1]/Lbb_lam.value)

def pure_absorb_plot():
    lam_start,lam_end = tar.tardis_config.spectrum.start,tar.tardis_config.spectrum.end
    bins       = tar.runner.spectrum.wavelength[::-1]
    noint_mask = tar.runner.virt_packet_last_interaction_type == -1
    noint_lam  = (co.c.cgs / (tar.runner.virt_packet_nus[noint_mask] * un.Hz)).to("AA")
    noint_ws   = tar.runner.virt_packet_energies[noint_mask] * un.erg / tar.time_of_simulation
    Lnorm      = noint_ws[(noint_lam >= lam_start)*(noint_lam<=lam_end)].sum()

    pl.plot(lam,in_L_lam,label="Source",color="red")
    ax  = pl.gca()
    ret = ax.hist(noint_lam, weights=noint_ws, bins=bins, normed=True)
    for reti in ret[-1]:
        reti.set_facecolor("blue")
        reti.set_edgecolor("blue")
        reti.set_linewidth(0)
        reti.set_height(reti.get_height() * Lnorm.value)
    pl.plot(lam,in_L_lam,label="Source",color="red")
    pl.show()

def compare():
    pl.plot(lam,in_L_lam,label="Source")
    pl.step(lam,mc_L_lam,label="MC")
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

def lam_p_vary_summary():
    fnames = ["10p15clam","20p15clam","50p15clam","100p15clam","200p15clam","20pr55clam","200pr55clam",]

    for fil in fnames:
        dat = np.load(fil+".npz")
        pl.semilogy(dat["lam"],dat["in_L_lam"],label=fil)
    pl.semilogy(dat["lam"],dat["mc_L_lam"],label="MC")
    pl.xlabel("A")
    pl.ylabel("erg/A/s")
    pl.title("Varying number of wavelengths and ps")
    pl.legend(loc="best")
    pl.show()

def lam_p_vary_20():
    dat15 = np.load("20pr55clam.npz")
    dat55 = np.load("200pr55clam.npz")
    pl.plot(dat15["lam"],dat15["in_L_lam"],label="20p 5500 lam")
    pl.plot(dat55["lam"],dat55["in_L_lam"],label="200p 5500 lam")
    pl.xlabel("A")
    pl.ylabel("erg/A/s")
    pl.title("Varying number of wavelengths and ps")
    pl.legend(loc="best")
    pl.show()

def BBcomp():
    pl.plot(lam,in_L_lam,label="Source")
    pl.plot(lam,mc_L_lam,label="MC")
    pl.plot(lam,Lbb_lam,label="BB")
    pl.ylabel("erg/A/s")
    pl.xlabel(u"A")
    pl.title("Comparison")
    pl.legend(loc="best")
    pl.show()

def BBnorm():
    an = np.loadtxt("analytic_solution.dat")
    pl.plot(lam,in_L_lam/Lbb_lam,label="Source")
    pl.plot(lam,mc_L_lam/Lbb_lam,label="MC")
    pl.plot(an[:,0],an[:,1],label="Exact")
    pl.ylabel("erg/A/s")
    pl.xlabel(u"A")
    pl.title("Comparison")
    pl.legend(loc="best")
    pl.ylim(0,2)
    pl.show()

def Hsellnorm():
    an = np.loadtxt("analytic_solution.dat")
    pl.plot(lam,in_L_lam/an[:,1]/Lbb_lam)
    pl.show()

def geotest():
    r1 = np.load("4sh40p.npy")
    r2 = np.load("10sh40p.npy")
    r3 = np.load("40sh40p.npy")
    r4 = np.load("70sh40p.npy")
    r5 = np.load("40sh100p.npy")
    pl.plot(lam,r1,label="4 sh 40 p")
    pl.plot(lam,r2,label="10 sh 40 p")
    pl.plot(lam,r3,label="40 sh 40 p")
    pl.plot(lam,r4,label="70 sh 40 p")
    pl.plot(lam,r5,label="40 sh 100 p")
    pl.xlabel(u"A")
    pl.title("Comparison")
    pl.legend(loc="best")
    pl.show()
