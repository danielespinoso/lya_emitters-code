import pylab
import numpy
import cosmolopy.luminosityfunction as lf
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import scipy.integrate as si
import math

#z=2.2:
#limflux=9.81201246187e+43
#limflux=2.04482425923e+43
#limflux=3.25905921321e+43 #F378
limflux=3.02391290526e+43 #F395
limvol=0.14*10**(-9)

lumar=numpy.arange(41,47.5,0.05)

lstar=42.59
fi=-3.09
alpha=-1.75
a=59.4
b=-1.48

greens_x=numpy.array([42.30,42.55,42.8,43.05,43.3,43.55])
greens_y=10**numpy.array([-2.76,-3.11,-3.6,-4.4,-4.62,-5.15])

def sch_lum(lum,fi,alpha,lstar):
    return fi*(lum/lstar)**alpha*numpy.e**(-lum/lstar)

def sch_lum_exp(lum,fi,alpha,lstar):
    return fi*numpy.log(10)*10**((alpha+1)*(lum-lstar))*numpy.exp(-10**(lum-lstar))

def lin_lum(lum,a,b):
    return 10**(a+b*numpy.log10(lum))

def lin_lum_exp(lum,a,b):
    return 10**(a+b*lum)

def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return b

x1=numpy.log10(limflux)
x2=numpy.interp(limvol,sch_lum_exp(lumar,10**fi,alpha,lstar)[::-1],lumar[::-1])
x3=numpy.interp(limvol,lin_lum_exp(lumar,a,b)[::-1],lumar[::-1])

x_sc=numpy.arange(x1,x2,0.01)
y_sc=sch_lum_exp(x_sc,10**fi,alpha,lstar)

x_lin=numpy.arange(x1,x3,0.01)
y_lin=lin_lum_exp(x_lin,a,b)

tot_sc=numpy.trapz(y_sc,x=x_sc)
tot_lin=numpy.trapz(y_lin,x=x_lin)

vol=8.0958517311*10**9

print "J-PLUS galaxies, Schechter assumption:", tot_sc*vol
print "J-PLUS galaxies, linear assumption:", tot_lin*vol



import numpy as np
from scipy.interpolate import interp1d
from tools import *
filt_set = ['J0395', 'gJAVA', 'rJAVA']
zmin, zmax = zl(band=filt_set[0], wlength=1215.67)
method = '2FM'

n_tiles = 53
jplus_area = (1.4**2.)*n_tiles #now in sq.degrees
sphere_area = 4*np.pi/((np.pi/180.)**2.) #now in sq.degrees
fsky = jplus_area/sphere_area
dLogL = 0.12

name = '/home/CEFCA/dspinoso/works/lya_emitters/tools/comoving_distance.txt'
zeta, eta = np.genfromtxt(name, skip_header=1, unpack=True)
line = interp1d(zeta, eta, kind='nearest', bounds_error=None, fill_value='extrapolate')
dc_min = line(zmin);  dc_max = line(zmax) #now in Mpc
dp_min = line(zmin)/(1+zmin);  dp_max = line(zmax)/(1+zmax) #now in Mpc
dVc = (4.*np.pi*(dc_max**3. - dc_min**3.)/3.)*fsky
#dVp = abs((4.*np.pi*(dp_max**3. - dp_min**3.)/3.))*fsky
dV = dVc


final_list = '/home/CEFCA/dspinoso/works/lya_emitters/results/candidates_allparams_'+method+'.txt'
lum = np.genfromtxt(final_list, skip_header=1, usecols=(6), unpack=True)
Loglum = np.log10(lum)
notinf = (Loglum > 0)
Loglum = Loglum[notinf]
Lmin = min(Loglum) ; Lmax = max(Loglum)
centers, histo = histogram(Loglum, Lmin, Lmax, dLogL)
phi = histo*(1./(dV*dLogL))#*0.1
errs = 1/np.sqrt(histo)*(1./(dV*dLogL))

old_final_list = '/home/CEFCA/dspinoso/works/lya_emitters/results/old_candidates_allparams_'+method+'_zSelected.txt'
old_lum = np.genfromtxt(old_final_list, skip_header=1, usecols=(6), unpack=True)
old_Loglum = np.log10(old_lum)
old_notinf = (old_Loglum > 0)
old_Loglum = old_Loglum[old_notinf]
old_Lmin = min(old_Loglum) ; old_Lmax = max(old_Loglum)
old_centers, old_histo = histogram(old_Loglum, old_Lmin, old_Lmax, dLogL)
old_phi = old_histo*(1./(dV*dLogL))
old_errs = 1/np.sqrt(old_histo)*(1./(dV*dLogL))


violet = (0.6, 0., 0.6)
fig, ax = plt.subplots()

plt.plot(greens_x,greens_y,'go',markersize=10)
plt.plot(lumar,sch_lum_exp(lumar,10**fi,alpha,lstar),'b-',linewidth=2)
plt.plot(lumar,lin_lum_exp(lumar,a,b),'b--',linewidth=2)
#plt.plot(old_centers, old_phi, 'or', markersize=7) ### old_daniele
plt.plot(centers, phi, 'o', c=violet, markersize=7) ### daniele
plt.axhline(y=1.37*10**(-6),ls='dotted',color='k',lw=2)
plt.axhline(y=0.12*10**(-9),ls='--',color='k',lw=2)
plt.axvline(x=numpy.log10(limflux),color='k',ls='--',lw=2)
plt.axvline(x=43.92,color='r',ls='--',lw=2)
plt.fill_between(x_sc,y_sc,y2=0.14*10**(-9),facecolor='black',alpha=0.8,edgecolor='None')
plt.fill_between(x_lin,y_lin,y2=0.14*10**(-9),facecolor='grey',alpha=0.8,edgecolor='None')
ax.set_yscale('log')
plt.ylim([10**(-10),10**(-2)])
plt.xlim([42,47])
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
#ax.set_xscale('log')
#ax.xaxis.set_major_formatter(tick.FuncFormatter(fmt))
ax.yaxis.set_major_formatter(tick.FuncFormatter(fmt))
plt.xlabel("$\mathrm{log} L_{\mathrm{Ly}\\alpha} [\mathrm{erg s}^{-1}]$",size=20)
plt.ylabel("$\mathrm{log}\Phi[\mathrm{Mpc}^{-3}(\mathrm{dlog}L)^{-1}]$",size=20) 
plt.savefig("be_jplus_z2.pdf",bbox_inches='tight')
plt.show()
