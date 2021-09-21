#===================================================================
#Computes fluxes and heating rates for the grey gas model.
#The fluxes are computed as a function of p/ps, given net optical
#thickness of the atmosphere tauinf .
#
#Since the OLR is just the upward flux at p=0, this can also
#be used to do grey gas OLR computations.  Different temperature
#profiles can be treated by just editing the
#functions T(pps) and dTdp(pps)
#
#The moist adiabat function in phys.py can be used as well,
#if one employs the option to create an interpolation function that
#returns temperature at an arbitrary pressure.
#===================================================================

#Data on section of text which this script is associated with
Chapter = '4.**'
Figure = 'fig:AllTropNetIRFluxGrey'
#
#This is also the solution script for
Problem = '{Workbook:RadBalance2:PressBroadenedHeating}'
#This script can also be modified to use for the problem
# '{Workbook:RadBalance2:StratTropOLRGrey}'

#**ToDo:  
#
#         *The optically thick approximation breaks down near the top
#          of the atmosphere, especially with pressure broadening included.
#          This makes plotting difficult. Currently it's handled with a value
#          cutoff, but there's probably a better way.

import math,phys
from ClimateUtilities import *

#Specify temperature as a function of p/ps
#pps is pressure divided by surface pressure. Note
#temperature gradient is written as dT/d(p/ps).
#Tstrat and Ts are set as globals
def T(pps):
    #This if block takes care of round off error, which
    #can make pressure go slightly negative at the top of
    #the atmosphere sometimes
    if pps > 0.:
        Tadiabat = Ts*(pps)**Rcp
    else:
        Tadiabat = 0.
    #Return the adiabatic temperature, or the stratospheric
    #temperature, if it is less
    return max(Tstrat,Tadiabat)
def dTdp(pps):
    #This if block takes care of round off error, which
    #can make pressure go slightly negative at the top of
    #the atmosphere sometimes
    if abs(T(pps) - Tstrat) < 1.e-6:
        dTadiabat = 0.
    else:
        dTadiabat = Rcp*Ts*(pps)**(Rcp-1)
    return dTadiabat


#Grey gas transmission function.
#tauinf is a global
def Trans(tau1,tau2):
    return math.exp(-abs(tau1-tau2))

#Integrand for upward or downward flux. Note that
#the Schwartzschild integral is written here as an integral
#over p/ps, and correspondingly the gradient of T is written as
#dT/d(p/ps). The solution is written in the form of
#Eq. (4.13) (in First Edition).
#
#Change log:
#     *5/20/2012: I changed the definition of tau1 and tau2
#      to correspond to the definition in the text. This doesn't
#      change the result, since Trans just depends on |tau1-tau2|
#
#     *5/20/2012: Fixed the boundary terms in Iplus and Iminus.
#      These terms didn't affect any results shown in the text but
#      make a difference if Tg differs from Ts, or (for Iminus)
#      if Tstrat is nonzero
#
def integrand(ppsp,params):
    #Without pressure broadening
    if PressureBroadening:
        tau1 = params.tauinf*(1.-ppsp**2)
        tau2 = params.tauinf*(1.-params.pps**2)
    else:
        tau1 = params.tauinf*(1.-ppsp)
        tau2 = params.tauinf*(1. - params.pps)    
    return Trans(tau1,tau2)*4.*phys.sigma*T(ppsp)**3*dTdp(ppsp)

def Iplus(pps,tauinf):
    params = Dummy()
    params.pps = pps
    params.tauinf = tauinf
    limit = min(1.,pps+10./tauinf)
    quad = romberg(integrand,10)
    if PressureBroadening:
        tau = tauInf*(1.-pps**2)
    else:
        tau = tauInf*(1.-pps)
    BddTerm = (phys.sigma*Tg**4 - phys.sigma*Ts**4)*Trans(0.,tau)
    return quad([pps,limit],params,.1)+ phys.sigma*T(pps)**4 +BddTerm

def Iminus(pps,tauinf):
    params = Dummy()
    params.pps = pps
    params.tauinf = tauinf
    limit = max(0.,pps-10./tauinf)
    quad = romberg(integrand,10)
    if PressureBroadening:
        tau = tauInf*(1.-pps**2)
    else:
        tau = tauInf*(1.-pps)
    return quad([pps,0.],params,.1)+ phys.sigma*T(pps)**4 - phys.sigma*Tstrat**4*Trans(tau,tauInf)

#This function returns a curve object containing
#the net upward flux as a function of p/ps (i.e. pressure
#relative to surface pressure) , and also gives
#the optically thick approximation to the net upward flux
#
#Note that the heating in the optically thick approximation becomes
#very large in the upper atmosphere, where, the approximation breaks
#down.  To keep this divergence from messing up the axes of the graph,
#in this routine, the heating rate is cut off at a maximum value, that
#is chosen to be large enough that one can see the divergence between
#the asymptotic and exact (numerical) result.  The flat regions
#of heating in the graph thus have no physical meaning. 
def heatList(tauinf):
    c = Curve()
    c.addCurve(p,'p')
    Ip = [Iplus(pps,tauinf) for pps in p]
    Im = [Iminus(pps,tauinf) for pps in p]
    h = [Ip[i]-Im[i] for i in range(len(p))]
    #**ToDo: Find some better way to keep the optically thick curve from messing up the plots
    if PressureBroadening:
        h1 = [(.5/(pps+1.e-10))*2.*4.*5.67e-8*T(pps)**3*dTdp(pps)/tauinf for pps in p]
    else:
        h1 = [2.*4.*5.67e-8*T(pps)**3*dTdp(pps)/tauinf for pps in p]
    #
    #Cut off the maximum of h1 so it doesn't wreck the plot of h
    maxh = max(h)
    h1 = [min(hh1,2.*maxh) for hh1 in h1]
    #
    c.addCurve(h,'NetFlux','Net Flux, Computed')
    c.addCurve(h1,'NetFluxThick','Net Flux, optically thick approx')
    #Set up options to plot as a profile
    c.switchXY = c.reverseY = True
    #Labels and title
    c.PlotTitle = 'tauInf = %f'%tauinf
    c.Ylabel = 'Flux, W/m**2'
    c.Xlabel = 'Normalized Pressure'
    return(c)

#These are all globals
p = [.01*i for i in range(101)]
Rcp = 2./7.
Tstrat = 0. #Stratospheric temperature
Ts = 300. #Surface air temperature
Tg = 300. #Ground temperature
#Say whether you want pressure broadening or not
PressureBroadening = False

#Do the plots
c1 = heatList(1.)
c10 = heatList(10.)
c50 = heatList(50.)

plot(c1)
#Note: In the text, we suppressed the plotting of
#the optically thick limit for the tauInf = 1 case, because
#the optically thick approximation is very inaccurate in this
#case and expands the scale of the graph so much it's hard to
#see the variation in the correct result. 
plot(c10)
plot(c50)

#This script can also be used to plot OLR vs. Tg or tauinf,
#as illustated below. The OLR is just Iplus(0,tauinf).
tauList = [.1*i for i in range(1,500)]
OLRList = [Iplus(0.,tauInf) for tauInf in tauList]
cOLR=Curve()
cOLR.PlotTitle ='OLR vs optical depth'
cOLR.addCurve(tauList,'TauInf')
cOLR.addCurve(OLRList,'OLR')
plot(cOLR)





