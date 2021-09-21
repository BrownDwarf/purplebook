#This script illustrates the runaway greenhouse phenomenon
#for a grey gas.  The atmosphere consists of a one-component
#saturated condensible atmosphere, with T(p) determined by the
#Clausius-Clapeyron relation.  The surface pressure (and hence
#optical thickness of the atmosphere) increases with T.  The
#script generates OLR(Tsurf) for this case. It assumes constant
#specific absorption cross-section kappa, but more general
#cases can be treated by editing the integrand function f
#and the expression for tauInf
#
#In this case, we find the OLR by carrying out the definite
#integral giving the solution to the Schwarzschild equations.
#The integral is carried out by integrating from the top of the
#atmosphere down to the ground, using delTau = tauInf-tau as the
#integration variable.
#
#
#
#**ToDo:  Turn the OLR computation into a function that
#         returns a Curve object. That will make it easier
#         for the students to modify the script to make graphs
#         exploring the parameter dependence
#
#**ToDo:  Replace the ODE integrator with call to romberg quadrature
#
#


#Data on section of text which this script is associated with
Chapter = '4.**'
Figure = '**'
#

from ClimateUtilities import *
import phys
from math import *


#Temperature profile. The argument is the ratio of
#pressure to the reference pressure
def T(pp0):
    pp0 = max(pp0,1.e-20) #Avoids math range error at top of atmosphere   
    return T0/(1. - RTL*log(pp0))

#Integrand function, for OLR computation
#Function has a single parameter, which is tau0.
#"dummy" is an unused argument, but has to be there for
#the ODE integrator
def f(delTau,dummy,tau0):
	delTau = min(delTau,100.)
	pp0 = delTau/tau0
	return phys.sigma*T(pp0)**4.*exp(-delTau)

#Function to compute OLR. This uses the differential
#equation integrator, but that could be replaced by
#quadrature, since we're just doing a definite integral
def OLR(ps):
    tauInf = kappa*ps/g
    tau0 = kappa*p0/g
    #Determine number of subdivisions
    n = int(10*tauInf) #Guarantee 10 subdivisions per unit optical depth
    n = max(n,50) #Don't let n get too small when we're optically thin
    ddelTau = tauInf/n
    #
    #Use the differential equation integrator to evaluate the
    #definite integral of f(delTau) over the atmosphere
    fi = integrator(f,0.,0.,ddelTau)
    #Pass the parameter of the function to the integrator
    #The integrator doesn't need to know tauInf, since that comes
    #in only through the limit of integration and the boundary
    #radiation term we add in at the end.
    fi.setParams(tau0)
    #
    #Integration loop
    #**ToDo: Could stop integrating when the contribution of the
    #deeper atmosphere becomes too small. That would save some time.
    for i in range(n-1):
	ans = fi.next()
    #
    #Now we have to add in the contribution from the surface
    tauInf = min(tauInf,100.) #Avoids underflow in optically thick case
    return ans[1] + phys.sigma*T(ps/p0)**4.*exp(-tauInf)


#--------Main script starts here-------------------------------
#Constants for Clausius-Clapeyron. Thes are globals
T0 = 300. #Reference temperature
p0 = phys.satvpg(T0) #Water vapor saturation pressure
RTL = phys.water.R*T0/phys.water.L_vaporization

#Absorption constants. Globals (used in OLR function)
g = 10.
kappa = .1 #Very roughly appropriate for water vapor

#Note: If you want to set p0 so that kappa p0/g = 1, as in the text,
#and then compute the corresonding T(p0), you can do the following instead
#The answer should be the same as using the fixed reference pressure above.
#p0 = g/kappa
#The next line estimates the temperature at the pressure p0, using
#                         the simplified form of Clausius-Clapeyron
#T0 = 300./(1.-(phys.water.R*300./phys.water.L_vaporization)*log(p0/phys.satvpg(300.)))
#tau0 = 1.
#RTL= phys.water.R*T0/phys.water.L_vaporization

#List of surface pressures we want to do the computation for
psList = [10.*(1.06)**i for i in range(-100,120)]
TsList = [T(ps/p0) for ps in psList] #Corresponding temperatures

#The following statement computes the OLR for each surface pressure ps.
#By Clausius-Clapeyron, ps also determines Ts.  It's more
#convenient to loop over ps than Ts, because we already have
#a function that computes T(p).
OLRList = [OLR(ps) for ps in psList]

#Plot results
c = Curve()
c.addCurve(TsList)
c.addCurve(OLRList)
w = plot(c)

