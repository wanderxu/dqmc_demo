#!/usr/bin/env python
"""
Author:          Xiao Yan Xu <wanderxu@gmail.com>

Description:

  get sigma(k,iwn)

"""
import math
import numpy as np
import scipy.integrate as integrate

# read L
fp = open( "in.para", 'r' )
lines = fp.readlines()
L = float(lines[0].split()[0])
ltrot = int(lines[0].split()[1])
dtau = 0.05
mu = -0.5

# define some temp array for storage kpoint and splope
k2 = []
####slope = []
####slerr = []

# wn
nuse=1
wn = []
for n in xrange(nuse):
    wn.append( (2.0*n)*np.pi*ltrot*dtau )


# read g(k,t)
tau = []
gk = []
gkerr = []
fp = open( "chitau.dat" )
lines = fp.readlines()
for line in lines[:]:
    tau.append(float(line.split()[0]))
    gk.append(float(line.split()[1]))
    gkerr.append(float(line.split()[2]))

# exp(iwnt)
##expiwnt=np.zeros((nuse,ltrot),dtype=complex)
coswnt=np.zeros((nuse,ltrot+1),dtype=float)
sinwnt=np.zeros((nuse,ltrot+1),dtype=float)
for nt in xrange(ltrot+1):
    for n in xrange(nuse):
        ##expiwnt[n,nt]=np.exp(1j*wn[n]*tau[nt])
        coswnt[n,nt]=np.cos(wn[n]*tau[nt])
        sinwnt[n,nt]=np.sin(wn[n]*tau[nt])

# g(k,t) to g(k,iwn)
gkwn = []
gkwnerr = []
for n in xrange(nuse):
    ztmp=complex(0.0,0.0)
    rretmp = integrate.simps(gk*coswnt[n], tau)
    rimtmp = integrate.simps(gk*sinwnt[n], tau)
    ztmp = complex(rretmp,rimtmp)
    gkwn.append(ztmp)

    rretmp = integrate.simps(gkerr*abs(coswnt[n]), tau)
    rimtmp = integrate.simps(gkerr*abs(sinwnt[n]), tau)
    gkwnerr.append(complex(rretmp,rimtmp))
print  gkwn[0].real, gkwnerr[0].real
