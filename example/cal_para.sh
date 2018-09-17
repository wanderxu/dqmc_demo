#!/bin/bash

# model_para
Larray=$(echo '4')
betaarray=$(echo '0.50')
dtau=0.05
mu=-1.11856
muA=0
muB=0
rhub=1
rj=0
js=-1
hxarray=$(awk 'BEGIN{for(i=0.50;i<=3.001;i+=0.50) printf("%6.2f",i)}')
xmag=1
flux_x=0
flux_y=0

#ctrl_para

nwrap=10
nsweep=100
nbin=20
llocal=T
nsw_stglobal=1
lsstau=T
ltau=F
nuse=0
nublock=16
