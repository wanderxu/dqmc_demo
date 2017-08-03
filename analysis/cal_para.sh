#!/bin/bash
#Larray=$(echo '9 12 15 21 24') !The lattice size must be multiple of 3
Larray=$(echo '15')
#betaarray=$(echo '0.5 0.6 0.8 1 1.2 1.5 2 2.5 3 4 5 6 8 10')
betaarray=$(echo '30')
nwrap=15
rhub=1.0
#hxarray=$(awk 'BEGIN{for(i=2.00;i<=5.15;i+=0.10) printf("%6.2f",i)}')
#hxarray=$(echo '1.00 1.50 2.00 2.50 2.60 2.70 2.80 2.90 3.20 3.40 3.60 3.80 4.00 4.50 5.00 6.00')
#hxarray=$(echo '3.16 3.18 3.22 3.24 3.26 3.28 3.30 3.32 3.35 3.45 3.50 3.55')
#hxarray=$(echo '3.00')
hxarray=$(echo '1.72 1.73 1.74 1.75 1.76 1.77 1.78 1.79 1.80 1.81 1.82 1.83')
#hxarray=$(echo '1.795')
mu=-0.5
nsw_stglobal=1
nsweep=100
nbin=30
xmag=0.0
shiftx=0.0
shifty=0.0
lsstau=F
lsstau0r=F
ltau=F
ltauall=F
nuse=0
js=-1
dtau=0.05

echo " L = " $Larray
echo " beta = " $betaarray
echo " nwrap = " $nwrap
echo " hx = " $hxarray
declare -i num_hx
num_hx=0
for hxtmp in $hxarray; do
    let num_hx=num_hx+1
done
echo " num_hx = " $num_hx
