#!/bin/bash
#Larray=$(echo '9 12 15 21 24') !The lattice size must be multiple of 3
Larray=$(echo '36')
#betaarray=$(echo '0.5 0.6 0.8 1 1.2 1.5 2 2.5 3 4 5 6 8 10')
#betaarray=$(echo '2 4 6 8 10 14 18 22 26 30 33 36 40')
betaarray=$(echo '20')
nwrap=30
rhub=0
#hxarray=$(awk 'BEGIN{for(i=2.00;i<=5.15;i+=0.10) printf("%6.2f",i)}')
#hxarray=$(echo '1.00 1.50 2.00 2.50 2.60 2.70 2.80 2.90 3.20 3.40 3.60 3.80 4.00 4.50 5.00 6.00')
#hxarray=$(echo '3.16 3.18 3.22 3.24 3.26 3.28 3.30 3.32 3.35 3.45 3.50 3.55')
#hxarray=$(echo '3.00')
#hxarray=$(echo '2.80 2.90 3.00 3.10 3.20 3.27 3.30 3.40 3.50 3.60 3.70 3.80 3.90 4.00')
hxarray=$(echo '0.30')
mu=-0.5
nsw_stglobal=1
nsweep=250
nbin=3
xmag=0.0
shiftxarray=$(echo '0.0001 0.25 0.5')
shiftyarray=$(echo '0.0002 0.25 0.5')
lsstau=F
lsstau0r=F
ltau=T
ltauall=T
nuse=12
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
