#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
datadir=$WORKDIR/dat/
echo $WORKDIR
cd $WORKDIR
for h in $hxarray; do
  for L in $Larray; do
    cd $WORKDIR
    if [ -f chi00_h${h}L${L}.dat ]; then
        rm  chi00_h${h}L${L}.dat
    fi
    for beta in $betaarray; do
        cd $datadir
        jobdir=h${h}/L${L}b${beta}
        if [ -d h${h}/L${L}b${beta} ]; then
            cd $jobdir
            awk '{print $0}' trim.bin \
            |awk '{ if(NR>0) print $0 }' > chi00.tmp
            awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.8f %12.8f \t", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR) )} }'\
            chi00.tmp | awk '{if(NR==1) print rhx, $0 }' rhx=$beta >> $WORKDIR/chi00_h${h}L${L}.dat
            rm chi00.tmp
            #mv isingzztau_corrlt.bin isingzztau_corrlt.bin.tmp
            #awk '{ if((NR-1)%21 > 9 && (NR-1)%21 < 21) print $0}' isingzztau_corrlt.bin.tmp > isingzztau_corrlt.bin
        fi
    done
  done
done
