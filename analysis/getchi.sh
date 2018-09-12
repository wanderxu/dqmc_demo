#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
datadir=$WORKDIR/../xi1muo5/
echo $WORKDIR
cd $WORKDIR
for beta in $betaarray; do
  for L in $Larray; do
    cd $WORKDIR
    if [ -f chi00_b${beta}L${L}.dat ]; then
        rm  chi00_b${beta}L${L}.dat
    fi
    for h in $hxarray; do
        cd $datadir
        jobdir=b${beta}L${L}/h${h}
        if [ -d b${beta}L${L}/h${h} ]; then
            cd $jobdir
            awk '{if($1==0.0 && $2==0.0 && $3==0.0) print $4}' isingzztau_corrlt.bin \
            |sort -n |awk '{ if(NR>0) print $0 }'|tac | awk '{ if(NR>0) print $0 }' > chi00.tmp
            awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR) )} }'\
            chi00.tmp | awk '{if(NR==1) print rhx, $0 }' rhx=$h >> $WORKDIR/chi00_b${beta}L${L}.dat
            rm chi00.tmp
            #mv isingzztau_corrlt.bin isingzztau_corrlt.bin.tmp
            #awk '{ if((NR-1)%21 > 9 && (NR-1)%21 < 21) print $0}' isingzztau_corrlt.bin.tmp > isingzztau_corrlt.bin
        fi
    done
  done
done
