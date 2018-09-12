#!/bin/bash

source cal2.sh

for beta in $betaarray; do
    for Lsize in $Larray; do
        for hx in $hxarray; do
            WORKDIR="$PWD"
            datadir=$WORKDIR/dat/b${beta}L${Lsize}h${hx}/
            echo $WORKDIR
            cd $WORKDIR
            ii=0
            jj=0
            for bx in $shiftxarray; do
                for by in $shiftyarray; do
                    cd $WORKDIR
                    if [ -f gtauu${ii}${jj}.dat ]; then
                         rm gtauu${ii}${jj}.dat
                    fi
                    if [ -f gtaud${ii}${jj}.dat ]; then
                         rm gtaud${ii}${jj}.dat
                    fi
                    cd $datadir
                    jobdir=bx${bx}by${by}
                    if [ -d bx${bx}by${by} ]; then
                        cd $jobdir
                        awk -F ',' '{print $1}' gtau_up.bin | awk -F '(' '{if(NF==1) {print $1} else {print $2}}' > gtauu${ii}${jj}.dat
                        #awk '{if(NR%2==0) print$0}' gtau_up.bin > ct1
                        #awk '{if(NR%2==1) print$0}' gtau_up.bin > ct0
                        #paste ct0 ct1 > gtauu${ii}${jj}.dat
                        mv gtauu${ii}${jj}.dat $WORKDIR
                        #rm ct1 ct0
                        awk -F ',' '{print $1}' gtau_dn.bin | awk -F '(' '{if(NF==1) {print $1} else {print $2}}' > gtaud${ii}${jj}.dat
                        #awk '{if(NR%2==0) print$0}' gtau_dn.bin > ct1
                        #awk '{if(NR%2==1) print$0}' gtau_dn.bin > ct0
                        #paste ct0 ct1 > gtaud${ii}${jj}.dat
                        mv gtaud${ii}${jj}.dat $WORKDIR
                        #rm ct1 ct0
                    fi
                    ((jj=jj+1))
                done
                jj=0
                ((ii=ii+1))
            done
            cd $WORKDIR
            if [ -d b${beta}L${Lsize}h${hx} ]; then
               rm -r -f b${beta}L${Lsize}h${hx}
            fi
            mkdir b${beta}L${Lsize}h$hx
            mv gtau*.dat b${beta}L${Lsize}h$hx
        done
    done
done
