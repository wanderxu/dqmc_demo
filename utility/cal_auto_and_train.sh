#!/bin/bash
# NOTE: please set codedir first
codedir=
autoexe=$codedir/utility/auto/auto.py
x_confc=$codedir/utility/ana_confc/x_confc
#trainexe=$codedir/utility/trainning/chi-square.py
trainexe=$codedir/utility/trainning/chi-square.rsquare.py

# calculate auto-correlation time
####echo " calculating auto-correlation time ... "
#####awk '{print $3}' ener1.bin >in.dat
####awk '{if(NR>500) print $0}' totsz.bin > in.dat
####python $autoexe

    nnimax_array=$( echo "1" )
#   nnimax_array=$(awk 'BEGIN{for(i=1;i<11;i+=1) printf("%4i",i)}')
#   nntmax_array=$( echo "10" )
    nntmax_array=$(awk 'BEGIN{for(i=1;i<11;i+=1) printf("%4i",i)}')
 nnimax_hyb_array=$( echo "0" )
 nntmax_hyb_array=$( echo "0" )
#nntmax_hyb_array=$(awk 'BEGIN{for(i=1;i<11;i+=1) printf("%4i",i)}')

for nnimax in $nnimax_array; do
for nntmax in $nntmax_array; do
for nnimax_hyb in $nnimax_hyb_array; do
for nntmax_hyb in $nntmax_hyb_array; do

    if [ "$nnimax_hyb" -le "$nnimax" ]; then
        echo " "
        echo " ### processing nnimax = $nnimax nntmax = $nntmax nnimax_hyb = $nnimax_hyb nntmax_hyb = $nntmax_hyb "
        echo " "
        echo " analysis confout.bin, generate data for tranning ... "
        L=$( grep -w L ftdqmc.out|awk '{print $3}' )
        ltrot=$( grep  "^ ltrot" ftdqmc.out|awk '{print $3}' )
        echo "$L $ltrot $nnimax $nntmax $nnimax_hyb $nntmax_hyb" > in.para
        $x_confc > train.dat.tmp
        awk '{if(NR>500) print $0}' train.dat.tmp >train.dat

        
        echo " tranning ... "
        python $trainexe > train.log.$nnimax.$nntmax.$nnimax_hyb.$nntmax_hyb

    fi
done
done
done
done
