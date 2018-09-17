#!/bin/bash

nskip=5

source cal_para.sh

WORKDIR="$PWD"
datadir=$WORKDIR/../
echo $WORKDIR
cd $WORKDIR
for beta in ${betaarray}; do
  for L  in ${Larray}; do
    #pretag=N${Nflavor}L${L}bX${bratio}
    pretag=beta${beta}dtau${dtau}L${L}
    rm $WORKDIR/${pretag}_m.dat
    rm $WORKDIR/${pretag}_binder.dat
    rm $WORKDIR/${pretag}_ekint.dat
    rm $WORKDIR/${pretag}_ecoup.dat
    rm $WORKDIR/${pretag}_ejs.dat
    rm $WORKDIR/${pretag}_ehx.dat
    rm $WORKDIR/${pretag}_etot.dat
    rm $WORKDIR/${pretag}_ising_Spipi.dat
    rm $WORKDIR/${pretag}_ising_Spipi_ccratio.dat

    for h in ${hxarray}; do
        cd $datadir
        maindir=b${beta}L${L}/h${h}
        if [ -f $maindir/ener1.bin ]; then
            cd $maindir
            echo "nskip = $nskip, processing $maindir ... "

            # analysis ener1.bin
            awk  '{if(NR>nskipv) print $0}' nskipv=$nskip ener1.bin > ener.tmp
            awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt((sumsq[i]/NR-sum[i]^2/NR^2)/(NR-1)) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==2) print hv, $0}' hv=$h ener.tmp2 >> $WORKDIR/${pretag}_m.dat
            awk '{if(NR==3) print hv, $0}' hv=$h ener.tmp2 >> $WORKDIR/${pretag}_binder.dat
            awk '{if(NR==6) print hv, $1/Lv/Lv, $2/Lv/Lv}' hv=$h Lv=$L ener.tmp2 >> $WORKDIR/${pretag}_ekint.dat
            awk '{if(NR==7) print hv, $1/Lv/Lv, $2/Lv/Lv}' hv=$h Lv=$L ener.tmp2 >> $WORKDIR/${pretag}_ecoup.dat
            awk '{if(NR==8) print hv, $1/Lv/Lv, $2/Lv/Lv}' hv=$h Lv=$L ener.tmp2 >> $WORKDIR/${pretag}_ejs.dat
            awk '{if(NR==9) print hv, $1/Lv/Lv, $2/Lv/Lv}' hv=$h Lv=$L ener.tmp2 >> $WORKDIR/${pretag}_ehx.dat
            rm ener.tmp2
            rm ener.tmp

            # analysis Ising spin correlation, S(pi,pi), from isingzztau_corrlt.bin
            # also calculate correlation ratio
            awk '{if(NR%(Lv*Lv)==1) print $4}' Lv=$L isingzztau_corrlt.bin |awk '{if(NR>nskipv) print $0}' nskipv=$nskip > Spipi.tmp
            awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt((sumsq[i]/NR-sum[i]^2/NR^2)/(NR-1)) )} }'\
            Spipi.tmp > Spipi.tmp2
            # output Spipi
            awk '{if(NR==1) print hv, $1, $2}' hv=$h Lv=$L Spipi.tmp2 >> $WORKDIR/${pretag}_ising_Spipi.dat

            awk '{if(NR%(Lv*Lv)==2) print $4}' Lv=$L isingzztau_corrlt.bin |awk '{if(NR>nskipv) print $0}' nskipv=$nskip > Spipidq.tmp
            awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt((sumsq[i]/NR-sum[i]^2/NR^2)/(NR-1)) )} }'\
            Spipidq.tmp > Spipidq.tmp2

            paste Spipi.tmp2 Spipidq.tmp2 |awk '{if(NR==1) print h, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1) }' h=$h >> $WORKDIR/${pretag}_ising_Spipi_ccratio.dat

            #rm Spipi.tmp2
            #rm Spipi.tmp
            #rm Spipidq.tmp2
            #rm Spipidq.tmp

        fi
    done
    cd $WORKDIR
    paste ${pretag}_ekint.dat ${pretag}_ecoup.dat ${pretag}_ejs.dat ${pretag}_ehx.dat |awk '{print $1, $2+$5+$8+$11, sqrt($3**2+$6**2+$9**2+$12**2)}' > ${pretag}_etot.dat
 done
done
