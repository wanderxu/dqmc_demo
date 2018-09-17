#!/bin/bash

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
EXE=../../../src/ftdqmc
cd $WORKDIR
for beta in ${betaarray}; do
  for L  in ${Larray}; do
    for h in ${hxarray}; do
        cd $WORKDIR

        maindir=b${beta}L${L}
        if [ ! -d $maindir ]; then
            mkdir $maindir
        fi
        cd $maindir

        jobdir=h${h}
        if [ ! -d $jobdir ]; then
            mkdir $jobdir
        fi
        cd $jobdir

        if [ -f confout ]; then
            cp confout confin
        fi
		#cp $WORKDIR/Heff.para .

cat>ftdqmc.in<<endin
&model_para
L = $L,
beta = $beta,
dtau = $dtau,
mu =  $mu,
muA = $muA,
muB = $muB,
rhub = $rhub,
rj = $rj,
js = $js,
hx = $h,
xmag = $xmag,
flux_x = $flux_x,
flux_y = $flux_y,
/
&ctrl_para
nwrap = $nwrap,
nsweep = $nsweep,
nbin = $nbin,
llocal = $llocal,
nsw_stglobal = $nsw_stglobal,
lsstau = $lsstau,
ltau = $ltau,
nuse = $nuse,
nublock = $nublock 
/
endin

cat>run.sub<<endsub
#!/bin/bash
#PBS -N $beta.$L.$h
#PBS -o jobrun.out
#PBS -e jobrun.err
#PBS -l nodes=1:ppn=24
#PBS -l walltime=240:00:00

ncpu=\`cat \$PBS_NODEFILE | wc -l\`

cd \$PBS_O_WORKDIR

mpirun -machinefile \$PBS_NODEFILE -np \$ncpu $EXE

endsub
        qsub run.sub
    done
  done
done
