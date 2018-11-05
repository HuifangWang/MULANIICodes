#!/bin/bash

resource="walltime=8:00:00"
property="host>='n02'"

#export MCR_CACHE_ROOT=$TMPDIR

for basicdir in ./mlnData/
do


wins=1500

dirname=$basicdir'Win'$wins/
#if [ ! -e $dirname ]; then
#mkdir $dirname;
#cp -r $basicdir'data' $dirname;
#fi

a=`ls $dirname'data'`

#degree='h'

for prename in $a
do


program="./run_mln_opt_para_vcs_FPN_G25_main.sh /soft/mcr/v85 $basicdir $prename $wins"
oarsub -p "$property" -l "$resource" "$program" 

done
done
