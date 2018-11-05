#!/bin/bash

resource="walltime=02:00:00"
property="host>='n13'"

basicdir='./mlnData'
parafile=$basicdir'opt_para_HCNS30.mat'

wins=1500
inwins=-1
dirname=$basicdir/Win$wins/


a=`ls $dirname'data'`

#degree='h'

for prename in $a
do


program="./run_mln_adp_seeg_pbg_main.sh /soft/mcr/v85 $dirname $prename $wins $inwins $parafile"
oarsub -p "$property" -l "$resource" "$program" 

done
