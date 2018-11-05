#!/bin/bash
resource="walltime=00:30:00"
property="host>='n02'"



flagTF=T
flagData=M

for basicdir in mlnData/
do

wins=1500

dirname=$basicdir'Win'$wins/
if [ ! -e $dirname ]; then
mkdir $dirname;
cp -r $basicdir'data' $dirname;
fi

a=`ls $dirname'data'`


for prename in $a
do

program="./run_mln_Cal_MULAN_fs_mlnvcs_wins_BM.sh /soft/mcr/v85 $dirname $prename $wins"
oarsub -p "$property" -l "$resource" "$program" 

done
done
