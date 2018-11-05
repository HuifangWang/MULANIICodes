dirname='/Users/huifangwang/NormalSEEG_CL/PS1001/mln/';
ind=3;
reference='arr';
data_filename='WI041222014';
its=0;
wins=20;
Pname='PS1001';
prenom=['ND',num2str(ind),'_',reference,'_',Pname,data_filename,'ts',num2str(its),'.mat'];

flagTF='T';
flagData='T';
mln_Cal_MULAN_fs_mlnvcs_wins(dirname,prenom,flagTF,flagData,wins)