function mln_Save_LFP4MULAN()

dirname='/Users/huifangwang/MULAN3/mlnData/LFP/Slowtheta/20130516_CTR/raw/';
prenom='20130516_CTR';
Nds=1;
%Nds=[2,5,10];
load([dirname,prenom,'.mat'],'eegHpc','eegEC','Par')
oldData=[eegHpc';eegEC'];
ts=0;
periodt=12*60;
fs=Par.lfpSampleRate;
chanchosen=[9:16,17:24,34:2:64];
cdata=oldData(chanchosen,ts*fs+1:(ts+periodt)*fs);
str={};
neegH=16;
for ieegH=1:neegH
str{ieegH}=['H',num2str(chanchosen(ieegH))];
end
for ieegEC=17:32
str{ieegEC}=['EC',num2str(chanchosen(ieegEC)-32)];
end
Params.str=str;

output_directory=[dirname,'data/'];
if ~exist(output_directory, 'dir')
  mkdir(output_directory);
end
for ids=1:length(Nds)
Data=downsample(cdata',Nds(ids))';
Params.fs=fs/Nds(ids);

save([output_directory,prenom,'ds',num2str(Nds(ids)),'.mat'],'Data','Params')
end