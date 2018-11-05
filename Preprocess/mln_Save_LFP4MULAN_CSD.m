function mln_Save_LFP4MULAN_CSD()

dirname='/Users/huifangwang/MULAN3/mlnData/LFP/SlowTheta/';
prenom='20130516_CTR';
Nds=[2,5];
load([dirname,prenom,'.mat'],'eegHpc','eegEC','Par')
csdEC=mln_lfp2csd(eegEC,1,32);
[nt,nchan]=size(eegHpc);
csdHpc=NaN(nt,nchan-2*4);
Nleed=8;
for ip=1:4
    indexip=(ip-1)*Nleed+1:ip*Nleed;
    indexipn=(ip-1)*(Nleed-2)+1:ip*(Nleed-2);
    csdHpc(1:nt,indexipn)=mln_lfp2csd(eegHpc(:,indexip),1,8);
end

oldData=[csdHpc';csdEC'];
ts=0;
periodt=12*60;
fs=Par.lfpSampleRate;

chanchosen=[7:18,25:2:54];
cdata=oldData(chanchosen,ts*fs+1:(ts+periodt)*fs);
str={};
neegH=12;
for ieegH=1:neegH
str{ieegH}=['HD',num2str(chanchosen(ieegH))];
end


for ieegEC=13:26
str{ieegEC}=['ED',num2str(chanchosen(ieegEC)-24)];
end
Params.str=str;

output_directory=[dirname,'CSD/data/'];
if ~exist(output_directory, 'dir')
  mkdir(output_directory);
end
for ids=1:length(Nds)
Data=downsample(cdata',Nds(ids))';
Params.fs=fs/Nds(ids);

save([output_directory,'CSD_',prenom,'ds',num2str(Nds(ids)),'.mat'],'Data','Params')
end



function csd=mln_lfp2csd(lfp,ns,nf)
csd=2*lfp(:,ns+1:nf-1)-lfp(:,ns:nf-2)-lfp(:,ns+2:nf);