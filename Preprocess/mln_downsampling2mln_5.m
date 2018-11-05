function mln_downsampling2mln_5(dirname,savedir,prenom,flagMR,Nds)

% downsampling and save to mln format
%sname='nmmN30L10CS100S5N40100';
%dirname='/Users/huifangwang/MULAN3/mlnData/SEEG/Patient00/';
%prenom='Patient00P3ch53.mat';
%dirname='/Users/huifangwang/MULANa/SeegData/PS1001/data/';
%prenom='WB_PS1001WI041222014ts0.mat';
%flagMR='R';
%Nds='4'

%Nds=str2num(Nds);
indmat=strfind(prenom,'.mat');
if ~isempty(indmat)
    sname=prenom(1:indmat-1);
else
    sname=prenom;
end

filename=[dirname,'/',sname,'.mat'];
old=load(filename);
cdata=old.Data;


Params=old.Params;
Params.fs=Params.fs/Nds;
sfile=[savedir,'/',sname,'ds',num2str(Nds),'.mat'];

Data=downsample(cdata',Nds)';
if flagMR=='M'
    Connectivity=old.Connectivity;
    save(sfile,'Data','Params','Connectivity');
else
    save(sfile,'Data','Params');
end
