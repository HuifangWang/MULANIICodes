
%dirname='/Users/huifangwang/MULANa/SeegData/PS1001/hidNode/data/';
%dirname='/Users/huifangwang/NormalSEEG_CL/PS1001/mln/data/';
basicname='/Users/huifangwang/CAMULAN/Data/PS1001/mln/tsignals/arr/CA634EOL_1-1+/ND1L20/';
rawdir = [basicname, 'rdata/'];
savedirname = [basicname, 'data/'];
fnames=dir([rawdir,'*.mat']);
numfids = length(fnames);
if ~exist(savedirname,'dir')
    try mkdir(savedirname);
    end
end
flagMR='R';
Nds=4;
for is=1:numfids
    prenom=fnames(is).name;
    mln_downsampling2mln_5(rawdir,savedirname,prenom,flagMR,Nds)
end