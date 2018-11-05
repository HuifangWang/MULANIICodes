%%% this function is main funtion to generate the datasets with random but
%%% few links

nsn=1;
ns=num2str(nsn);
mincs='0.1';
maxcs='1';
nc='30';
slength='20000';
fs='1000';
nl='35';

dirname=['./Examples/data/FRNS',nc,'CS',num2str(1),'S',num2str(nsn),'/'];



for is=1:nsn
    istru=num2str(is);
    prenom='nmm';
    mln_Gen_NMM_fs_vs_ns_few(dirname,prenom, nc,istru, slength, mincs,maxcs, nl,fs)
end