function mln_testRule()
%data=MisInput;
data=[0.56 0.6 0.75 0.8 0 1 1];
numMFs=2;
inmftype='gbellmf';
outmftype='constant';

threshold=[0.5,0.5,0.75,0.75];
kforMF=5;

methodlog={'BCorrD','PCorrD','GC','oPDCF'};

fis=mln_initialfis(data,methodlog,numMFs, inmftype, outmftype,threshold,kforMF);
mln_ruleview(fis,data,threshold);

