function fis=initialMLNfisV1(data, numMFs, inmftype, outmftype,threshold,kforMF)

[Ndatrow,Ndatcol]=size(data);
NinputG1=Ndatcol-3;
NinputG2=2;
fis.Ninput=NinputG1;
fis.Ncondi=NinputG2;
default_numMFs=2;
default_inmftype='gbellmf';
default_outmftype='constant';
if nargin <= 3,
	outmftype = default_outmftype;
end
if nargin <= 2,
	inmftype = default_inmftype;
end
if nargin <= 1,
	numMFs = default_numMFs;
end
 if length(numMFs)==1,
   numMFs=numMFs*ones(1, NinputG1);
 end
fis.name='MLN_fis_Hn1';


fis.type = 'sugeno';

fis.andMethod = 'prod';
fis.orMethod = 'max';
fis.defuzzMethod = 'wtaver';
fis.impMethod = 'prod';
fis.aggMethod = 'max';
if size(inmftype, 1) ==1 &&  NinputG1>1
   for i=2:NinputG1
      inmftype(i,:)=inmftype(1,:);
   end
end

range = [min(data,[],1)' max(data,[],1)'];
in_mf_param = MLN_genparam(threshold, numMFs, inmftype,kforMF);

%% initialize the input
k=1;
for i = 1:NinputG1,
 fis.input(i).name = ['CCM' num2str(i)];
 fis.input(i).range=range(i,:);
 for j=1:numMFs(i)
  MFType = deblank(inmftype(i, :));
  fis.input(i).mf(j).name = ['in' num2str(i) 'mf' num2str(j)];
  fis.input(i).mf(j).type = MFType;
  if strcmp(MFType,'gaussmf') | strcmp(MFType,'sigmf') ...
                | strcmp(MFType,'smf'),
     fis.input(i).mf(j).params= in_mf_param(k,1:2);
  elseif strcmp(MFType,'trimf') | strcmp(MFType,'gbellmf'),
     fis.input(i).mf(j).params= in_mf_param(k,1:3);
  else
     fis.input(i).mf(j).params= in_mf_param(k,1:4);
  end  
  k=k+1;
 end
end

%% initiailize the output
rule_n=20;
fis.output(1).name='output';

fis.output(1).range=range(end,:); 
%   if strcmp(outmftype, 'linear')
%    fis.output(1).mf(i).params=zeros(1, in_n+1);
%   else
   
%   end
% end

%% initialize the rule
rule_list = zeros(rule_n, NinputG1+NinputG2+1);

irule=1;
rule_list(irule,:)=[1,1,1,1,0,0,1];
fis.output(1).mf(irule).params=0;
fis.output(1).mf(irule).range=0;

irule=irule+1;
rule_list(irule,:)=[1,2,1,1,0,0,1];
fis.output(1).mf(irule).params=0;
fis.output(1).mf(irule).range=[0, 0.1];

irule=irule+1;
rule_list(irule,:)=[1,1,2,2,0,0,1];
fis.output(1).mf(irule).params=0;
fis.output(1).mf(irule).range=[0, 0.2];

irule=irule+1;
rule_list(irule,:)=[1,0,1,1,1,0,1];
fis.output(1).mf(irule).params=0;
fis.output(1).mf(irule).range=[0, 0.2];

irule=irule+1;
rule_list(irule,:)=[2,0,1,1,1,0,1];
fis.output(1).mf(irule).params=0;
fis.output(1).mf(irule).range=[0, 0.2];

irule=irule+1;
rule_list(irule,:)=[2,0,1,1,0,0,1];
fis.output(1).mf(irule).params=0.1;
fis.output(1).mf(irule).range=[0, 0.2];

irule=irule+1;
rule_list(irule,:)=[2,0,1,1,0,1,1];
fis.output(1).mf(irule).params=0;
fis.output(1).mf(irule).range=[0, 0.2];

irule=irule+1;
rule_list(irule,:)=[2,0,2,2,1,0,1];
fis.output(1).mf(irule).params=0.9;
fis.output(1).mf(irule).range=[0.5, 1];

irule=irule+1;
rule_list(irule,:)=[2,0,2,2,0,1,1];
fis.output(1).mf(irule).params=0.9;
fis.output(1).mf(irule).range=[0.5, 1];

irule=irule+1;
rule_list(irule,:)=[2,0,2,0,0,0,1];
fis.output(1).mf(irule).params=0.8;
fis.output(1).mf(irule).range=[0.5, 1];


irule=irule+1;
rule_list(irule,:)=[2,2,2,2,0,0,1];
fis.output(1).mf(irule).params=1;
fis.output(1).mf(irule).range=1;

rule_n=irule;
rule_list=rule_list(1:rule_n,:);

for i = 1:rule_n,
  fis.output(1).mf(i).name=['out1mf', num2str(i)];  
 fis.output(1).mf(i).type=outmftype;
end



%rule_list = [rule_list  ones(rule_n, 1) ones(rule_n, 1)];
fis.rule=[];
fis=MLN_setfis(fis, 'rulelist', rule_list);

if length(fis.rule)> 250
    wmsg = sprintf('genfis1 has created a large rulebase in the FIS. \nMATLAB may run out of memory if this FIS is tuned using ANFIS.\n');
    warning(wmsg);
end



