% Generate the initial ML_fis structure
% Huifang Wang oct. 17,13 update rules
% Huifang Wang, April 1, 15 update rules
function fis=mln_initialfis(data, methodlog,numMFs, inmftype, outmftype,para,kforMF)
threshold=para.pc(1:4);
[Ndatrow,Ndatcol]=size(data);
NinputG1=length(methodlog);
NinputG2=2;
fis.Ninput=NinputG1;
fis.Ncondi=NinputG2;
default_numMFs=2;
default_inmftype='gbellmf';
default_outmftype='constant';
fis.kforMF=kforMF;
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
fis.name='Mulan_fis';


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

%range = [min(data,[],1)' max(data,[],1)'];
range =repmat([0 1],7,1);
in_mf_param = mln_genparam(threshold, numMFs, inmftype,kforMF);

%% initialize the input
k=1;
for i = 1:NinputG1,
 fis.input(i).name = methodlog{i};
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
%rule_n=20;
fis.output(1).name='output';

fis.output(1).range=range(end,:); 
%   if strcmp(outmftype, 'linear')
%    fis.output(1).mf(i).params=zeros(1, in_n+1);
%   else
   
%   end
% end

%% initialize the rule
rule_list = mln_rule_list();
indr1=find(strcmp(para.name,'r1'));
%rule_list(:,end)=1;
Nrule=size(rule_list,1);
for irule=1:Nrule

%rule_list(irule,:)=[1,1,1,1,0,0,1];
fis.output(1).mf(irule).params=para.pc(indr1+irule-1);
fis.output(1).mf(irule).range=0;
fis.output(1).mf(irule).name=['out1mf', num2str(irule)];  
fis.output(1).mf(irule).type=outmftype;
end

rule_list(:,end)=1;
fis.condit.name={'C1','C2'};
%rule_list = [rule_list  ones(rule_n, 1) ones(rule_n, 1)];
fis.rule=[];
fis=mln_setfis(fis, 'rulelist', rule_list);

if length(fis.rule)> 250
    wmsg = sprintf('genfis1 has created a large rulebase in the FIS. \nMATLAB may run out of memory if this FIS is tuned using ANFIS.\n');
    warning(wmsg);
end



