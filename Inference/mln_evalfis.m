function [output,IRR,ORR,ARR,CRR]=mln_evalfis(input,fis)

%input=[0.326 0.5];
%condi=[0,0];

% This version only for demonstration
% Dispatch the inputs
%    IRR: the result of evaluating the input values through the membership
%        functions. This is a matrix of size Nr-by-N, where Nr is the number
%        of rules, and N is the number of input variables.
%      * ORR: the result of evaluating the output values through the membership
%        functions. This is a matrix of size NPts-by-Nr*L. The first Nr
%        columns of this matrix correspond to the first output, the next Nr
%        columns correspond to the second output, and so forth.
%      * ARR: the NPts-by-L matrix of the aggregate values sampled at NPts
%        along the output range for each output.
Nrule=length(fis.rule);
IRR=NaN(Nrule,fis.Ninput);
ORR=NaN(Nrule,1);
ARR=NaN(Nrule,1);
CRR=NaN(Nrule,2);
for iinput = 1: fis.Ninput
    fis.input(iinput).value = input(iinput);
    fis.input(iinput).mf_n=length(fis.input(iinput).mf);
    for jmf=1:fis.input(iinput).mf_n
        fis.input(iinput).mf(jmf).value=mln_evalmf(fis.input(iinput).value,fis.input(iinput).mf(jmf).params,fis.input(iinput).mf(jmf).type);
    end
end
rulelist=mln_getfis(fis,'rulelist');

n_rule=length(fis.rule);
rulevalue=zeros(n_rule,1);
mycase=input(fis.Ninput+1:fis.Ninput+fis.Ncondi);
for irule=1:n_rule
    % first check the conditions
    condition=rulelist(irule,fis.Ninput+1:fis.Ninput+fis.Ncondi);
    [isfit,crr1]=isfit_mln_cond_rule(condition,mycase);
    CRR(irule,:)=crr1;
    if ~isfit
        continue;
    end
    antecedent=rulelist(irule,1:fis.Ninput);
    ivalue=1;
    for i_input=1:fis.Ninput
        if antecedent(i_input)>0
            IRR(irule,i_input)=fis.input(i_input).mf(antecedent(i_input)).value;
            ivalue=ivalue*fis.input(i_input).mf(antecedent(i_input)).value;
        end
    end
    rulevalue(irule)=ivalue;
    ARR(irule,1)=ivalue;
end
output_mf=mln_getfis(fis,'outmfparams');
output=rulevalue'*output_mf/sum(rulevalue);
ORR=output_mf;


function [isfit,crr1]=isfit_mln_cond_rule(C1inrule,Cforilink)
isfit=1;
crr1=NaN(2,1);
NC1=length(C1inrule);
for iC1=1:NC1
    if C1inrule(iC1)==-1
        crr1(iC1,1)=1;
        continue;
    else
        if ~C1inrule(iC1)==Cforilink(iC1)
            crr1(iC1,1)=0;
            isfit=0;
        else
            crr1(iC1,1)=1;
        end
    end
end
