function rule_list=mln_rule_list()
irule=1;
rule_list(irule,:)=[1,1,1,1,-1,-1,0];
%irule=irule+1;
%rule_list(irule,:)=[1,2,1,1,-1,-1,0]; % deleted it for more general cases

irule=irule+1;
rule_list(irule,:)=[1,1,2,2,-1,-1,0];


irule=irule+1;
rule_list(irule,:)=[2,2,1,1,-1,1,0];

irule=irule+1;
rule_list(irule,:)=[2,2,1,1,1,-1,0];

%irule=irule+1;
%rule_list(irule,:)=[2,2,1,1,1,1,0];

irule=irule+1;
rule_list(irule,:)=[2,2,2,2,-1,-1,1];

%irule=irule+1;
%rule_list(irule,:)=[2,2,2,2,1,1,0.5];