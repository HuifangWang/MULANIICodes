
function varargout=mln_AUC(M1,RM,issystem,isample) 
 if nargin<3
   issystem=0;
 end
 if nargin<4
     isample=1;
 end
    [~,~,~,~,~,~,~,~,auc,~,thcuf] =mln_calc_FalseRate(M1,RM,issystem,isample);
    varargout{1}=auc;
    varargout{2}=thcuf(1);
    