function calTpr=mln_cal_tfp(iMat,gtMat,Tfr,thetaAUC)
% TprSW the true positive rate for Strong/W connection link
% 3 dimentations, true positive rate for s/W links and true negative rate
% for no connection
%%  for strong link
calTpr.TprS=-1;
calTpr.TprW=-1;
calTpr.Tnr=-1;
try
mlnScoreG=iMat(gtMat>min(Tfr.Gcs));
calTpr.TprS=sum(mlnScoreG>thetaAUC(1))/length(mlnScoreG);
end
%% for weak link 
try
mlnScoreP=iMat(gtMat>min(Tfr.Pcs));
calTpr.TprW=sum(mlnScoreP>=thetaAUC(2))/length(mlnScoreP);
end
%% for True negative link
try
mlnScoreN=iMat(gtMat==0);
calTpr.Tnr=sum(mlnScoreN<thetaAUC(2))/length(mlnScoreN);
end
%% for conditional green link
try
csGS=gtMat(iMat>=thetaAUC(1));
calTpr.ConGS=sum(csGS>=min(Tfr.Gcs))/length(csGS);
calTpr.ConGW=sum((min(Tfr.Pcs)<=csGS) .* (csGS<min(Tfr.Gcs)))/length(csGS);
catch
    calTpr.ConGS=-0.1;
    calTpr.ConGW=-0.1;
end
%% for pink link
try 
csPS=gtMat(logical((thetaAUC(2)<=iMat).* (iMat<thetaAUC(1))));
if length(csPS)==0
    
    calTpr.ConPS=-0.1;
    calTpr.ConPW=-0.1;
else
calTpr.ConPS=sum(csPS>=min(Tfr.Gcs))/length(csPS);
calTpr.ConPW=sum((min(Tfr.Pcs)<=csPS) .* (csPS<min(Tfr.Gcs)))/length(csPS);
end
catch
    calTpr.ConPS=-0.1;
    calTpr.ConPW=-0.1;
end

