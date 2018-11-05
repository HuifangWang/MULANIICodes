function trfGP=mln_find_thetaGS(iMat,gtMat,Tfr)
thetaW=0.9:0.001:0.999;
NthW=length(thetaW);
TprS=ones(NthW,NthW).*-0.1;
TprW=TprS;
Tnr=TprS;
ConGS=Tnr;
ConGW=Tnr;
ConPS=Tnr;
ConPW=Tnr;
for ithetaW=1:NthW
    for ithetaS=ithetaW:NthW
        trfGP=mln_cal_tfp(iMat,gtMat,Tfr,[thetaW(ithetaS),thetaW(ithetaW)]);
        TprS(ithetaW,ithetaS)=trfGP.TprS;
        TprW(ithetaW,ithetaS)=trfGP.TprW;
        Tnr(ithetaW,ithetaS)=trfGP.Tnr;
        ConGS(ithetaW,ithetaS)=trfGP.ConGS;
        ConGW(ithetaW,ithetaS)=trfGP.ConGW;
        ConPW(ithetaW,ithetaS)=trfGP.ConPW;
        ConPS(ithetaW,ithetaS)=trfGP.ConPS;
    end
end
figure
subplot(2,4,1)
imagesc(thetaW,thetaW,TprS)
title('TprS')
colorbar
subplot(2,4,2)
imagesc(thetaW,thetaW,TprW)
colorbar
title('TprW')
subplot(2,4,3)
imagesc(thetaW,thetaW,Tnr)
colorbar
title('Tnr')
subplot(2,4,5)
imagesc(thetaW,thetaW,ConGS)
title('ConGS')
colorbar
subplot(2,4,6)
imagesc(thetaW,thetaW,ConGW)
colorbar
title('ConGW')
subplot(2,4,7)
imagesc(thetaW,thetaW,ConPW)
colorbar
title('ConPW')
subplot(2,4,8)
imagesc(thetaW,thetaW,ConPS)
colorbar
title('ConPS')
