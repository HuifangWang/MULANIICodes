function [FPws,FNws,FPNws,thetaws]=mln_TPN_vcs(iMat,gtMat,wscs)

% this version we deleted the AUC infomration
% June 1, Huifang, H1s and H1w are not overlapped.
% FPws is false positive for CS=0
% FNw is false position for CS>wcs
% FNs is false position for CS>scs
% FPNws is false links =FPws+FNw+FNs

H0=iMat(gtMat==0);
%H1w=iMat(and(gtMat>=wscs(1),gtMat<=wscs(2)));
H1w=iMat(gtMat>=wscs(1));
H1s=iMat(gtMat>=wscs(2));
NH0=length(H0);
NH1w=length(H1w);
NH1s=length(H1s);

[thetaw,FPw,FNw,FPNw]=mln_thetav_ws(H1w,H0);
[thetas,FPs,FNs,FPNs]=mln_thetav_ws(H1s,H0);
        
  thetaws=[thetaw,thetas];
  FPws=[FPw,FPs];
  FNws=[FNw,FNs];
  FPNws=[FPNw,FPNs];

    function [theta,FP,FN,FPN]=mln_thetav_ws(H1,H0)
        FPv=NaN(1,3);
        FNv=NaN(1,3);
        thetav=[min(H1),max(H0),mean([min(H1),max(H0)])];
        for ith=1:3
                FPv(ith)=length(H0(H0>=thetav(ith)));
                FNv(ith)=length(H1(H1<thetav(ith)));
        end
        [FPN,bith]=min(FPv+FNv,[],2);
        FP=FPv(bith);
        FN=FNv(bith);
        theta=thetav(bith);



