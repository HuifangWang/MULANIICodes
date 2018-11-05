%Huifang Wang, Jan, 30, 2015, output the connection matrix of MULAN
function [mdMat,mnMat,fMat]=mln_output_sta_CM(allGMat,mmRule)
    Nchan=size(allGMat,1);
    mdMat=zeros(Nchan,Nchan);
    mnMat=zeros(Nchan,Nchan);
    fMat=zeros(Nchan,Nchan);
    
    for ichan= 1:Nchan
        for jchan= 1:Nchan
            if jchan==ichan
                continue
            end
            
            pij=allGMat(ichan,jchan,:);
            mdMat(ichan,jchan)=median(pij);
            mnMat(ichan,jchan)=mean(pij);
            if mmRule(1,2)<=mdMat(ichan,jchan) && mmRule(1,3)<=mnMat(ichan,jchan)
                fMat(ichan,jchan)=mmRule(1,1);
            else if mmRule(2,2)<=mdMat(ichan,jchan) && mmRule(2,3)<=mnMat(ichan,jchan)
                    fMat(ichan,jchan)=mmRule(2,1);
                else
                    fMat(ichan,jchan)=0;
                end
            end
            
        end
    end