function [mlnMat,sim]=mln_Sta_finalNetwork_Mat_gen(miMat,MatRef,Theta,ncof)

[Nnode,~,Nboots]=size(miMat);
xvalues=linspace(0,1,50);
xind=mln_v_find(xvalues,Theta,1);
mlnMat=NaN(Nnode);
for inode=1:Nnode
    for jnode=1:Nnode
        
        pij=squeeze(miMat(inode,jnode,:));
        [ph,xh]=hist(pij,xvalues);
         if mln_sumhist_more2(ph,xind,Nboots,ncof,'AND')
             mlnMat(inode,jnode)=1;
         else
             mlnMat(inode,jnode)=0;
         end
    end
end

sim=mln_Sim_Cos(MatRef,mlnMat);

function TorF=mln_sumhist_more2(ph,xind,Nboots,ncof,flaglogic)
switch flaglogic
    case 'OR'
        TorF=0;
    case 'AND'
        TorF=1;
end
Nth=length(xind);
for ith=1:Nth
    if sum(ph(xind(ith):end))>Nboots*ncof(ith)
        switch flaglogic
            case 'OR'
                TorF=1;
                return;
        end
    else
        switch flaglogic
            case 'AND'
                TorF=0;
                return;
        end
    end
    
end


function xind=mln_v_find(xvalues,Theta,nN)
Nth=length(Theta);
xind=NaN(Nth,nN);
for ith=1:Nth
    xind(ith,:)=find(xvalues>Theta(ith),nN);
end