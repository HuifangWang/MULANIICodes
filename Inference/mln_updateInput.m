function Input_AL=mln_updateInput(Input_AL,MatNet)
Norder=2;
[~,NinputC]=size(Input_AL);
Ninput=NinputC-2;
C1=mln_IndirectedLink(MatNet,Norder);
C2=mln_CommonSource(MatNet,Norder);
C1=mln_removediag(C1);
C2=mln_removediag(C2);
C3=mln_removediag(MatNet');
Input_AL(:,Ninput+1)=~~(C1+C2);
Input_AL(:,Ninput+2)=~~C3;
                
        
