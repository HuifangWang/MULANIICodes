%Huifang Wang, Jan 17, 14, MULAN parameter validataion.
function mln_rule_Validation(dirname,prenom)

%dirname='nmm';
%prenom='nmmnmmCS100S20N200';
neededfile=[dirname,'/ToutResults','/mlnRes_',prenom,'.mat'];
tosavefile=[dirname,'/ToutResults','/ruleVal_',prenom,'.mat'];
if ~exist(tosavefile,'file')
    
    if exist(neededfile,'file')
        load(neededfile,'param');
        
        G12=param.methodlog.G12;
        G3=param.methodlog.G3;
        G4=param.methodlog.G4;
        th1=param.Thres.th1;
        th2=param.Thres.th2;
        mlnth=param.Thres.mth;
        kforMF=5;
        Simth=0.8;
        NiG1=size(G12,1);
        NiG3=length(G3);
        NiG4=length(G4);
        Nith1=length(th1);
        Nith2=length(th2);
        Nmlnth=length(mlnth);
        nboot1=1000;
        nboot2=100;
        dataprenom=[dirname,'/ToutResults','/Tout_',prenom,'.mat'];
        
        methods_used={G12{:},G3{:},G4{:}};
        rule_list=mln_rule_list();
        [Nrule,nouput]=size(rule_list);
        ctrule_list=rule_list;
        ctrule_list(:,nouput)=1-(rule_list(:,nouput));
        Norder=2;
        RuleScoreAllrule=cell(Nith1,Nith2,Nmlnth,Nrule);
        ctRuleScoreAllrule=cell(Nith1,Nith2,Nmlnth,Nrule);
        %% bootstraping the results for given methods
        Input=GenerateTDM(dataprenom,methods_used,nboot1,nboot2); % Input.values; Input.methods
        load(neededfile,'iMULAN');
        for iG1=1:NiG1
            for iG3=1:NiG3
                for iG4=1:NiG4
                    methodlog=[G12(iG1,:),G3{iG3},G4{iG4}];
                    display([methodlog]);
                    for ith1=1:Nith1
                        for ith2=1:Nith2
                            for imlnth=1:Nmlnth
                                
                                threshold=[th1(ith1),th1(ith1),th2(ith2),th2(ith2)];
                                Ninput=length(methodlog);
                                ind_methods=mln_xorder_cell(Input.methods,methodlog);
                                iInput4methods=Input.values(:,ind_methods,:);
                                mln_Mat=iMULAN.Mat(:,:,iG1,iG3,iG4,ith1,ith2,imlnth);
                                %% condition+output
                                cMat=mln_Mat.*(mln_Mat>Simth);
                                C1=mln_IndirectedLink(cMat,Norder);
                                C2=mln_CommonSource(cMat,Norder);
                                MatC=mln_removediag(cMat);
                                C1=mln_removediag(C1);
                                C2=mln_removediag(C2);
                                
                                res_mulan=[C1 C2 double(~~MatC)];
                                
                                for irule=1:Nrule
%                                     crule=rule_list(irule,:);
%                                     indC=mln_ind_RL(res_mulan,crule(1,end-2:end));
%                                     cInput4methods=iInput4methods(indC,:,:);
%                                     RuleScore=mln_RuleScore(crule,threshold,cInput4methods);
%                                     RuleScoreAllrule{ith1,ith2,imlnth,irule}=[RuleScoreAllrule{ith1,ith2,imlnth,irule};RuleScore(:)];
                                    
                                    ctrule=ctrule_list(irule,:);
                                    indC=mln_ind_RL(res_mulan,ctrule(1,end-2:end));
                                    cInput4methods=iInput4methods(indC,:,:);
                                    ctRuleScore=mln_RuleScore(ctrule,threshold,cInput4methods);
                                    ctRuleScoreAllrule{ith1,ith2,imlnth,irule}=[ctRuleScoreAllrule{ith1,ith2,imlnth,irule};ctRuleScore(:)];
                                    
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        EntropyRuleScore=NaN(Nith1,Nith2,Nmlnth,Nrule);
        medianRuleScore=NaN(Nith1,Nith2,Nmlnth,Nrule);
        for ith1=1:Nith1
                        for ith2=1:Nith2
                            for imlnth=1:Nmlnth
                    for irule=1:Nrule
                        EntropyRuleScore(ith1,ith2,imlnth,irule)=mln_entropy(ctRuleScoreAllrule{ith1,ith2,imlnth,irule},200);
                        medianRuleScore(ith1,ith2,imlnth,irule)=median(ctRuleScoreAllrule{ith1,ith2,imlnth,irule});
                    end
                            end
                        end
        end
        save(tosavefile,'ctRuleScoreAllrule','EntropyRuleScore','medianRuleScore');
                    
        
    end
end

    function Input=GenerateTDM(dataprenom,methods_used,nboot,nboot2)
        load(dataprenom,'Connectivity');
        Nchan=size(Connectivity,1);
        Nlink=Nchan*(Nchan-1);
        Nmethod=length(methods_used);
        Input.values=NaN(Nlink,Nmethod,nboot2);
        Input.methods=methods_used;
        
        nbins=100;
        Norder=2;
        
        %datafile=['Tout_',dataprenom,'.mat'];
        
        MatA=load(dataprenom,methods_used{:});
        for imethod=1:Nmethod
            Mat1=MatA.(methods_used{imethod});
            
            [Nnode,~,Nwin]=size(Mat1);
            isizeboot=Nlink;
            
            MatBoot=NaN( isizeboot*nboot,1);
            
            iMat=abs(mean(Mat1,3));
            iMat(1:size(iMat,1)+1:end) = NaN;
            a=reshape(iMat,Nnode*Nnode,1);
            MatBoot(1:isizeboot,1)=a(~isnan(a));
            
            for iboot=2:nboot
                iBoot_ind=randi(Nwin,[1,Nwin]);
                iBoot_Mat=Mat1(:,:,iBoot_ind);
                iMat=abs(mean(iBoot_Mat,3));
                iMat(1:size(iMat,1)+1:end) = NaN;
                a=reshape(iMat,Nnode*Nnode,1);
                MatBoot((iboot-1)*isizeboot+1:iboot*isizeboot,1)=a(~isnan(a));
            end
            xi=linspace(min(MatBoot),max(MatBoot),nbins);
            n_elements = histc(MatBoot,xi);
            c_elements = cumsum(n_elements);
            c_elements=c_elements./max(c_elements);
            %% seperate MatBoot to MBLink=M(Nlink,)
            %% the means for all windows
            MBLink=reshape(MatBoot,Nlink,nboot);
            for iboot=1:nboot2
                for ipairs=1:Nlink
                    Input.values(ipairs,imethod,iboot)=c_elements(find(xi>=MBLink(ipairs,iboot),1));
                end
            end
            %         %% the means for  N-1 windows
            %         Vjack=1:Nwin;
            %         for ijackknife=1:Nwin
            %             a=Vjack;
            %             a(ijackknife)=[];
            %
            %         lMat=abs(mean(Mat1(:,:,a),3));
            %         lMat(1:size(iMat,1)+1:end) = NaN;
            %         a=reshape(lMat,Nnode*Nnode,1);
            %         pairs=a(~isnan(a));
            %         for ipairs=1:Nlink
            %         iTDM(ipairs+ijackknife*Nlink,imethod)=c_elements(find(xi>=pairs(ipairs),1));
            %         end
            %
            %         end
            %
        end
        % load(datafile,'Connectivity');
        % MatC=Connectivity;
        % C1=MULAN_IndirectedLink(MatC,Norder);
        % C2=MULAN_CommonSource(MatC,Norder);
        % MatC=mln_removediag(MatC);
        % C1=mln_removediag(C1);
        % C2=mln_removediag(C2);
        %
        %         for iwin=1:Nwin+1
        %         iTDM((iwin-1)*Nlink+1:iwin*Nlink,Ninput+3)=MatC./max(MatC);
        %         iTDM((iwin-1)*Nlink+1:iwin*Nlink,Ninput+1)=C1;
        %         iTDM((iwin-1)*Nlink+1:iwin*Nlink,Ninput+2)=C2;
        %         end
        
        
        function iTDMR=mln_data2rules(iTDM,threshold)
            iTDMR=iTDM;
            Nth=length(threshold);
            for icol=1:Nth
                iTDMR(:,icol)=(iTDM(:,icol)>=threshold(icol))+1;
            end
            
            function yvalue=mln_weightbell(threshold,rule,idata)
                k=5;
                Nth=length(threshold);
                value=1;
                for i=1:Nth
                    
                    a = [threshold(i),1-threshold(i)]';
                    b = k*2.*a;
                    c =[0,1]';
                    
                    params=  [a b c];
                    value = gbellmf(idata(i), params(rule(i),:))*value;
                end
                yvalue=value;
                
                function RuleScore=mln_RuleScore(rule,threshold,InputBS)
                    
                    [NindC,~,Nboot]=size(InputBS);
                    RuleScore=NaN(NindC,Nboot);
                    for ilink=1:NindC
                        
                        for iboot=1:Nboot
                            RuleScore(ilink,iboot)=mln_weightbell(threshold,rule(1:4),InputBS(ilink,:,iboot));
                        end
                        
                    end
                    
                    
                    function indC=mln_ind_RL(Cmat,L)
                        % L is the condition (C1,C2,L)
                        NL=length(L);
                        [Nrow,~]=size(Cmat);
                        indC=1:Nrow;
                        for iL=1:NL
                            if L(iL)==-1
                                continue;
                            else
                                cindC=find(Cmat(:,iL)==L(iL));
                                indC=indC(ismember(indC,cindC));
                            end
                        end
                        
