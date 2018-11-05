function mln_plt_paper_sup_f2()

x=0:0.001:1;
input=0.6;

figure;

theta=[0.5,0.75];
for itheta=1:length(theta)
ha=subplot(1,2,itheta);    
k=5;
a = [theta(itheta),1-theta(itheta)]';
        b = k*2.*a;
        c =[0,1]';
clh={[0,0,230]/255.,[0,100,0]/255.};

hold on
output=NaN(2,1);
for iLH=1:2
params=[a(iLH),b(iLH),c(iLH)];

y=gbellmf(x, params);
output(iLH)=gbellmf(input,params);
plot(x,y,'LineWidth',2,'color',clh{iLH})
set(ha,'FontSize',18)
end

plot(input,output(1),input,output(2),'marker','.','MarkerSize',30)
output
hold off
end