function tstfwins=mln_tstf_2wins(switchTimes,Params)
Ndy=length(switchTimes);
tstfwins=NaN(Ndy,2);
tswins=1;
tstfwins(1,1)=tswins;
for iwins=2:Ndy
    tfwins=floor((switchTimes(iwins)*Params.fs-Params.wins*(1-Params.overlap))/(Params.wins*(1-Params.overlap)));
    tstfwins(iwins-1,2)=tfwins;
    tswins=tfwins+1;
    tstfwins(iwins,1)=tswins;
end
    