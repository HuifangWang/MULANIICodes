function [ts,tf]=mln_wins2time(iwins,Params)
% given iwins, ruturn time involved;
% return start time and end time
ts=(Params.wins*(iwins-1)*(1- Params.overlap))/Params.fs;
tf=ts+Params.wins/Params.fs;