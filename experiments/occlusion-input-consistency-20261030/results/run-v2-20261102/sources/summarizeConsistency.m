function s=summarizeConsistency(GO,LO,window)
% Combine at each time, then average over the fixed terminal window.
assert(isvector(GO)&&isvector(LO)&&numel(GO)==numel(LO));
assert(window>=1&&window<=numel(GO)); ix=(numel(GO)-window+1):numel(GO);
g=GO(ix); l=LO(ix); c=min(g,l);
s.finalConsistency=mean(c); s.finalGO=mean(g); s.finalLO=mean(l);
s.finalConsistencyQ10=quantile(c,.1); s.finalJointHighFraction=mean(g>=.9 & l>=.9);
end
