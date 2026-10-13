function stats=compareSegments(a,b,startStep,endStep)
%COMPARESEGMENTS Exact state/data/cache/RNG comparison, including NaN identity.
ia=a.data.step>=startStep & a.data.step<=endStep;
ib=b.data.step>=startStep & b.data.step<=endStep;
assert(isequal(a.data.step(ia),b.data.step(ib)));
stats.maxDiff=0; stats.mismatchCount=0; stats.fields=0;
for name=fieldnames(a.data).'
    x=a.data.(name{1}); y=b.data.(name{1});
    if ndims(x)==3 || any(strcmp(name{1},{'inputField','selfFreshField','edgeWeight','age','lastBearing','lastWeight','observedBearing','observedWeight','injectedBearing','trueBearing','trueDistance','selected','fresh','hidden','active'}))
        x=x(:,:,ia); y=y(:,:,ib);
    else
        x=x(ia,:); y=y(ib,:);
    end
    stats=accumulate(stats,x,y);
end
fa=find(ia).'; fb=find(ib).'; fa=[fa(1),fa+1]; fb=[fb(1),fb+1];
for name={'position','heading','u'}
    x=a.state.(name{1}); y=b.state.(name{1});
    if strcmp(name{1},'heading'), x=x(:,fa); y=y(:,fb); else, x=x(:,:,fa); y=y(:,:,fb); end
    stats=accumulate(stats,x,y);
end
stats=accumulate(stats,a.state.alpha,b.state.alpha);
if endStep==a.endStep && endStep==b.endStep
    stats=accumulate(stats,a.finalState,b.finalState);
    stats=accumulate(stats,a.finalRng,b.finalRng);
end
end
function s=accumulate(s,x,y)
s.fields=s.fields+1;
if ~isequaln(x,y)
    s.mismatchCount=s.mismatchCount+1;
    if isnumeric(x)&&isnumeric(y)
        z=abs(x-y); z=z(isfinite(z));
        if ~isempty(z), s.maxDiff=max(s.maxDiff,max(z(:))); end
    end
end
end
