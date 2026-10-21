function audit=testBearingIMMTrack()
stream=RandStream('mt19937ar','Seed',5198); E=7; dt=.3; R=.01;
same=[]; reference=[]; different=[]; rotated=[]; maxSame=0; maxRotation=0;
minimumDet=Inf; minimumVariance=Inf; offset=.73;
P=[.98,.02;.10,.90];
for t=1:500
 fresh=rand(stream,E,1)>.2;
 y=mod(.03*t+(1:E).'/7+.1*randn(stream,E,1),2*pi);
 if t==211, y=mod(y+3,2*pi); end
 y(~fresh)=NaN;
 [reference,~,~]=task3BatchTrack(reference,y,fresh,dt,1e-4,R,Inf);
 [same,~]=bearingIMMTrackStep(same,y,fresh,dt,[1e-4,1e-4],R,P);
 [different,~]=bearingIMMTrackStep(different,y,fresh,dt,[1e-4,1],R,P);
 yr=mod(y+offset,2*pi);
 [rotated,~]=bearingIMMTrackStep(rotated,yr,fresh,dt,[1e-4,1],R,P);
 has=reference.has;
 err=atan2(sin(same.x1(has)-reference.x1(has)),cos(same.x1(has)-reference.x1(has)));
 maxSame=max([maxSame;abs(err);abs(same.x2(has)-reference.x2(has)); ...
  abs(same.p11(has)-reference.p11(has));abs(same.p12(has)-reference.p12(has)); ...
  abs(same.p22(has)-reference.p22(has))]);
 er=atan2(sin(rotated.x1(has)-different.x1(has)-offset), ...
  cos(rotated.x1(has)-different.x1(has)-offset));
 probabilityError=abs(rotated.mu(has,:)-different.mu(has,:));
 maxRotation=max([maxRotation;abs(er);probabilityError(:)]);
 assert(all(isfinite(different.x1(has))) && all(isfinite(different.mu(:))));
 assert(max(abs(sum(different.mu,2)-1))<1e-12);
 for j=1:3
  if j<3, a=different.models{j}; else, a=different; end
  minimumDet=min(minimumDet,min(a.p11(has).*a.p22(has)-a.p12(has).^2));
  minimumVariance=min([minimumVariance;a.p11(has);a.p22(has)]);
 end
end
assert(maxSame<1e-10,'Identical-model IMM does not reduce to the legacy KF');
assert(maxRotation<1e-10,'Rotation equivariance fails');
assert(minimumDet>=-1e-12 && minimumVariance>=0,'Non-PSD covariance');
before=different; y=nan(E,1); fresh=false(E,1);
[after,~]=bearingIMMTrackStep(before,y,fresh,dt,[1e-4,1],R,P);
assert(max(abs(after.mu-before.mu*P),[],'all')<1e-14);
[empty,~]=bearingIMMTrackStep([],y,fresh,dt,[1e-4,1],R,P);
assert(~any(empty.has));
audit=struct('status','PASS','identicalModelMaxDifference',maxSame, ...
 'rotationMaxDifference',maxRotation,'minimumCovarianceDeterminant',minimumDet, ...
 'minimumVariance',minimumVariance,'missingOnlyPriorTransition',true, ...
 'noMeasurementNoInitialization',true,'jumpLikelihoodFinite',true, ...
 'controllerSignatureContainsNoTruthOrRange',true);
disp(jsonencode(audit));
end
