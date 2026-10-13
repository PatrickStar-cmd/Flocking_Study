function d=decodeNeuralActivity(y,alpha)
d.neuralActivityMass=sum(y,1).';
cx=sum(y.*cos(alpha).',1).'; cy=sum(y.*sin(alpha).',1).';
d.bumpConcentration=zeros(size(cx)); positive=d.neuralActivityMass>0;
d.bumpConcentration(positive)=hypot(cx(positive),cy(positive))./d.neuralActivityMass(positive);
d.rawDecodedAngle=mod(atan2(cy,cx),2*pi);
d.decoded_direction_defined=positive & d.bumpConcentration>1e-8;
d.decodedAnalysisAngle=d.rawDecodedAngle;
d.decodedAnalysisAngle(~d.decoded_direction_defined)=nan;
end
