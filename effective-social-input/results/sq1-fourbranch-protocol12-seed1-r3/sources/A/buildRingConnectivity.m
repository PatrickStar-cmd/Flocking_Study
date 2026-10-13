function W = buildRingConnectivity(nNeurons, exponent)
%BUILDRINGCONNECTIVITY Construct the modified cosine ring connectivity.

validateattributes(nNeurons, {'numeric'}, {'scalar','integer','>=',3});
validateattributes(exponent, {'numeric'}, {'scalar','positive'});

alpha = linspace(0, 2*pi, nNeurons + 1);
alpha = alpha(1:end-1);
delta = abs(alpha(:) - alpha(:).');
delta = min(delta, 2*pi - delta);
W = cos((delta/pi).^exponent * pi);
end

