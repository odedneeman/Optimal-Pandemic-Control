function EffectiveR = getEffectiveR(xList, betaList, tV)
% Returns the series of effective R in each step.
global gamma;
EffectiveR = (betaList / gamma) .* xList(1, :);
EffectiveR(tV*100:end) = 0; % Forcing R-Effective to be 0 after vaccine arrival
end