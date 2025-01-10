function [betasWithV, xWithV, nListWithV] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep)
% Calculates new list of states (x) and betas, with values in tV and later
% corresponding to the state with the vaccine (there are only S, C, D and
% max employment)
global maxN;

% Adjust finale step to be consistent for continuous and discrete:
finalStepAdj = finalStep;
if mod(finalStep, 10) == 0
    finalStepAdj = finalStep + 1;
end

tvStep = round(tV / simulationDt);
Etv = xList(2, tvStep);
Itv = xList(3, tvStep);    
Rtv = xList(4, tvStep);    
Cs = Etv + Itv + Rtv + xList(6, tvStep);    
Ds = xList(5, tvStep);    
Ss = xList(1, tvStep);
xSteady = [Ss; 0; 0; 0; Ds; Cs];
betasWithV = betaList;
nListWithV = nList;
xWithV = xList;
for curStep = tvStep+1 : finalStepAdj
    curTime = (curStep - 1) * simulationDt;
    nListWithV(curStep) = maxN;
    betasWithV(curStep) = getBetaFromN(maxN, curTime);
    xWithV(:, curStep) = xSteady;
end
end