clear;
addpath("..\Common\")
load("continuous_beta_N_0.53.mat", "xList", "nList", "betaList", "simulationDt", "finalStep");
global tV;
tV = 540;
%[betaList, xList, nList] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
referenceXList = xList;
referenceDeathToll = xList(5,:);
referenceResolving = xList(4,:);
referenceSusceptible = xList(1,:);
%load("continuous_beta_N_0.53.mat", "nList");
referenceNList = nList; 
load("continuous_beta_N_0.53.mat", "costList");
costList = costList(costList>0);
referenceCost = min(real(costList));
%load("continuous_beta_N_0.53.mat", "betaList");
referenceBetaList = betaList;
costFinalStep = 63501;
% Parameter set 1: model parameters
global sigma;
global gamma;
global theta;
global delta;
sigma = 1/3;
theta = 1/11;
delta = 0.008;
gamma = 1/4;
% Gumbel distribution
global muTv;
global sigmaTv;
muTv = 565.83;
sigmaTv = 44.74;

% Parameter set 2: beta function (transmission rate)
global betaW;
global betaN;
global betaLambda;
global lambda;
global alpha;

betaW = 0.376;
betaN = 0.53;
betaLambda = 0.339;
lambda = 0.12;
alpha = 0.69; 


% Parameter set 3: cost function
global phi;
global chi;
global r;
global w;
phi = 1;
chi = 24000;
r = 0.04/365;
w = 158.9;
%nu = 1/(1.5*365);
%r = r + nu;

% Parameter set 4: input constraints 
global maxN;
global minN;
maxN = 1;
minN = 0.68;

% Parameter set 5: endo response
global kappa;
global muF;
global sigmaF;
global phiKappa;
kappa = 30500;
muF = 245;
sigmaF = 27.5;
phiKappa = 0.18;

% Define of initial states
systemDimension = 6; % CONSTANT
x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
% Define of simulation parameters
simulationDt = 0.01;
finalTime = 635;
finalStep = round(finalTime / simulationDt) + 1;
timeList = 0 : simulationDt : finalTime;

feedbackNList = zeros(1, finalStep);
feedbackNList(1) = referenceNList(1);

% change the parameter
betaN = 0.6;

% first simulate the system 

xList = zeros(systemDimension, finalStep);
xList(:,1) = x0;
betaList = zeros(1, finalStep);
trueNList = zeros(1, finalStep);
for curStep = 1 :finalStep
    curTime = (curStep - 1) * simulationDt;
    curState = xList(:, curStep);

    curKappa = getTVKappa(curTime);
    R = curState(4);
    curEndo = curKappa * delta * theta * R;
    if curEndo > 1 - minN
        curEndo = 1 - minN;
    end
    maxNWithEndo = maxN - curEndo;
    curN = referenceNList(curStep);
    %curN = 1;
    if curN > maxNWithEndo
        curN = maxNWithEndo;
    end
    if curN < minN
        curN = minN;
    end
    trueNList(curStep) = curN;
    
    curBeta = getBetaFromN(curN, curTime);
    betaList(curStep) = curBeta;
    curDynamics = seirdDynamics(curState, curBeta);
    nextState = curState + curDynamics * simulationDt;
    if curStep < finalStep
        xList(:, curStep + 1) = nextState;
    end
end
disp("Cost for the original system = " + referenceCost)
inaccurateCost = costFunctionIntegral(xList(:,1:costFinalStep), betaList(1:costFinalStep), simulationDt);
disp("Cost for the true system = " + inaccurateCost)
figure(1);
plot(timeList, xList(5,:) * 100000)
hold on;
plot(timeList, referenceXList(5,:) * 100000)


figure(2);
plot(timeList, trueNList)
hold on;
plot(timeList, referenceNList)
figure(3);
hold on;
figure(4);
hold on;
plot(timeList, xList(1,:) .* betaList / gamma);
plot(timeList, referenceSusceptible .* referenceBetaList / gamma);



% test the controller
xList = zeros(systemDimension, finalStep);
betaList = zeros(1, finalStep);
xList(:,1) = x0;
trueNList = zeros(1, finalStep);
nList = zeros(1, finalStep);
nList(1) = referenceNList(1);
controllerKp = 10000;
controllerKd = 80000000;
controllerKi = 1;
lastError = 0;
last2Error = 0;
dList = zeros(1, finalStep);
pList = zeros(1, finalStep);
iList = zeros(1, finalStep);
for curStep = 1 :finalStep
    curTime = (curStep - 1) * simulationDt;
    curState = xList(:, curStep);

    curKappa = getTVKappa(curTime);
    R = curState(4);
    curEndo = curKappa * delta * theta * R;
    if curEndo > 1 - minN
        curEndo = 1 - minN;
    end
    maxNWithEndo = maxN - curEndo;
    curN = nList(curStep);
    if curN > maxNWithEndo
        curN = maxNWithEndo;
    end
    trueNList(curStep) = curN;

    curBeta = getBetaFromN(curN, curTime);
    betaList(curStep) = curBeta;
    curDynamics = seirdDynamics(curState, curBeta);
    nextState = curState + curDynamics * simulationDt;
    
    curDeathToll = curState(5);
    curDeathTollRef = referenceDeathToll(curStep);
    outputError = curDeathTollRef - curDeathToll;
    du_p = controllerKp * (outputError - lastError);
    du_d = controllerKd * ((outputError - lastError) - (lastError - last2Error));
    du_i = controllerKi * (outputError);
    du = du_p + du_d + du_i;
    
    last2Error = lastError;
    lastError = outputError;
    if curStep < finalStep
        xList(:, curStep + 1) = nextState;
        nList(curStep + 1) = nList(curStep) + du;
        nList(curStep + 1) = min(maxN, nList(curStep +1));
    end
end
figure(1);
plot(timeList, xList(5,:)* 100000)
xlabel("Time (days)")
ylabel("Death toll per 100,000")
figure(4);
plot(timeList, xList(1,:) .* betaList / gamma);
figure(2);
plot(timeList, trueNList)
xlabel("Time (days)")
ylabel("Actual employment rate")

controlledCost = costFunctionIntegral(xList(:,1:costFinalStep), betaList(1:costFinalStep), simulationDt);
disp("Cost for the controlled model (death toll) = " + controlledCost);


% approximation of the input with average
xList = zeros(systemDimension, finalStep);
betaList = zeros(1, finalStep);
xList(:,1) = x0;
nList = zeros(1, finalStep);
pidNList = zeros(1, finalStep);
trueNList = zeros(1, finalStep);


discreteControlPeriod = 1; % control period, in days
discreteControlPeriodStep = fix(discreteControlPeriod / simulationDt);

% initialize for the first control period
for curIndex = 1 : discreteControlPeriodStep
    nList(curIndex) = referenceNList(1);
    pidNList(curIndex) = referenceNList(1);
end
controllerKpDiscrete = 1000;
controllerKdDiscrete = 300000;
controllerKiDiscrete = 1;
lastError = 0;
last2Error = 0;

errorThreshold = 0.000008;
inputThreshold = 0.05; 
dwellThreshold = 14;
policyChangeList = 0;
% run simulation
for curStep = 1 :finalStep
    curTime = (curStep - 1) * simulationDt;
    curState = xList(:, curStep);

    curKappa = getTVKappa(curTime);
    R = curState(4);
    curEndo = curKappa * delta * theta * R;
    if curEndo > 1 - minN
        curEndo = 1 - minN;
    end
    maxNWithEndo = maxN - curEndo;
    curN = nList(curStep);
    if curN > maxNWithEndo
        curN = maxNWithEndo;
    end
    if curN < minN
        curN = minN;
    end
    trueNList(curStep) = curN;
    curBeta = getBetaFromN(curN, curTime);
    betaList(curStep) = curBeta;
    curDynamics = seirdDynamics(curState, curBeta);
    nextState = curState + curDynamics * simulationDt;
    if mod(curStep, discreteControlPeriodStep) == 0
        curDeathToll = curState(5);
        curDeathTollRef = referenceDeathToll(curStep);
        outputError = curDeathTollRef - curDeathToll;

        %curSusceptible = curState(1);
        %curEffectiveR = curBeta / gamma * curSusceptible;
        %curSusceptibleRef = referenceSusceptible(curStep);
        %curBetaRef = referenceBetaList(curStep);
        %curEffectiveRRef = curBetaRef / gamma * curSusceptibleRef;
        %outputError = curEffectiveRRef - curEffectiveR;
        

        du_p = controllerKpDiscrete * (outputError - lastError);
        du_d = controllerKdDiscrete * ((outputError - lastError) - (lastError - last2Error));
        du_i = controllerKiDiscrete * (outputError);
        du = du_p + du_d + du_i;
       
        pList(curStep) = du_p;
        dList(curStep) = du_d;
        iList(curStep) = du_i;
        last2Error = lastError;
        lastError = outputError;
        
        nextN = pidNList(curStep) + du;
        nextN = max(minN, nextN);
        nextN = min(maxN, nextN);
        for nextIndex = 0 : discreteControlPeriodStep
            if curStep + nextIndex > finalStep
                break;
            end
            pidNList(curStep + nextIndex) = nextN;
        end
        
        for nextIndex = 0 : discreteControlPeriodStep
            if curStep + nextIndex > finalStep
                break;
            end
            if abs(outputError) >= errorThreshold && ...
                    abs(nList(curStep) - pidNList(curStep)) >= inputThreshold && ...
                curStep - policyChangeList(end) >= dwellThreshold / simulationDt
                nList(curStep + nextIndex) = pidNList(curStep + nextIndex);
                policyChangeList = [policyChangeList, curStep];
            else
                nList(curStep + nextIndex) = nList(curStep);
            end
        end
        
    end

    if curStep < finalStep
        xList(:, curStep + 1) = nextState;
    end
end
policyChangeGap = [];
for curIndex = 2 : length(policyChangeList)
    gap = policyChangeList(curIndex) - policyChangeList(curIndex - 1);
    policyChangeGap = [policyChangeGap, gap * simulationDt];
end
disp(policyChangeGap)
discreteXList = xList;
discreteNList = trueNList;
discreteBetaList = betaList;
figure(1);
plot(timeList, xList(5, :) * 100000);
figure(2);
plot(timeList, trueNList);
figure(3);
plot(pList)
hold on;
plot(dList)
plot(iList)
xlabel("Time (days)")
ylabel("Control input")
legend(["Proportional input", "Derivative input", "Integral input"])
figure(4);
plot(timeList, xList(1,:) .* betaList / gamma);
xlabel("Time (days)")
ylabel("Effective R value")

discreteControlCost = costFunctionIntegral(xList(:,1:costFinalStep), betaList(1:costFinalStep), simulationDt);
disp("Cost function with discrete control = " + string(discreteControlCost));
% add the reference optimal value to the figure

load("..\ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_" + string(betaN) + ".mat")
figure(1);
plot(timeList, xList(5,:) * 100000);
legend(["Real system \beta_N =" + string(betaN), "Original system \beta_N = 0.53", "Controlled system (continuous, death toll), \beta_N = " + string(betaN), ...
   "Controlled system (discrete, death toll), \beta_N = " + string(betaN), "Optimal death for \beta_N = " + string(betaN)])
figure(2);
plot(timeList, nList);
legend(["Real system \beta_N = " + string(betaN), "Original system \beta_N = 0.53", "Controlled system (continuous, death toll), \beta_N = " + string(betaN), ...
    "Controlled system (discrete, death toll), \beta_N = " + string(betaN), "Optimal employment for \beta_N =" + string(betaN)])
figure(4);
plot(timeList, betaList .* xList(1,:) / gamma);
legend(["Real system \beta_N = " + string(betaN), ...
    "Original system \beta_N = 0.53", ...
    "Controlled system (continuous, death toll) \beta_N = " + string(betaN), ...
    "Controlled system (discrete, death toll) \beta_N = " + string(betaN),...
    "Optimal effective R for \beta_N = " + string(betaN)])

xList = discreteXList;
betaList = discreteBetaList;
nList = discreteNList;
fileName = "./matResults/betaNDeathToll_" + betaN + ".mat";
save(fileName);


