clear;
addpath("..\Common\")
load("continuous_beta_N_0.53.mat", "xList");
referenceXList = xList;
referenceDeathToll = xList(5,:);
referenceResolving = xList(4,:);
referenceSusceptible = xList(1,:);
load("continuous_beta_N_0.53.mat", "nList");
referenceNList = nList; 
load("continuous_beta_N_0.53.mat", "costList");
costList = costList(costList>0);
referenceCost = min(real(costList));
load("continuous_beta_N_0.53.mat", "betaList");
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
betaN = 0.65;

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
controllerKd = 50000000;
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
    pList(curStep) = du_p;
    dList(curStep) = du_d;
    iList(curStep) = du_i;
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
legend(["Real system betaN =" + string(betaN), "Original system betaN = 0.53", "Controlled system (death toll) betaN =" + string(betaN)])
xlabel("Time (days)")
ylabel("Death toll per 100,000")
figure(2);
plot(timeList, trueNList)
legend(["N for real system betaN =" + string(betaN), "N for original system betaN = 0.53", "Controlled system (death toll) betaN =" + string(betaN)])
xlabel("Time (days)")
ylabel("Actual employment rate")
figure(3);
plot(timeList, pList)
hold on;
plot(timeList, dList)
plot(timeList, iList)
xlabel("Time (days)")
ylabel("Control input")
legend(["Proportional input", "Derivative input", "Integral input"])
controlledCost = costFunctionIntegral(xList(:,1:costFinalStep), betaList(1:costFinalStep), simulationDt);
disp("Cost for the controlled model (death toll) = " + controlledCost);


% approximation of the input with average
xList = zeros(systemDimension, finalStep);
betaList = zeros(1, finalStep);
xList(:,1) = x0;
nList = zeros(1, finalStep);
trueNList = zeros(1, finalStep);


discreteControlPeriod = 14; % control period, in days
discreteControlPeriodStep = fix(discreteControlPeriod / simulationDt);

% initialize for the first control period
for curIndex = 1 : discreteControlPeriodStep
    nList(curIndex) = referenceNList(1);
end
controllerKpDiscrete = 1;
controllerKdDiscrete = 20000;
controllerKiDiscrete = 0;
lastError = 0;
last2Error = 0;

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

        du_p = controllerKpDiscrete * (outputError - lastError);
        du_d = controllerKdDiscrete * ((outputError - lastError) - (lastError - last2Error));
        du_i = controllerKiDiscrete * (outputError);
        du = du_p + du_d + du_i;

        last2Error = lastError;
        lastError = outputError;
        nextN = nList(curStep) + du;
        nextN = max(minN, nextN);
        nextN = min(maxN, nextN);
        for nextIndex = 0 : discreteControlPeriodStep
            if curStep + nextIndex > finalStep
                break;
            end
            nList(curStep + nextIndex) = nextN;
        end
    end

    if curStep < finalStep
        xList(:, curStep + 1) = nextState;
    end
end
figure(1);
plot(timeList, xList(5, :) * 100000);
figure(2);
plot(timeList, trueNList);
figure(4);
plot(timeList, xList(1,:) .* betaList / gamma);
xlabel("Time (days)")
ylabel("Effective R value")

discreteControlCost = costFunctionIntegral(xList(:,1:costFinalStep), betaList(1:costFinalStep), simulationDt);
disp("Cost function with discrete control = " + string(discreteControlCost));
% add the reference optimal value to the figure

load("C:\Users\Klaus\Documents\Graduate\OptimizationRefactored\ContinuousTime\ContinuousSensitivityResultsWithEndo\beta_N_" + string(betaN) + ".mat")
figure(1);
plot(timeList, xList(5,:) * 100000);
legend(["Real system betaN =" + string(betaN), "Original system betaN = 0.53", "Controlled system betaN =" + string(betaN), ...
   "Controlled system (discrete, death toll), betaN =" + string(betaN), "Optimal solution for betaN =" + string(betaN)])
figure(2);
plot(timeList, nList);
legend(["N for real system betaN =" + string(betaN), "N for original system betaN = 0.53", "Controlled system (death toll) betaN =" + string(betaN), ...
    "Controlled system (discrete, death toll), betaN = " + string(betaN), "Optimal N for betaN =" + string(betaN)])
figure(4);
plot(timeList, betaList .* xList(1,:) / gamma);
legend(["Effective R for real system betaN = " + string(betaN), ...
    "Effective R for original system betaN = 0.53", ...
    "Controlled system (discrete, death toll), betaN" + string(betaN), ...
    "Optimal effective R for betaN =" + string(betaN)])


