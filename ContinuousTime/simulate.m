

% Here we run simulations for many parameters values, either for betaW,
% betaN, or employment. In the case of employment, the simulation is for
% fixed employment policies. In the case of betaWs, the policy is the
% baseline optimal policy.

clear;    

paramType = 1; % 1 - N (employment policy), 2 - betaW 3 - betaN
fileSuff = ""; % if you want to add a suffix to the file name

global betaW
global betaN
fixedNList = 0.68:0.005:1;
betaWList = 0.3:0.002:0.45;
betaNList = 0.4:0.0025:0.6;
paramString = "N";
betaW = 0.376;
betaN = 0.4;

if paramType == 2
    load("..\DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat");
    nInput = nList;
    paramString = "betaW";
end

if paramType == 3
    load("..\DiscreteTime\DiscreteSensitivityResultsWithEndo\beta_N_0.53.mat");
    nInput = nList;
    paramString = "betaN";
end

optionList = {fixedNList, betaWList, betaNList};
parameterList = optionList{paramType};
statNames = {"Death Toll", "GDP Loss", "Harm", "Span", "J Value"};
statsPerParam = zeros(numel(statNames) + 1, numel(parameterList));
statsPerParam(1, :) = parameterList;


for i = 1:numel(parameterList)
    param = parameterList(i);

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
    % if newParams.isKey("muTv")
    %     muTv = newParams("muTv");
    % end
    sigmaTv = 44.74;
    tV = 540;
    % Parameter set 2: beta function (transmission rate)
    global betaN;
    global betaLambda;
    global lambda;
    global alpha;
    
    if paramType == 2
        betaW = param;
    end
    if paramType == 3
        betaN = param;
    end
    % if newParams.isKey("betaN")
    %     betaN = newParams("betaN");
    % end
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
    % To disable the endogenous response, set kappa = 0. 
    kappa = 30500;
    muF = 245;
    sigmaF = 27.5;
    phiKappa = 0.18;
    
    
    global maxNRIter;
    maxNRIter = 10;
    
    % Define of initial states
    systemDimension = 6; % CONSTANT
    x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
    
    inputDimension = 1; % CONSTANT
    armijoAlpha = 0.01;
    armijoBeta = 0.5;
    
    
    % Define of simulation parameters
    simulationDt = 0.01;
    finalTime = 635;
    % if newParams.isKey("finalTime")
    %     finalTime = newParams("finalTime");
    % end
    finalStep = round(finalTime / simulationDt) + 1;
    tvStep = round(tV / simulationDt);
    timeList = 0 : simulationDt : finalTime;
    fixedStepsize = 0.01;
    
    % List of the system states during simulation (time horizon)
    xHist = zeros(systemDimension, finalStep);
    uHist = zeros(1, finalStep);
    nHist = zeros(1, finalStep);
    
    % Initialize the system states
    xList = zeros(systemDimension, finalStep);
    xList(:,1) = x0;
    % Initialize system costates and inputs
    % Assume the initial input is the minimum, i.e. n = minN for all t.
    % This could somewhat avoid the infeasibility by endogenous response.
    costateList = zeros(systemDimension, finalStep);
    betaList = ones(inputDimension, finalStep);

    if paramType == 1
        nInput = ones(inputDimension, finalStep) * param;
    end

    nList = nInput;
    fminTimer = 0.0;
    nrTimer = 0.0;
   

     % initialize beta
    for curStep = 1 : finalStep
        curTime = (curStep - 1) * simulationDt;
        curBeta = getBetaFromN(nList(curStep), curTime);
        betaList(curStep) = curBeta;
    end
    % initialize x
    for curStep = 2 : finalStep
        curTime = (curStep - 2) * simulationDt;
        curX = xList(:, curStep - 1);
        curBeta = betaList(curStep - 1);
        % calculate the value of beta corresponding to the endo
        betaRange = getBetaRange(curTime, curX);
        curBeta = max(curBeta, betaRange(1));
        curBeta = min(curBeta, betaRange(2));
        betaList(curStep - 1) = curBeta;
        nList(curStep - 1) = getNFromBeta(curBeta, curTime);
        dX = seirdDynamics(curX, curBeta);
        nextX = curX + dX * simulationDt;
        if curStep <= finalStep
            xList(:, curStep) = nextX;
        end
    end
    
    [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    % allList = (vertcat(timeList, xListVac, betaListVac, nListVac))';
    deathToll = xListVac(5, tvStep) * 100000; % Death toll
    GDPLoss = 100 * getGDPLoss(xListVac, nListVac, simulationDt, finalTime); % GDP loss
    span = 100 * (1- xListVac(1, tvStep)); %  Span of the disease
    harm = 100 * getEquivalentLoss(xListVac, betaListVac, simulationDt, finalStep) / finalTime; % Harm
    jValue = costFunctionIntegral(xList, betaList, simulationDt); % J
    statsPerParam(2:6, i) = [deathToll, GDPLoss, harm, span, jValue];
    disp("param = " + param + " done")
end

statsPerParam = real(statsPerParam);
save("statsPer" + paramString + "_betaN_" + betaN + ".mat", "statsPerParam")


