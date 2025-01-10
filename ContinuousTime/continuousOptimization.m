function ret = continuousOptimization(newParams, saveFileName, loadFileName)

    % Define global constants    
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
    if newParams.isKey("muTv")
        muTv = newParams("muTv");
    end
    sigmaTv = 44.74;
    
    % Parameter set 2: beta function (transmission rate)
    global betaW;
    global betaN;
    global betaLambda;
    global lambda;
    global alpha;
    
    betaW = 0.376;
    if newParams.isKey("betaW")
        betaW = newParams("betaW");
    end
    betaN = 0.53;
    if newParams.isKey("betaN")
        betaN = newParams("betaN");
    end
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
    %kappa = 0;
    muF = 245;
    sigmaF = 27.5;
    phiKappa = 0.18;
    
    
    global maxNRIter;
    maxNRIter = 10;
    
    % Define of initial states
    systemDimension = 6; % CONSTANT
    x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
    
    inputDimension = 1; % CONSTANT
    armijoAlpha = 0.01; %def
    %armijoAlpha = 0.1;
    armijoBeta = 0.5;
    
    
    % Define of simulation parameters
    simulationDt = 0.01;
    finalTime = 635;
    if newParams.isKey("finalTime")
        finalTime = newParams("finalTime");
    end
    finalStep = round(finalTime / simulationDt) + 1;
    timeList = 0 : simulationDt : finalTime;
    fixedStepsize = 0.01; %def
    maxIter = 7;
    
    % List of the system states during simulation (time horizon)
    xHist = zeros(systemDimension, finalStep);
    uHist = zeros(1, finalStep);
    nHist = zeros(1, finalStep);
    costHist = zeros(1, maxIter);
    
    
    %%%%%%%%%%%%%%%%%%%%%%% BEGIN OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%
    % Initialize the system states
    xList = zeros(systemDimension, finalStep);
    xList(:,1) = x0;
    % Initialize system costates and inputs
    % Assume the initial input is the minimum, i.e. n = minN for all t.
    % This could somewhat avoid the infeasibility by endogenous response.
    costateList = zeros(systemDimension, finalStep);
    betaList = ones(inputDimension, finalStep);
    %nList = ones(inputDimension, finalStep) * minN; %def
    
    %nList = ones(inputDimension, finalStep) * 0.905; %for betaW=0.35
    nList = ones(inputDimension, finalStep) * 0.84; %for betaW=0.4
    initialNList = nList;
    fminTimer = 0.0;
    nrTimer = 0.0;
    
    %%% continue optimization from previous results
    %load(".\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat", "betaList")
    if loadFileName ~= ""
        load(loadFileName, "nList")
        disp("Successfully loaded initial guess from " + loadFileName)
    end
    %initialNList = movmean(nList, 25000); %for benchmark smoothed
    % initialNList = movmean(nList, 70000); 
    % figure;
    % plot(timeList, initialNList);

    % nList = initialNList;
    
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

    costList = zeros(1, maxIter);
    costList(1) = costFunctionIntegral(xList, betaList, simulationDt);
    disp("Initial cost = " + string(costList(1)) + ", this may come from the previous opt")
    armijoStepsize = 0.5;
    for curIteration = 1:maxIter
        disp(curIteration)
        for curStep = 2 : finalStep + 1
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
        % Calculate the cost function with current input and state
        curCost = costFunctionIntegral(xList, betaList, simulationDt);
        disp(curCost)
        %if curIteration > 1
        %if abs(costList(curIteration) - curCost) < 0.00005 %change is small
        if curCost > costList(curIteration)
            break;
        end
        %end
        costList(curIteration + 1) = curCost;
    
        % Use backward Euler to calculate the costate
        
        for curStep = finalStep - 1 : -1 : 1
            curTime = (curStep) * simulationDt;
            curBeta = betaList(:, curStep + 1);
            curX = xList(:, curStep + 1);
            
            dfDx = getDfDx(curX, curBeta);
            dLDx = getDLDx(curX, curBeta, curTime);
            curP = costateList(:, curStep + 1);
            dP = - dfDx' * curP - dLDx';
            prevP = curP - dP * simulationDt;
            costateList(:, curStep) = prevP;
        end
        
        minimizerList = zeros(1, finalStep);
        
        % Calculate the minimizer of the Hamiltonian
        % Since we don't know the range of beta
        % We have to use the ranges at the last iteration.
        for curStep = 1 : finalStep 
            curTime = (curStep - 1) * simulationDt;
            curP = costateList(:, curStep);
            curX = xList(:, curStep);
            betaRange = getBetaRange(curTime, curX);
        
            target = @(beta)(hamiltonian(curX, beta, curP, curTime));
            minBeta = betaRange(1);
            maxBeta = betaRange(2);
    
            % May use different minimizers to solve the PMP
            curMinimizer = fminbnd(target, minBeta, maxBeta); %this was
            %the default
            %curMinimizer = getMinimizerNewtonRaphson(curX, curP, curTime, minBeta, maxBeta);
            minimizerList(curStep) = curMinimizer;
        end
        % Uncomment this block to use Armijo stepsize
        armijoStepsize = getArmijoStepsize(xList, betaList, costateList, minimizerList, simulationDt, finalTime, armijoAlpha, armijoBeta, ...
               armijoStepsize);
        betaList = betaList + armijoStepsize * (minimizerList - betaList);
        % betaList = betaList + fixedStepsize * (minimizerList - betaList);
    end
    
    
    if saveFileName ~= ""
        save(saveFileName);
    end
    ret = 0;
end

