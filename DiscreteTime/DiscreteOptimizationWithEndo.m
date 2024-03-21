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
    sigmaTv = 44.74;
    
    % Parameter set 2: beta function (transmission rate)
    global betaW;
    global betaN;
    global betaLambda;
    global lambda;
    global alpha;
    
    betaW = 0.376;
    betaN = 0.53;
    if newParam ~= 0
        betaW = newParam;
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
    kappa = 30500;
    muF = 245;
    sigmaF = 27.5;
    phiKappa = 0.18;
    
    
    global maxNRIter;
    maxNRIter = 10;
    
    
    x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
    simulationDt = 0.01;
    finalTime = 675;
    finalStep = round(finalTime / simulationDt);
    timeList = 0 : simulationDt : finalTime;
    fixedStepsize = 0.01;
    maxIter = 90;
    numericalVerification = 0;
    
    
    
    xList = zeros(6, finalStep + 1); 
    xList(:,1) = x0;
    % List of costates
    pList = zeros(6, finalStep + 1);
    
    betaList = ones(1, finalStep + 1);
    nList = ones(1, finalStep + 1);
    
    
    t1 = [200, 450];
    %t1 = [38,  224];
    n1 = [0.6, 0.7, 0.8];
    costList = zeros(1, maxIter);
    
    dt1List = zeros(1, maxIter);
    dt1NumericalList = zeros(1, maxIter);
    stepSize = 1;
    armijoBeta = 0.5;
    armijoAlpha = 0.1;
    armijoIncrease = 0;
    stepsizeList = zeros(1, maxIter);
    
    inputArmijoBeta = 0.5;
    inputArmijoStepsize = 0.001;
    inputArmijoAlpha = 0.01;
    inputArmijoIncrease = 0;
    inputStepsizeList = zeros(1, maxIter);
    
    
    zigzagCount = 0;
    zigzagState = 0;
    for curIteration = 1 : maxIter
        disp(n1)
        disp(t1)
        zigzagCount = zigzagCount + 1;
        if zigzagCount == 30
            zigzagCount = 0;
            zigzagState = 1 - zigzagState;
        end
        
        disp(curIteration)
        disp(zigzagState)
       
        % First step: calculate the state according to t1
        stage = 1;
        for curStep = 1 : finalStep + 1
            % The input to the system should be adjusted according to
            % endogenous response.
            curTime = (curStep - 1) * simulationDt;
            if stage <= length(t1) && curTime >= t1(stage)
                stage = stage + 1;
            end
            
            curState = xList(:, curStep);
            
            % calculate endogenous response
            curKappa = getTVKappa(curTime);
            R = curState(4);
            curEndo = curKappa * delta * theta * R;
            maxNWithEndo = maxN - curEndo;
            
            % adjust the input according to endogenous response
            curN = n1(stage);
            if curN > maxNWithEndo
                curN = maxNWithEndo;
            end
            if curN < minN
                curN = minN;
            end
            
            nList(curStep) = curN;
            curBeta = getBetaFromN(curN, curTime);
            betaList(curStep) = curBeta;
            if curStep <= finalStep
                curDynamics = seirdDynamics(curState, curBeta);
                nextState = curState + curDynamics * simulationDt;
                xList(:, curStep + 1) = nextState;
            end
        end
        
    
        curCost = costFunctionIntegral(xList, betaList, simulationDt);
        curCost = real(curCost);
        disp(curCost);
        costList(curIteration) = curCost;
        for curStep = finalStep + 1 : -1: 2
            curTime = (curStep - 1) * simulationDt;
            curP = pList(:, curStep);
            curX = xList(:, curStep);
            curBeta = betaList(:, curStep); % curStep - 1 here? 
    
            dfDx = getDfDx(curX, curBeta);
            dLDx = getDLDx(curX, curBeta, curTime);
    
            dPDt = - dfDx' * curP - dLDx';
            
            nextP = curP - dPDt * simulationDt;
            pList(:, curStep - 1) = nextP;
        end
       
        if zigzagState == 0
        
            % calculate the derivative of the cost with respect to the input
            dJDn1 = zeros(1, length(n1));
            startTime = 0;
            for curStage = 1:length(t1) + 1
                if curStage ~= length(t1) + 1
                    endTime = t1(curStage);
                else
                    endTime = finalTime;
                end
                startStep = int32(startTime / simulationDt) + 1;
                endStep = int32(endTime / simulationDt) + 1;
                curDerivative = 0;
                for curStep = startStep + 1 : endStep
                    curTime = double(curStep - 1) * simulationDt;
                    curState = xList(:, curStep);
                    curBeta = betaList(curStep);
                    curN = n1(curStage);
                    integrand = pList(:, curStep)' * getDfDn(curState, curBeta, curTime, curN) ...
                            + getDLDn(curState, curBeta, curTime, curN);
                    curDerivative = curDerivative + integrand * simulationDt;
                end
                dJDn1(curStage) = curDerivative;
                startTime = endTime;
            end
            
                
            % calculate the armijo stepsize for current input
            flag = 0;
            while flag == 0
                nextN1 = n1 - dJDn1 * inputArmijoStepsize;
                nextN1 = real(nextN1);
                nextN1 = max(nextN1, [minN, minN, minN]);
                nextN1 = min(nextN1, [maxN, maxN, maxN]);
                nextXList = xList;
                nextNList = nList;
                nextBetaList = betaList;
        
        
                stage = 1;
                for curStep = 1 : finalStep + 1
                    % The input to the system should be adjusted according to
                    % endogenous response.
                    curTime = (curStep - 1) * simulationDt;
                    if stage <= length(t1) && curTime >= t1(stage)
                        stage = stage + 1;
                    end
                    
                    curState = nextXList(:, curStep);
                    
                    % calculate endogenous response
                    curKappa = getTVKappa(curTime);
                    R = curState(4);
                    curEndo = curKappa * delta * theta * R;
                    maxNWithEndo = maxN - curEndo;
                    
                    % adjust the input according to endogenous response
                    curN = nextN1(stage);
                    if curN > maxNWithEndo
                        curN = maxNWithEndo;
                    end
                    if curN < minN
                        curN = minN;
                    end
                    
                    nextNList(curStep) = curN;
                    curBeta = getBetaFromN(curN, curTime);
                    nextBetaList(curStep) = curBeta;
                    if curStep <= finalStep
                        curDynamics = seirdDynamics(curState, curBeta);
                        nextState = curState + curDynamics * simulationDt;
                        nextXList(:, curStep + 1) = nextState;
                    end
                end
                nextCost = costFunctionIntegral(nextXList, nextBetaList, simulationDt);
                nextCost = real(nextCost);
                if real(curCost - nextCost) >= inputArmijoAlpha * inputArmijoStepsize * (dJDn1 * dJDn1')
                    flag = 1;
                    disp(inputArmijoStepsize)
                else
                    inputArmijoStepsize = inputArmijoStepsize * inputArmijoBeta;
                    disp(inputArmijoStepsize)
                    if inputArmijoStepsize < 1e-8
                        flag = 1;
                    end
                end
            end
            n1 = n1 - inputArmijoStepsize * dJDn1;
            disp(dJDn1)
            n1 = max(minN, n1);
            n1 = min(0.99, n1);
        else
            % calculate the armijo stepsize for current derivative
                 %t1 = t1 + 0.01;
            % The gradient becomes k-dimensional, where k is the number of
            % switches.
            t1Step = ceil(t1 * 100) + 1;
            %curX = xList(:, t1Step);
            totalSwitch = length(t1);
            dt1 = zeros(1, totalSwitch);
            for switchIndex = 1 : length(t1)
                curStep = t1Step(switchIndex);
                curX = xList(:, curStep);
                curT1 = t1(switchIndex);
                % calculate the derivative to time
                % note that the input before time 0 is invalid
                % try to make the minimal value 0.01 and maximal value tf - 0.01.
                % How to find the corresponding two inputs???
                nBefore = n1(switchIndex);
                nAfter = n1(switchIndex + 1);
                
                %betaBefore = betaList(curStep - 1);
                %betaAfter = betaList(curStep);
                betaBefore = getBetaFromN(nBefore, curT1);
                betaAfter = getBetaFromN(nAfter, curT1);
                curDt1 = pList(:, curStep)' * (seirdDynamics(curX, betaBefore) - seirdDynamics(curX, betaAfter) );
             
        
                dt1(switchIndex) = pList(:, curStep)' * (seirdDynamics(curX, betaBefore) - seirdDynamics(curX, betaAfter) )+ ...
                    costFunction(curX, betaBefore, curT1) - costFunction(curX, betaAfter, curT1);
        
            end
            if (t1(1) == t1(2)) && (dt1(1) < 0) && (dt1(2) > 0)
                disp("t1 = t2")
                zigzagCount = 0;
                zigzagState = 1 - zigzagState;
                continue;
            end
            
            flag = 0;
            stepSize = 160;
            
            while flag == 0
                nextT1 = t1 - dt1 * stepSize;
                nextT1 = real(nextT1);
                nextT1 = min(nextT1, [finalTime, finalTime]);
                nextT1 = max(nextT1, [0, 0]);
                if nextT1(2) < nextT1(1)
                    nextT1(2) = nextT1(1);
                end
                nextXList = xList;
                nextNList = nList;
        
                stage = 1;
                for curStep = 1 : finalStep + 1
                    curTime = (curStep - 1) * simulationDt;
                    if stage <= length(nextT1) && curTime >= nextT1(stage)
                        stage = stage + 1;
                    end
                    
                    curState = nextXList(:, curStep);
                    
                    % calculate endogenous response
                    curKappa = getTVKappa(curTime);
                    R = curState(4);
                    curEndo = curKappa * delta * theta * R;
                    maxNWithEndo = maxN - curEndo;
                    
                    % adjust the input according to endogenous response
                    curN = n1(stage);
                    if curN > maxNWithEndo
                        curN = maxNWithEndo;
                    end
                    if curN < minN
                        curN = minN;
                    end
                    
                    nextNList(curStep) = curN;
                    curBeta = getBetaFromN(curN, curTime);
                    nextBetaList(curStep) = curBeta;
                    if curStep <= finalStep
                        curDynamics = seirdDynamics(curState, curBeta);
                        nextState = curState + curDynamics * simulationDt;
                        nextXList(:, curStep + 1) = nextState;
                    end
                end
                nextCost = costFunctionIntegral(nextXList, nextBetaList, simulationDt);
                nextCost = real(nextCost);
                if real(curCost - nextCost) >= armijoAlpha * stepSize * (dt1 * dt1')
                    flag = 1;
                else
                    stepSize = stepSize * armijoBeta;
                    if stepSize < 1e-2
                        flag = 1;
                    end
                end
        
            end
        
            nextT1 = t1 - dt1 * stepSize;
            nextT1 = real(nextT1);
            nextT1 = min(nextT1, [finalTime - 0.01, finalTime - 0.01]);
            nextT1 = max(nextT1, [0.01, 0.01]);
            if nextT1(2) < nextT1(1)
                nextT1(2) = nextT1(1);
            end
            
            t1 = nextT1;
            disp(stepSize);
            stepsizeList(curIteration) = stepSize;
        end
        
        % t1 = t1 - dt1 * stepSize;
        % t1 = real(t1);
        % t1 = min(t1, [finalTime, finalTime]);
        % t1 = max(t1, [0, 0]);
        % if t1(2) < t1(1)
        %     t1(2) = t1(1);
        % end
    end
    saveFileName = "./test_results/beta_W_" + newParam +'.mat';
    save(saveFileName);