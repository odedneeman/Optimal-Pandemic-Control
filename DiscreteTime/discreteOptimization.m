function ret = discreteOptimization(newParams, saveFileName)
    % Define global constants 
    % Parameter set 1: model parameters
    global sigma;
    global gamma;
    global theta;
    global delta;
    sigma = 1/3;
    theta = 1/11;
    delta = 0.008;
    if newParams.isKey("delta")
        delta = newParams("delta");
    end
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
    kappa = 30500;
    muF = 245;
    sigmaF = 27.5;
    phiKappa = 0.18;
    
    
    global maxNRIter;
    maxNRIter = 10;
    
    
    x0 = [0.9999; 0.00005; 0.00005; 0; 0; 0];
    simulationDt = 0.01;
    finalTime = 635;
    % changing finalTime to match confident bound for gumble dist.
    % if newParams.isKey("muTv")
    %     confidentialEqn = @(x)(getGumbelCDF(x) - 0.99);
    %     confidentialBound = fsolve(confidentialEqn, 500);
    %     finalTime = ceil(confidentialBound);
    % end
    finalStep = round(finalTime / simulationDt);
    timeList = 0 : simulationDt : finalTime;
    fixedStepsize = 0.01;
    maxIter = 90; %def 90
    
    xList = zeros(6, finalStep + 1); 
    xList(:,1) = x0;
    % List of costates
    pList = zeros(6, finalStep + 1);
    
    betaList = ones(1, finalStep + 1);
    nList = ones(1, finalStep + 1);
    
    
    %t1 = [200, 450]; %was def
    %n1 = [0.6, 0.7, 0.8]; %was def
    %t1 = [300, 400]; %def for benchmark
    %t1 = [200, 400]; %new def
    t1 = [150, 350];
    %t1 = [150, 400];
    %t1 = [300, 450];
    %n1 = [0.9, 0.7, 0.9]; %new def
    n1 = [0.8, 0.8, 0.9]; %for betaN=0.4, 0.6, tV=500,520,600
    %n1 = [0.871167022022223, 0.872933608119833, 0.910999137883876];
    %n1 = [0.8, 0.75, 0.85]; %for betaW=0.45
    %n1 = [0.75, 0.75, 0.85];

    initialN1 = n1;
    initialT1 = t1;
    costList = zeros(1, maxIter);
    
    switchingArmijoBeta = 0.5;
    switchingArmijoAlpha = 0.1;
  
    inputArmijoBeta = 0.5;
    inputArmijoStepsize = 0.001;
    inputArmijoAlpha = 0.01;
    
    zigzagCount = 0;
    zigzagState = 0;
    for curIteration = 1 : maxIter
        disp(n1)
        disp(t1)
        zigzagCount = zigzagCount + 1;
        if zigzagCount == 30 %def = 30
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
            if curEndo > 1 - minN
                curEndo = 1 - minN;
            end
            maxNWithEndo = maxN - curEndo;
            % adjust the input according to endogenous response
            curN = n1(stage);
            if curN > maxNWithEndo
                curN = maxNWithEndo;
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
            curBeta = betaList(:, curStep); 
    
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
                    if curEndo > 1 - minN
                        curEndo = 1 - minN;
                    end
                    maxNWithEndo = maxN - curEndo;
         
                    % adjust the input according to endogenous response
                    curN = nextN1(stage);
                    if curN > maxNWithEndo
                        curN = maxNWithEndo;
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
            t1Step = ceil(t1 * 100) + 1;
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
                    if curEndo > 1 - minN
                        curEndo = 1 - minN;
                    end
                    maxNWithEndo = maxN - curEndo;
                    % adjust the input according to endogenous response
                    curN = n1(stage);
                    if curN > maxNWithEndo
                        curN = maxNWithEndo;
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
                if real(curCost - nextCost) >= switchingArmijoAlpha * stepSize * (dt1 * dt1')
                    flag = 1;
                else
                    stepSize = stepSize * switchingArmijoBeta;
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
        end
        
    end
    if saveFileName ~= ""
        save(saveFileName);
    end
end

