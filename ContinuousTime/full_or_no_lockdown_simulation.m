clear; 
addpath("..\Common\")
% Define global constants    

    global fullLockdown;
    fullLockdown = 0;

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
    
    % Parameter set 2: beta function (transmission rate)
    global betaW;
    global betaN;
    global betaLambda;
    global lambda;
    global alpha;
    
    betaW = 0.3;
    % if newParams.isKey("betaW")
    %     betaW = newParams("betaW");
    % end
    betaN = 0.53;
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
    fullLockdownN = 0.75;
    
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
    timeList = 0 : simulationDt : finalTime;
    fixedStepsize = 0.01;
    maxIter = 500;
    
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
    nList = ones(inputDimension, finalStep) * maxN;
    if fullLockdown == 1
        nList = ones(inputDimension, finalStep) * fullLockdownN;
    end
    fminTimer = 0.0;
    nrTimer = 0.0;
   
    %%% continue optimization from previous results
    %load(".\ContinuousSensitivityResultsWithEndo\beta_N_0.53.mat", "betaList")
    % if loadFileName ~= ""
    %     load(loadFileName, "nList")
    %     disp("Successfully loaded initial guess from " + loadFileName)
    %end
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

    lock_str = "No_Intervention";
    if fullLockdown == 1
        lock_str = "Full_Lockdown";
    end
    figure(1);
    plot(timeList, xList(5,:)* 100000)
    xlabel("Time (days)")
    ylabel("Death toll per 100,000")
    title("Death toll over time with " + lock_str);

    figure(2);
    Re = (betaList * 4) .*  xList(1, :);
    plot(timeList, Re)
    plot(timeList, xList(1,:) .* betaList / gamma);
    xlabel("Time (days)")
    ylabel("Effective R value")
    title("Effective R over time with " + lock_str);

    saveName = 'ContinuousSensitivityResultsWithEndo\betaW_0.3_' + lock_str + '.mat';
    save(saveName)
    % tV=540;
    % [betaListVac, xListVac, nListVac] = getVaccineState(xList, betaList, nList, tV, simulationDt, finalStep);
    % allList = (vertcat(timeList, xListVac, betaListVac, nListVac))';
    % runNames = {'Step', 'S', 'E', 'I', 'R', 'D', 'C', 'Beta', 'Employment'};
    % dataTable = array2table(allList, 'VariableNames', runNames);
    % excelFileName = "" + lock_str + ".xlsx";
    % writetable(dataTable, excelFileName, 'Sheet', 'Sheet1');


