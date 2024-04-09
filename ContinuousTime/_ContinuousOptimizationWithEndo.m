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
finalStep = round(finalTime / simulationDt) + 1;
timeList = 0 : simulationDt : finalTime;
fixedStepsize = 0.01;
maxIter = 200;

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
nList = ones(inputDimension, finalStep) * minN;
fminTimer = 0.0;
nrTimer = 0.0;

for curStep = 1 : finalStep
    curTime = (curStep - 1) * simulationDt;
    curBeta = getBetaFromN(minN, curTime);
    betaList(curStep) = curBeta;
end
costList = zeros(1, maxIter);
costList(1) = costFunctionIntegral(xList, betaList, simulationDt);
armijoStepsize = 0.5;
simulationBeginTime = tic;
for curIteration = 1:maxIter
    disp(curIteration)
    % Simulate the system to obtain the system states
    for curStep = 2 : finalStep
        curTime = (curStep - 2) * simulationDt;
        
        curX = xList(:, curStep - 1);
        curBeta = betaList(curStep - 1);
        dX = seirdDynamics(curX, curBeta);
        
        nextX = curX + dX * simulationDt;
        xList(:, curStep) = nextX;
    end
    

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
    for curStep = 1 : finalStep 
        curTime = (curStep - 1) * simulationDt;
        curP = costateList(:, curStep);
        curX = xList(:, curStep);
        betaRange = getBetaRange(curTime, curX);
    
        target = @(beta)(hamiltonian(curX, beta, curP, curTime));
        minBeta = betaRange(1);
        maxBeta = betaRange(2);

        %minimizerBeginTime = tic;
        %curMinimizer = fminbnd(target, minBeta, maxBeta);
        %minimizerRunTime = toc(minimizerBeginTime);
        %fminTimer = fminTimer + minimizerRunTime;
        
        minimizerBeginTime = tic;
        curMinimizer = getMinimizerNewtonRaphson(curX, curP, curTime, minBeta, maxBeta);
        minimizerRunTime = toc(minimizerBeginTime);
        nrTimer = nrTimer + minimizerRunTime;

        minimizerList(curStep) = curMinimizer;
    end
    % Uncomment this block to use Armijo stepsize
    %armijoStepsize = getArmijoStepsize(xList, betaList, costateList, minimizerList, simulationDt, finalTime, armijoAlpha, armijoBeta, ...
    %        armijoStepsize);
    %betaList = betaList + armijoStepsize * (minimizerList - betaList);
    
    betaList = betaList + fixedStepsize * (minimizerList - betaList);


    for curStep = 1 : finalStep
        curTime = (curStep - 1) * simulationDt;
        nList(curStep) = getNFromBeta(real(betaList(curStep)), curTime);
    end
    % Calculate the cost function with current input and state
    curCost = costFunctionIntegral(xList, betaList, simulationDt);
    disp(curCost)
    if curCost > costList(curIteration)
        break;
    end
    costList(curIteration + 1) = curCost;

end


%simulationTotalTime = toc(simulationBeginTime);
%saveFileName = "./sensitivity_results_continuous/beta_W_" + param +'.mat';
%save(saveFileName);