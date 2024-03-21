function ret = getArmijoStepsize(x, beta, p, minimizerBeta, simulationDt, finalTime, amijoAlpha, amijoBeta, ...
            lastStepsize)
    finalStep = int32(finalTime/simulationDt);
    if lastStepsize < 0.25
        ret = lastStepsize * 4;
    else 
        if lastStepsize < 0.5
            ret = lastStepsize * 2;
        else 
            ret = 1;
        end
    end
    optimal = 0;
    while optimal == 0
        curCost = costFunctionIntegral(x, beta, simulationDt);
        nextBeta = beta + ret *  (minimizerBeta - beta);
        nextX = x;
        for curStep = 2:finalStep
            dx = seirdDynamics(nextX(:, curStep - 1), nextBeta(:, curStep - 1));
            nextX(:, curStep) = nextX(:, curStep - 1) + dx * simulationDt;
        end
        
        nextCost = costFunctionIntegral(nextX, nextBeta, simulationDt);
        gradNorm = 0;
        for curStep = 1:finalStep
            curTime = double(curStep - 1) * simulationDt;
           
            gradNorm = gradNorm + (hamiltonian(x(:, curStep), minimizerBeta(:, curStep), p(:, curStep), curTime) ...
                - hamiltonian(x(:, curStep), beta(:, curStep), p(:, curStep), curTime)) * simulationDt;
                
        end
        if real(nextCost - curCost) < real(amijoAlpha * ret * gradNorm)
            optimal = 1;
        else
            ret = ret * amijoBeta;
            disp(ret)
            if ret < 1e-10
                optimal = 1;
            end
        end
    end
end