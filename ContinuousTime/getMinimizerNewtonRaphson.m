function ret = getMinimizerNewtonRaphson(x, p, curTime, minBeta, maxBeta)
    global maxNRIter;
    curMinimizer = minBeta + 0.01;
    for curIter = 1:maxNRIter
        derivative = hamiltonianDerivative(x, curMinimizer, p, curTime);
        if abs(real(derivative))<0.001
            break
        end
        doubleDerivative = hamiltonianDoubleDerivative(x, curMinimizer, p, curTime);
        direction = derivative / doubleDerivative;
        
        curMinimizer = curMinimizer - direction;
        if curMinimizer < minBeta
            curMinimizer = minBeta;
            break;
        else 
            if curMinimizer > maxBeta
                curMinimizer = maxBeta;
                break;
            end
        end
    end
    ret = curMinimizer;
    curMinimum = hamiltonian(x, ret, p, curTime);
    boundaryMinimum = hamiltonian(x, maxBeta, p, curTime);
    if boundaryMinimum < curMinimum
        ret = maxBeta;
        %disp("***minimum at boundary***")
    end
end