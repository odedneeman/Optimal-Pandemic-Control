function ret = getDfDx(x, beta)
    global sigma;
    global gamma;
    global theta;
    global delta;
    
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    C = x(6);

    ret = [-beta * I, 0, -beta * S, 0, 0, 0;
            beta * I, -sigma, beta * S, 0,0,0;
            0, sigma, -gamma, 0, 0, 0;
            0, 0, gamma, -theta, 0, 0;
            0, 0, 0, delta * theta, 0, 0;
            0, 0, 0, (1 - delta) * theta, 0, 0];

end