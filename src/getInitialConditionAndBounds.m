function [initQ,lBound,uBound] = getInitialConditionAndBounds(N_diverse,interpSteps)
    %Creates an initial condition for the wavetable synth
    %initQ contains initial values for each parameter
    % lBound and uBound contain the domain boundaries for each 

    %Bounds for all envelopes
    envInit = [0.0001 , 0.25 , 1, 0.25, 1.0, 0.5];
    envLow = [0.0001 , 0.0001 , 0, 0.0001, 0, 0];
    envUp = [2.0 , 2.0 , 1.0 , 2.0, 1.0, 1.0];

    %VOLUMES
    initQ(1:2) = [0.5, 0.5];

    %ENVELOPES
    initQ(3:8) = envInit;
    initQ(9:14) = envInit;
    initQ(15:20) = envInit;
    initQ(21:26) = envInit;

    %FILTER
    initQ(27) = 20000;
    initQ(28) = 20000;
    initQ(29) = 0;
    initQ(30) = 0;

    initQ(N_diverse:N_diverse+interpSteps-1) = 2 * (0.5 - rand(interpSteps,1)); %sin(2*pi*(0:1/interpSteps:1-1/interpSteps));
    initQ(N_diverse+interpSteps) = 1.0;
    initQ(N_diverse+interpSteps+1:N_diverse+(2*interpSteps)) = 2 * (0.5 - rand(interpSteps,1)); %sin(2*pi*(0:1/interpSteps:1-1/9));
    initQ(N_diverse+(2*interpSteps)+1) = 0.0;

    lBound(1:2) = [0.0, 0.0];
    lBound(3:8) = envLow;
    lBound(9:14) = envLow;
    lBound(15:20) = envLow;
    lBound(21:26) = envLow;
    lBound(27) = 50;
    lBound(28) = 50;
    lBound(29) = 0;
    lBound(30) = 0;

    lBound(N_diverse:N_diverse+interpSteps-1) = -1;
    lBound(N_diverse+interpSteps) = 0.0;
    lBound(N_diverse+interpSteps+1:N_diverse+(2*interpSteps)) = -1;
    lBound(N_diverse+(2*interpSteps)+1) = 0.0;

    uBound(1:2) = [1.0, 1.0];
    uBound(3:8) = envUp;
    uBound(9:14) = envUp;
    uBound(15:20) = envUp;
    uBound(21:26) = envUp;
    uBound(27) = 20000;
    uBound(28) = 20000;
    lBound(29) = 1.0;
    lBound(30) = 1.0;

    uBound(N_diverse:N_diverse+interpSteps-1) = 1;
    uBound(N_diverse+interpSteps) = 1.0;
    uBound(N_diverse+interpSteps+1:N_diverse+(2*interpSteps)) = 1;
    uBound(N_diverse+(2*interpSteps)+1) = 1.0;

end
