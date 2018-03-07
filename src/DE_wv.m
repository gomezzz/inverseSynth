function [snd,err,mem] = DE_wv(target,max_evals,Fs,tbSize,f,N_diverse,interpSteps)
    %An implementation  of a Differential Evolution Optimization, based on the original paper by Storn and Price

    fevals = 0;

    if (size(target,2) == 2)
      target(:,2) = [];
      target = target';
    end

    popsize = 100;

    CR = 0.8;
    F = 0.8;

    N = length(target);
    NinS = length(target) / Fs;

    fft_target = fft(target);

    [initQ,lBound,uBound] = getInitialConditionAndBounds(N_diverse,interpSteps);

    if nargin < 2
        max_evals = 5000;
    end

    nochangeiter = 10;

    noChangeCounter = 0;

    pop = zeros(popsize,length(initQ));
    val = zeros(popsize,1);
    tempval = zeros(popsize,1);

    iteration = 0;
    max_iterations = 1000;

    reinit = 0;
    pop(1,:) = initQ;

    snd = memberToWav(pop(1,:),tbSize,Fs,NinS,f,N_diverse,interpSteps);

    val(1) = objFunc(target,snd);

    fevals = fevals + 1;
    bestmem = pop(1,:);
    bestval = val(1);
    istart = 2;

    pop = computeRandomPop(pop,istart:popsize,lBound,uBound);

    for ind = 2:popsize
    %eval
    [pop(ind,:),val(ind),fevals] =  evaluate(lBound,uBound,pop(ind,:),fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);

        if val(ind) < bestval
            bestmem = pop(ind,:);
            bestval = val(ind);
        end
    end

    while iteration < max_iterations && max_evals > fevals
      %REINIT
      if reinit
        disp('Reinitializing...');
        [~, index] = sort(val);

        istart = max(1, round(0.2 * popsize)); % keep best twenty percent of the population members
        index = index(1:istart-1);

        pop(1:istart-1,:) = pop(index,:);
        val(1:istart-1)   = val(index);

        pop = computeRandomPop(pop,istart:popsize,lBound,uBound);

        for ind = istart:popsize
          %eval
          [pop(ind,:),val(ind),fevals] =  evaluate(lBound,uBound,pop(ind,:),fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);

          if val(ind) < bestval
              bestmem = pop(ind,:);
              bestval = val(ind);
              % disp(['New best at        ' num2str(bestmem)]);
          end
        end
        reinit = 0;
      end

      popold = pop;

      %compute new pop
      popnew = computeNewPop(pop,bestmem,CR,F);

      %check survivors
      for ind = 1:popsize
        [popnew(ind,:),tempval(ind),fevals] = evaluate(lBound,uBound,popnew(ind,:),fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);

        if tempval(ind) < val(ind) || (isnan(val(ind)) && ~isnan(tempval(ind)))
          val(ind) = tempval(ind);
          pop(ind,:) = popnew(ind,:);

          if val(ind) < bestval
              bestmem = pop(ind,:);
              bestval = val(ind);
              % disp(['New best at        ' num2str(bestmem)]);
          end
        end
      end

      % check if population has changed
      if isequal(pop, popold)
        noChangeCounter = noChangeCounter + 1;
      else
        noChangeCounter = 0;
      end

      if noChangeCounter >= nochangeiter
         reinit = 1;
      end

      disp(['Iteration = ' num2str(iteration) '   best = ' num2str(bestval) '   fevals = ' num2str(fevals) '    avg = '  num2str(sum(val(:))/popsize)]);


      iteration = iteration +1;
    end

    disp(['Done. Best Val was ' num2str(bestval) ' at ' num2str(bestmem)]);

    mem = bestmem;
    err = bestval;
    snd = memberToWav(bestmem,tbSize,Fs,NinS,f,N_diverse,interpSteps);
end

function popnew = computeRandomPop(pop,indices,lBound,uBound)
    if isempty(indices)
        return
    end

    stddev = 1;

    D = size(pop,2);
    r = rand(1,2);

    for ind = indices
        pop(ind,:) = lBound + stddev * rand(1,D) .* (uBound - lBound);
    end

    popnew = pop;
end

function popnew = computeNewPop(popold,bestmem,CR,F)
    % get dimensions
    [NP, D] = size(popold);

    index = randperm(4);

    rot  = 0:1:NP-1;                 % rotation index array (size NP)
    a1  = randperm(NP);              % shuffle locations of vectors
    rt = rem(rot + index(1), NP);    % rotate indices by index(1) positions
    a2  = a1(rt + 1);                % rotate vector locations

    pm1 = popold(a1,:);              % shuffled population 1
    pm2 = popold(a2,:);              % shuffled population 2

    bm = bestmem(ones(NP, 1), :);    % population filled with the best
                                   % member of the last iteration

    mui = double(rand(NP, D) < CR);  % all random numbers < CR are 1, 0 otherwise

    rotd = 0:1:D-1;                  % rotation index array (size D)

    mui = sort(mui, 2)';           % transpose, collect 1's in each column
    for ind = 1:NP
    n = floor(rand * D);
    if n > 0
      rtd = rem(rotd + n, D);
      mui(:,ind) = mui(rtd + 1,ind); % rotate column ind by n
    end
    end
    mui = mui';                    % transpose back

    mpo = mui < 0.5;                 % inverse mask to mui

    ui = bm  + F * (pm1 - pm2);    % differential variation

    popnew = (popold .* mpo) + (ui .* mui); % crossover
end

function Q = inDomain(lBound,uBound,Q)
    s = length(lBound);
    for ind = 1:s
        if (Q(ind) < lBound(ind))
          Q(ind) = lBound(ind);
        elseif (Q(ind) > uBound(ind))
          Q(ind) = uBound(ind);
        end
    end
end

function [Q,val,fevals] = evaluate(lBound,uBound,Q,fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps)
    Q = inDomain(lBound,uBound,Q);
    snd = memberToWav(Q,tbSize,Fs,NinS,f,N_diverse,interpSteps);
    val = objFunc(target,snd,fft_target);
    fevals = fevals+1;
end
