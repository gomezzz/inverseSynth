function [snd,err,mem] = PSO_wv(target,max_evals,Fs,tbSize,f,N_diverse,interpSteps)
    %An implementation  of a Particle Swarm Optimization, SPSO 2011 is used

    fevals = 0;

    if (size(target,2) == 2)
      target(:,2) = [];
      target = target';
    end

    N = length(target);
    NinS = length(target) / Fs;

    fft_target = fft(target);

    [initQ,lBound,uBound] = getInitialConditionAndBounds(N_diverse,interpSteps);

    NP = 100;
    dim = length(lBound);
    K = 3;

    omega = 0.7;
    c1 = 1.2;
    c2 = 1.2;

    Xi = zeros(NP,dim);
    Vi = zeros(NP,dim);
    pi = zeros(NP,dim);
    li = zeros(NP,dim);
    fi = zeros(1,NP);
    fi_cur = zeros(1,NP);
    fNeighbori = zeros(1,NP);
    neighborhood = zeros(NP,K);

    bestmem = initQ;
    [bestmem,bestval,fevals] = evaluate(lBound,uBound,bestmem,fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);

    max_withoutchange = 10;
    it_withoutchange = 0;
    oldbest = bestval;

    %initialization
    for i=1:NP
    Xi(i,:) = lBound + rand(1,dim) .* (uBound-lBound);

    Vi(i,:) = 0.5 * (rand(1,dim) .* (uBound-lBound) - Xi(i,:));

    neighborhood(i,:) = round(1+(NP-1)*rand(1,K));

    [Xi(i,:),fi(i),fevals] = evaluate(lBound,uBound,squeeze(Xi(i,:)),fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);

    if fi(i) < bestval
        bestmem = Xi(i,:);
        bestval = fi(i);
        disp(['New best at        ' num2str(bestmem)]);
    end
    end

    for i=1:NP
    pi(i,:) = Xi(i,:);
    fNeighbori(i) = fi(i);
    for j=1:K
      neigh = neighborhood(i,j);
      if (fNeighbori(i) > fi(neigh))
        li(i,:) = Xi(neigh,:);
        fNeighbori(i) = fi(neigh);
      end
    end
    end


    iteration = 0;
    %main loop
    while (fevals < max_evals)
        iteration = iteration + 1;
        disp(['Iteration = ' num2str(iteration) '   best = ' num2str(bestval) '   fevals = ' num2str(fevals) '    av = '  num2str(sum(fi(:))/NP)]);

        %update positions and velocities
        for i=1:NP
          %get new velocities
          Vi(i,:) = updateVelocity(lBound,uBound,Vi(i,:),Xi(i,:),pi(i,:),li(i,:),dim,c1,c2,omega);

          %get new positions
          Xi(i,:) = updatePosition(Xi(i,:),Vi(i,:));
        end


        %evaluate new positions and update personal best
        for i=1:NP
          [Xi(i,:),fi_cur(i),fevals] = evaluate(lBound,uBound,squeeze(Xi(i,:)),fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);

          if (fi(i) > fi_cur(i))
              pi(i,:) = Xi(i,:);
              fi(i) = fi_cur(i);

              if fi(i) < bestval
                  bestmem = Xi(i,:);
                  bestval = fi(i);
                %   disp(['New best at        ' num2str(bestmem)]);
              end
          end
        end

        %update local best
        for i=1:NP
          for j=1:K
            neigh = neighborhood(i,j);
            if (fNeighbori(i) > fi(neigh))
              li(i,:) = Xi(neigh,:);
              fNeighbori(i) = fi(neigh);
            end
          end
        end

        if oldbest == bestval
          it_withoutchange = it_withoutchange + 1;
        else
          it_withoutchange = 0;
        end

        oldbest = bestval;

        if it_withoutchange == max_withoutchange
          disp('reinit neighborhood topology');
          for i=1:NP
            neighborhood(i,:) = round(1+(NP-1)*rand(1,K));
          end
          it_withoutchange = 0;
        end
    end

    disp(['Done. Best Val was ' num2str(bestval) ' at ' num2str(bestmem)]);

    mem = bestmem;
    err = bestval;
    snd = memberToWav(bestmem,tbSize,Fs,NinS,f,N_diverse,interpSteps);
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

function [Q,val,fevals] = evaluate(lBound,uBound,Q,fevals,target,N,tbSize,Fs,NinS,f,fft_target,N_diverse,interpSteps);
  Q = inDomain(lBound,uBound,Q);
  snd = memberToWav(Q,tbSize,Fs,NinS,f,N_diverse,interpSteps);
  val = objFunc(target,snd,fft_target);
  fevals = fevals+1;
end

function Vi = updateVelocity(lBound,uBound,Vi,Xi,pi,li,dim,c1,c2,omega)
  p = Xi + (c1 * rand(1,dim)) .* (pi - Xi);
  l = Xi + (c2 * rand(1,dim)) .* (li - Xi);

  G = 1/3 * (Xi(:) + p(:) + l(:));

  zeropointfive(1:dim) = 0.5;

  r = rand(1,dim) - zeropointfive;

  x = G' + (2 * r) .* norm(G'-Xi);

  Vi = omega * Vi + x - Xi;
end

function Xi = updatePosition(Xi,Vi)
  Xi = Xi + Vi;
end
