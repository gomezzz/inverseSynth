function val = objFunc(target,candidate,fft_target)
        %Computes the objective function for a passed solution candidate, several functions are implemented below
        type = 'transient';
        if strcmp(type,'time')
            val = timeErr(target,candidate);
        elseif strcmp(type,'freq')
            val = freqErr(target,candidate,fft(target));
        elseif strcmp(type,'combined')
            val = combinedFreqTimeErr(target,candidate,fft(target));
        elseif strcmp(type,'transient')
            val = transientAwareTimeErr(target,candidate);
        end
end

function val = transientAwareTimeErr(target,candidate)
    %Computes an error in the time domain with a weighted treatmeant of the transient phase

    N = length(target);

    %Transient Error
    transientLengthInSamples = min(N,250);
    tmp = target(1:transientLengthInSamples)-candidate(1:transientLengthInSamples);
    transientError = sum(tmp .* tmp);
    transientNorm = sum(target(1:transientLengthInSamples) .* target(1:transientLengthInSamples));
    tmp = target(transientLengthInSamples:end)-candidate(transientLengthInSamples:end);

    %After transient Error
    restError = sum(tmp .* tmp);
    restNorm = sum(target(transientLengthInSamples:end) .* target(transientLengthInSamples:end));

    if restError > 0 && transientNorm > 0 && mean(abs(candidate)) > 0.05
        val = 0.5 *  (transientError / transientNorm) + 0.5 * (restError / restNorm);
    else
        val = 10000000;
    end
end

function val = timeErr(target,candidate)
    %Simple average absolute error in the time domain

    tmp = 10000;
    N = length(target);
    scaling = 0.5 + (N / tmp : -1 / tmp : 0 + 1 / tmp) .* (N / tmp : -1 / tmp : 0 + 1 / tmp);
    err = abs((target - candidate)) .* scaling;
    val = sum(err) / N;
end

function val = freqErr(target,candidate,t_fft)
     %Simple average absolute error in the frequency domain

    c_fft = fft(candidate);

    t_fft = t_fft(1:floor(length(target)/2)+1);
    c_fft = c_fft(1:floor(length(target)/2)+1);

    t_fft(1) = 0;
    c_fft(1) = 0;

    t_fft = t_fft ./ max(t_fft);
    c_fft = c_fft ./ max(c_fft);

    err = abs(t_fft - c_fft);
    val = sum(err);
end

function val = combinedFreqTimeErr(target,candidate,t_fft)
    %Combination of average relative error in the frequency and time domain domain

    %Time domain error
    err = abs(target - candidate);
    val = sum(err) / sum(abs(target));
  
    val = val + abs(abs(max(target)-min(target))-abs(max(candidate)-min(candidate)));


    %%FFT Computations
    c_fft = fft(candidate);

    t_fft = t_fft(1:floor(length(target)/2)+1);
    c_fft = c_fft(1:floor(length(target)/2)+1);

    t_fft(1) = 0;
    c_fft(1) = 0;

    t_fft = t_fft ./ max(t_fft);
    c_fft = c_fft ./ max(c_fft);

    %Frequency domain error
    err = abs(t_fft - c_fft);
    val = val + sum(err) / sum(abs(t_fft));
    
end
