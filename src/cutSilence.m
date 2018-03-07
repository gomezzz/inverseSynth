function sig = cutSilence(signal)
    %Cuts intro silence from a signal
    cutoffIndex = length(signal);
    for i=length(signal):-1:1
        if abs(signal(i)) > 0.0001
            break;
        else
            cutoffIndex = i;
        end
    end

    sig = signal(1:cutoffIndex);
end
