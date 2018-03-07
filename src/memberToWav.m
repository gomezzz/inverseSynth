function snd = memberToWav(member,tbSize,Fs,lengthInS,f,N_diverse,interpSteps)
    %Converts an optimization population member to a 
    wv = memberToWaveOsc(member,tbSize,Fs,lengthInS,N_diverse,interpSteps);

    %get the resulting sound
    snd = wv.getSound(f,lengthInS); 
    
end
