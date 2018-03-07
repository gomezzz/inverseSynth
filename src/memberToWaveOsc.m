function [wv] = memberToWaveOsc(member,tbSize,Fs,lengthInS, N_diverse, interpSteps)
    %Creates a wavetable Synth using the passed parameters as inputs for all parameters in the synth and necessary envelopes
    vol1 = member(1);
    vol2 = member(2);

    env1Params = member(3:8); %ADSR
    env2Params = member(9:14);
    phaseEnv1Params = member(15:20);
    phaseEnv2Params = member(21:26);

    ftr1Cut = member(27);
    ftr2Cut = member(28);
    ftr1Res = member(29);
    ftr2Res = member(30);

    env1 = envelopeGenerator(Fs,env1Params(1),env1Params(2),env1Params(3),env1Params(4),env1Params(5), env1Params(6));
    env2 = envelopeGenerator(Fs,env2Params(1),env2Params(2),env2Params(3),env2Params(4),env2Params(5), env2Params(6));
    phaseEnv1 = envelopeGenerator(Fs,phaseEnv1Params(1),phaseEnv1Params(2),phaseEnv1Params(3),phaseEnv1Params(4),phaseEnv1Params(5), phaseEnv1Params(6));
    phaseEnv2 = envelopeGenerator(Fs,phaseEnv2Params(1),phaseEnv2Params(2),phaseEnv2Params(3),phaseEnv2Params(4),phaseEnv2Params(5), phaseEnv2Params(6));

    tb1 = member(N_diverse:N_diverse+interpSteps-1);
    smoothness1 = member(N_diverse+interpSteps);
    tb2 = member(N_diverse+interpSteps+1:N_diverse+2*interpSteps);
    smoothness2 = member(N_diverse+(2*interpSteps)+1);

    wv = wavetableOscillator(env1,env2,tbSize,Fs,tb1,tb2,vol1,vol2,phaseEnv1,phaseEnv2,ftr1Cut,ftr1Res,ftr2Cut,ftr2Res,interpSteps);

    wv = wv.setTableUsingSpline(tb1,tb2,smoothness1,smoothness2);
end
