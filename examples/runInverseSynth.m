clear all
close all

%Read a file
target = audioread('booka.wav');
target(:,1) = 0.5 * (target(:,1) + target(:,2)); %make mono
target(:,2) = [];
target = target';

%Cut the silence in the beginning
target = cutSilence(target);

%Setup the wavetable synth parameters
tbSize = 256;   %wavetable size
f=110.0;    %Target frequency (currently only one note at specified frequency is possible. Pitch estimation is not yet part of the project)
Fs = 44100;     %Sampling rate
N_diverse = 30;     %number of parameters
interpSteps = 15;   %Interpolation steps in the wavetable
evaluations = 20000; %how many optimization steps? You will need at least a few thousands

[snd,err,mem]= DE_wv(target,evaluations,Fs,tbSize,f,N_diverse,interpSteps);         %Use differential evolution
% [snd,err,mem]= PSO_wv(target,evaluations,Fs,tbSize,f,N_diverse,interpSteps);    %Use particle swarm optimization

signal(:,1) = snd;  %Make result stereo
signal(:,2) = snd;

audiowrite('result.wav',signal,44100); %write result to wav

%Plot the result!

wv = memberToWaveOsc(mem,tbSize,44100,1.0,N_diverse,interpSteps);

plotWave(snd,44100,target)
plotSpectrum(snd,44100,target)
plotRMS(snd,44100,target)
plotEnvelope(wv.wvTb1,'Wavetable 1',wv.wvTb2,'Wavetable 2');
plotEnvelope(wv.wv1Env.env,'Envelope 1',wv.wv2Env.env,'Envelope 2');
plotEnvelope(wv.phaseEnv1.env,'Phase Envelope 1',wv.phaseEnv2.env,'Phase Envelope 2');


%Play around and play some notes
scale = f / 440.0;
f1 = scale * 440;
f2 = scale * 659.25;
f3 = scale * 880;
f4 = scale * 1046.50;
f5 = scale * 1318.51;

ll = 0.25;
lInSamples = round(ll*Fs);

resultsPlay(1:lInSamples) = wv.getSound(f1,ll);
resultsPlay(lInSamples+1:2*lInSamples) =  wv.getSound(f2,ll);
resultsPlay(2*lInSamples+1:3*lInSamples) =  wv.getSound(f3,ll);
resultsPlay(3*lInSamples+1:4*lInSamples) =  wv.getSound(f4,ll);
resultsPlay(4*lInSamples+1:5*lInSamples) =  wv.getSound(f5,ll);
resultsPlay(5*lInSamples+1:6*lInSamples) =  wv.getSound(f4,ll);
resultsPlay(6*lInSamples+1:7*lInSamples) =  wv.getSound(f3,ll);
resultsPlay(7*lInSamples+1:8*lInSamples) =  wv.getSound(f2,ll);
resultsPlay(8*lInSamples+1:9*lInSamples) =  wv.getSound(f1,ll);

stereoSig(1,:) = resultsPlay;
stereoSig(2,:) = resultsPlay;
audiowrite('resultPlay.wav',resultsPlay,44100);
