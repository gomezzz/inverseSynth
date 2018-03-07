# inverseSynth
This project implements a wavetable Synthesizer that uses Evolutionary algorithms to determine its settings to replicate a given sound.
This was created as a proof of concept.

This project relies on some of the funcitons I implemented in [MDSP](https://github.com/gomeZZZ/MDSP), which is also on GitHub.

The script in the examples folder showcases usage. 

Differential Evolution and Particle Swarm Optimization are implemented and used as Evolutionary Algorithms. 
The code is still somewhat limited, a target frequency has to be given, it's only mono and not ready to use on "real" sounds. 
The examples folder contains some simple sounds, which can be replicated quite successfully. 

# Example - Recreating the bass sound in Booka Shade's - In White Rooms

The blue "filtered" signal here is the original sound, the red "signal" is the sound created by the synth.

![alt text](https://github.com/gomeZZZ/inverseSynth/blob/master/examples/ExampleResult3.png)
![alt text](https://github.com/gomeZZZ/inverseSynth/blob/master/examples/ExampleResult2.png)
![alt text](https://github.com/gomeZZZ/inverseSynth/blob/master/examples/ExampleResult1.png)
