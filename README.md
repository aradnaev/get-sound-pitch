# get-sound-pitch
extracts pitch from sound waveform
algorithm is based on: 
"[1] IFA Proceedings 17, 1993
ACCURATE SHORT-TERM ANALYSIS OF THE FUNDAMENTAL FREQUENCY AND THE HARMONICS-TO-NOISE RATIO OF A SAMPLED SOUND by 
Paul Boersma"

# installation:
-add AGRSoundPitch.h, AGRSoundPitch.m to your project

# how to use: 
1. initWithBufferSize:(UInt32)bufferSize;

2. processData:(float*)data WithBufferSize:(UInt32)bufferSize WithSamplingRate:(Float64)samplingRate WithGlobalAbsolutePeak:(float*)globalAbsolutePeak WithFFTplotPointer:(float*)FFTdataHolder WithACPlotPointer:(float*)ACPlotHolder
(FFTdataHolder, ACPlotHolder are optional for outputing of FFT and Autocorrelation data)

3. output proprties:
 * previousPitch: previous frame pitch (Hz)
 * currentPitch: current frame pitch (Hz)
 * volume: current volume (dB)

# parameters affect performace
(see [1] for definitions and descriptions):
```Objective-C
    maximumNumberOfCandidatesPerFrame = 2000;
    _voicingThreshold = 0.4; //%VoicingThreshold = 0.4 by default
    _silenceThreshold = 0.1; //%silence threshold = 0.05 by default
    _octaveCost = pow(0.15,2); //default is 0.01 (10% = sqrt(0.01)) // the higher octaveCost, the higher weight of higher frequencies
    _minimumPitch = 97.9989; // 97.9989 Hz is G2 //in Hz, minimum frequency required to detect
    _maximumPitchToDetect =  587.330; // 587.330 Hz is D5, 493.883 Hz is B4
    _voicedUnvoicedCost = 0.1; // default 0.2
    _octaveJumpCost = 1.5; // default 0.2.
```   
    
# what is implemented: 
 * Hanning windowing
 * pitch detection based on autocorrelation  
 * proper autocorrelation normalization due to windowing
 * Viterbi search algorithm to suppress false pitch identification

# what is not implemented (but would be great):
 * Gaussian window
 * sinc interpolation

# performance:
pitch is accurate to about 1% level (need sinc interploation for better accuracy)


