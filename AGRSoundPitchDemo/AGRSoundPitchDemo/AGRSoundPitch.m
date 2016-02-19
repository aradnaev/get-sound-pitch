
//
//  Created by Alex Radnaev on 2/16/15.
//  Copyright (c) 2015 Radnaev.com. All rights reserved.
//
//Algorithm based on IFA Proceedings 17, 1993
//ACCURATE SHORT-TERM ANALYSIS OF THE FUNDAMENTAL FREQUENCY AND THE HARMONICS-TO-NOISE RATIO OF A SAMPLED SOUND
//Paul Boersma

#import "AGRSoundPitch.h"

@implementation AGRSoundPitch {

    COMPLEX_SPLIT _A;
    COMPLEX_SPLIT _APSD;
    FFTSetup      _FFTSetup;
    FFTSetup      _FFTSetupPSD;
    vDSP_Length   _log2n;
    vDSP_Length   _log2nPSD;
    float *window; // e.g. Hanning window
    float *r_window; // normalized autocorrelation of the window
    int samplingWindowType;
    int nOver2;

    int maximumNumberOfCandidatesPerFrame;
    //int nOver2PSD;
    int frameLengthWithZeros;
    float *dataFrame;
    float *dataFramePSD;
    float *normalizedAutoCorrelation;
    float *normalizedAutoCorrelationManual;
    float *normalizedAutoCorrelationDiff;
    float *maximaValues;
    float *maximaFrequencies;
    float *mirrorredData;
    int newBufferSize;
    
}


-(instancetype)initWithBufferSize:(UInt32)bufferSize{
  self = [super init];
  NSLog(@"initializing VoiceQuanta with bufferSize = %i",bufferSize);
    //settings derived experimentally for iPhone 6:
    //samping rate 44100 Hz
    //buffer size 1024
    //23.22 ms buffer duration (43.07 Hz)
    //11.61 ms autocorrelation length (86.13 Hz)
    //minimum pitch 90 Hz
    //silenceThreshold = 0.02
    //voicingThreshold = 0.05
    //octaveCost = 0.15^2
    //voiceUnvoiceCost = 0.1
    //octabeJumpCost = 1
    
    //PARAMETERS:
    _hearingThreshold = 0.00001;
    _hearingThresholdDB = [self amplitudeToDecibel:_hearingThreshold];
    samplingWindowType = 1; //1 - Hanning. to-do: implement Gaussian filter
    maximumNumberOfCandidatesPerFrame = 2000;
    _voicingThreshold = 0.4; //%VoicingThreshold = 0.4 by default
    _silenceThreshold = 0.1; //%silence threshold = 0.05 by default
    _octaveCost = pow(0.15,2); //default is 0.01 (10% = sqrt(0.01)) // the higher octaveCost, the higher weight of higher frequencies
//    The OctaveCost parameter favours higher fundamental frequencies. One of the reasons for the existence of this parameter is that for a perfectly periodic signal all the peaks are equally high and we should choose the one with the lowest lag. Other reasons for this parameter are unwanted local downward octave jumps caused by additive noise (section 6). Finally, an important use of this parameter lies in the difference between the acoustic fundamental frequency and the perceived pitch. 
    _minimumPitch = 97.9989; // 97.9989 Hz is G2 //in Hz, minimum frequency required to detect
    _maximumPitchToDetect =  587.330; // 587.330 Hz is D5, 493.883 Hz is B4
    _voicedUnvoicedCost = 0.1; // default 0.2
    _octaveJumpCost = 1.5; // default 0.2.
    _minimumVolume = 0.001;
    _maximumVolume = 1;
    
    
    //INIT:
    _pitchCandidatesFrequencies = (float*) malloc(maximumNumberOfCandidatesPerFrame*sizeof(float));
    _pitchCandidatesStrengths = (float*) malloc(maximumNumberOfCandidatesPerFrame*sizeof(float));
    
    //to-do: move the init code to init section
    nOver2 = bufferSize/2;

    frameLengthWithZeros = round(bufferSize*1.5);
    frameLengthWithZeros = pow(2,(ceil(log2(frameLengthWithZeros))));

    dataFrame = (float *) malloc(frameLengthWithZeros*sizeof(float));
    dataFramePSD = (float *) malloc(nOver2*sizeof(float));
    normalizedAutoCorrelation = (float *) malloc(nOver2*sizeof(float));
    normalizedAutoCorrelationManual = (float *) malloc(nOver2*sizeof(float));
    normalizedAutoCorrelationDiff = (float *) malloc((nOver2-1)*sizeof(float));
    maximaValues = (float *) malloc(maximumNumberOfCandidatesPerFrame*sizeof(float));
    maximaFrequencies = (float *) malloc(maximumNumberOfCandidatesPerFrame*sizeof(float));
    newBufferSize = bufferSize*2 - 1;

    mirrorredData = (float *)malloc(sizeof(float)*newBufferSize);
    
    if (self) {
        [self createFFTWithBufferSize:bufferSize withAudioData:nil];
        [self createSamplingWindowWithBufferSize:bufferSize WithWindowType:samplingWindowType];
        _FFTisSetup = true;
        for (int i=0; i<maximumNumberOfCandidatesPerFrame; i++) {
            _pitchCandidatesFrequencies[i] = 0.0;
            _pitchCandidatesStrengths[i] = 0.0;
            _pitchCandidatesCount = 1;
        }
    }
      NSLog(@"initializing VoiceQuanta done");
    return self;
    
}

-(void)createSamplingWindowWithBufferSize:(UInt32)bufferSize WithWindowType:(int)windowType{
    window = (float *) malloc(bufferSize*sizeof(float));
    r_window = (float *) malloc(bufferSize*sizeof(float));
    //window calculaton

 
    for(int ti=0; ti<bufferSize; ti++) {
       // window[ti] = (exp(-12*pow(2, ti/bufferSize - 0.5)) )- exp(-12))/(1- exp(-12))  ;
        if (windowType == 1) {//Hanning
            window[ti] = 0.5-0.5*cos(2*M_PI*ti/bufferSize);
            r_window[ti] = (1.0-ti/bufferSize)*( 2.0/3.0+(1.0/3.0)*cos(2.0*M_PI*ti/bufferSize)) +(1.0/(2.0*M_PI))*sin(2.0*M_PI*ti/bufferSize);
            //NSLog(@"r_w: %1.5f, ti: %d, test: %1.5f\n",r_window[ti],ti,1.0-(ti/bufferSize));
        } else {
            window[ti] = 1.0;
            r_window[ti] = 1.0;
        }
    }
    NSLog(@"initializing Window done");
}

-(void)createFFTWithBufferSize:(UInt32)bufferSize withAudioData:(float*)data {
    NSLog(@"Creating FFT with bufferSize = %i",bufferSize);
    // Setup the length
    _log2n = log2f(bufferSize);
    
    // Calculate the weights array. This is a one-off operation.
    _FFTSetup = vDSP_create_fftsetup(_log2n, FFT_RADIX2);
    
    // For an FFT, numSamples must be a power of 2, i.e. is always even

    // Populate *window with the values for a hamming window function
    //    float *_window = (float *)malloc(sizeof(float)*bufferSize);
    //  vDSP_hamm_window(_window, bufferSize, 0);
    // Window the samples
    //vDSP_vmul(data, 1, window, 1, data, 1, bufferSize);
    // free(window);
    
    // Define complex buffer
    _A.realp = (float *) malloc(nOver2*sizeof(float));
    _A.imagp = (float *) malloc(nOver2*sizeof(float));
    
    _log2nPSD = log2f(nOver2);
    _FFTSetupPSD = vDSP_create_fftsetup(_log2nPSD, FFT_RADIX2);

    _APSD.realp = (float *) malloc(nOver2*sizeof(float));
    _APSD.imagp = (float *) malloc(nOver2*sizeof(float));
    NSLog(@"initializing FFT done");
    
}

-(float)amplitudeToDecibel:(float)amplitude{
    return 20*log10f(amplitude/_hearingThreshold);
}


-(float*)mirrorData:(float*)data withBufferSize:(UInt32)bufferSize {
    //    int newBufferSize = bufferSize;

    for (int i = 0; i < bufferSize; i++) {
        mirrorredData[i] = data[bufferSize - i - 1];
        mirrorredData[i + bufferSize] = data[i];
        // mirrorredData[i] = data[i];
    }
    return mirrorredData;
}

-(float)processData:(float*)data WithBufferSize:(UInt32)bufferSize WithSamplingRate:(Float64)samplingRate WithGlobalAbsolutePeak:(float*)globalAbsolutePeak WithFFTplotPointer:(float*)FFTdataHolder WithACPlotPointer:(float*)ACPlotHolder {


    int maximumPitch = fmin(_maximumPitchToDetect,samplingRate/2);// %in Hz, maximum frequency required to detect, limited by Nyquist frequency
    float carrierFreq = 0;

    // volume calculation method=1
    float rms = 0.0;
    Float32 localAbsolutePeak = 0.0;
    Float32 sum = 0.0;
    for(int i=0; i<bufferSize; i++) {
        // Calculate the RMS magnitude
        rms = rms + (data[i])*(data[i]);
        
        // calculate absolute magnitude. to-do - vDSP_maxv a cleaner and faster way
        localAbsolutePeak = data[i] > localAbsolutePeak ? data[i] : localAbsolutePeak;
        
        sum = sum + data[i];
    }
    if (localAbsolutePeak > *globalAbsolutePeak ) {
        *globalAbsolutePeak = localAbsolutePeak;
     //   NSLog(@"higher peak detected \n");
    }
    float gain = 0.005;
    *globalAbsolutePeak = *globalAbsolutePeak + gain*(localAbsolutePeak - *globalAbsolutePeak);
   // NSLog(@"globalAbsolutePeak: %1.5f \n",*globalAbsolutePeak);
    float meanValue = sum/bufferSize;
    float localACAbsolutePeak = 0;
    for(int i=0; i<bufferSize; i++) {
        dataFrame[i] = (data[i] - meanValue)*window[i]; //to-do - use vDSP_vmul
        localACAbsolutePeak = dataFrame[i] > localACAbsolutePeak ? dataFrame[i] : localACAbsolutePeak;
    }
    self.volume = [self amplitudeToDecibel:sqrt(rms)/bufferSize];
    //NSLog([NSString stringWithFormat:@"volume: %1.5f",self.volume]);
    float mag =0.0;
    float maxMag = 0.0;
    int maximaCount=0;
    if (localAbsolutePeak == 0.0){
        //%no signal - zero pitch
        for (int m=0; m<maximumNumberOfCandidatesPerFrame; m++){
            _pitchCandidatesFrequencies[m] = 0;
            _pitchCandidatesStrengths[m] = 0.0; //_voicingThreshold + 2;            //NSLog(@"local max = %1.1f\n",localAbsolutePeak);
            _pitchCandidatesCount = 1;
        }
        
        carrierFreq = 0.0;
        for(int i=0; i<nOver2; i++) {
            FFTdataHolder[i] = 0.0;
            ACPlotHolder[i] = 0.0;
        }
        for(int i=nOver2; i<bufferSize; i++) {
            ACPlotHolder[i] = 0.0;
        }
    } else { // non zero signal - go!
        ////////////////////////////////////////////////
        // calculate autocorrelation function by FFT ////
        ////////////////////////////////////////////////
        
        // Pack samples:    // C(re) -> A[n], C(im) -> A[n+1]
        vDSP_ctoz((COMPLEX*)dataFrame, 2, &_A, 1, nOver2);
        // Perform a forward FFT using fftSetup and A     // Results are returned in A
        vDSP_fft_zrip(_FFTSetup, &_A, 1, _log2n, FFT_FORWARD);
        // Convert COMPLEX_SPLIT A result to magnitudes

       
        for(int i=0; i<nOver2; i++) {
            // Calculate the magnitude
            dataFramePSD[i] = _A.realp[i]*_A.realp[i]+_A.imagp[i]*_A.imagp[i];
            mag = _A.realp[i]*_A.realp[i]+_A.imagp[i]*_A.imagp[i];
            if (mag > maxMag) {carrierFreq = i;}
            maxMag = mag > maxMag ? mag : maxMag;
            if (FFTdataHolder) {
                FFTdataHolder[i] = 0.9*mag/maxMag;
            }


        }
        
        //   NSLog(@"FFT updated\n");
        //   %calculate auto correllation using inverse FFT of Power spectral
        //   %density (much faster than brute force calculation)
        
//        vDSP_ctoz((COMPLEX*)dataFramePSD, 2, &_APSD, 1, nOver2/2);
//        // Perform a forward FFT using fftSetup and A     // Results are returned in A
//        vDSP_fft_zrip(_FFTSetupPSD, &_APSD, 1, _log2nPSD, FFT_INVERSE);
//        
        
        //manual calculation of autocorrelation
        //TODO: calculate autocorreation by two FFT transforms
        maxMag = 0.0;
        Float32 acsum=0.0;
        for (int tau=0; tau<nOver2; tau++) {
            acsum = 0.0;
            for(int t=0;t<bufferSize-tau;t++){
                acsum = acsum + dataFrame[t]*dataFrame[t+tau];
            }
            normalizedAutoCorrelationManual[tau]=acsum;
            mag =acsum;
            maxMag = mag > maxMag ? mag : maxMag;
            
        }
        float normalizedACValue = 0;
        float firstACvalue =normalizedAutoCorrelationManual[0];
        for(int i=0; i<nOver2; i++) {  //nOver2 = bufferSize/2;
            normalizedACValue =normalizedAutoCorrelationManual[i]/firstACvalue;
            normalizedAutoCorrelationManual[i] = normalizedACValue/r_window[i];
            //normalizedAutoCorrelationManual[i] =  normalizedAutoCorrelationManual[i];///(maxMag*r_window[i]);
            
            normalizedAutoCorrelation[i] = normalizedAutoCorrelationManual[i];
            
            if (ACPlotHolder) {
                ACPlotHolder[i + nOver2] = 0.9*normalizedAutoCorrelationManual[i];
                ACPlotHolder[nOver2 - i - 1] = ACPlotHolder[i + nOver2];
            }
        }
 
        
        ///// DONE WITH AUTOCORRELATION ////

      /// PEAK SEARCH IN AUTOCORRELATION AND SORTING BY MAGNITUDE///
    /// OUTPUT: maximaCount, maximaValues[], maximaFrequencies[], maximaIndicies[] sorted by maximaValue
        bool positiveSlope = 1; //has to be 1 to include first datapoint as maxium

        float strongestCandidateAmplitude = 0;
        int strongestCandidateIndex = 0;
        float smallestMaximumValue = 1;
        int smallestMaximumIndexWithinCandidates = 0;
        for (int i=1; i<nOver2; i++){

                normalizedAutoCorrelationDiff[i-1]=normalizedAutoCorrelation[i]-normalizedAutoCorrelation[i-1];
                if(normalizedAutoCorrelationDiff[i-1] >= 0 ){
                    positiveSlope = 1;
                } else {
                    if (normalizedAutoCorrelationDiff[i-1] < 0 && positiveSlope == 1) { // slope changed from positive to negative --> maximum
                        maximaFrequencies[maximaCount] = samplingRate/(i+1);
                        if( maximaFrequencies[maximaCount]>_minimumPitch &&  maximaFrequencies[maximaCount]<maximumPitch) { // frequency within range
                            //found a new maximum within range but don't increase maxiumCount, because it's is used as zero-based index
                            maximaValues[maximaCount] = normalizedAutoCorrelation[i-1] - _octaveCost*log2(_minimumPitch/maximaFrequencies[maximaCount]);

                            if (maximaCount+1>maximumNumberOfCandidatesPerFrame-1){ //1 is the unvoiced candidate, will be added later
                                // if not within budget of candidates, check if found a stronger one
                                if(maximaValues[maximaCount]>smallestMaximumValue){ // stronger indeed --> replace
                                    maximaValues[smallestMaximumIndexWithinCandidates] = maximaValues[maximaCount];
                                    maximaFrequencies[smallestMaximumIndexWithinCandidates] =maximaFrequencies[maximaCount];
                                    if (maximaValues[maximaCount]>strongestCandidateAmplitude) { //also check if we found the strongest candidate
                                        strongestCandidateAmplitude = maximaValues[maximaCount]; //record strongest value
                                        strongestCandidateIndex = smallestMaximumIndexWithinCandidates; // and index
                                    }
                                }

                            } else { //still within budget of number of candidates, update smallest Max value
                                if (maximaValues[maximaCount]<smallestMaximumValue){
                                    smallestMaximumValue = maximaValues[maximaCount];
                                    smallestMaximumIndexWithinCandidates = maximaCount;
                                }
                                if (maximaValues[maximaCount]>strongestCandidateAmplitude) { //also check if found strongest candidate
                                    strongestCandidateAmplitude = maximaValues[maximaCount];
                                    strongestCandidateIndex = maximaCount;
                                }
                                maximaCount++; //increase number of maxima found only when still within budget
                            }
                            
                                
                        }
                        positiveSlope = 0;
                    }//end of maximum found block
                } // end of negative slope detection

        } // end of for loop for peak search in autocorrelation

        //add unvoiced candidate
        maximaValues[maximaCount] =_voicingThreshold + MAX(0,2-(localACAbsolutePeak/(*globalAbsolutePeak))/(_silenceThreshold/(1+_voicingThreshold)));
//        float unvoicedCandidateStrength = maximaValues[maximaCount];
//        float candidate1Strength = maximaValues[0];
//        float candidate2Strength = maximaValues[1];
//        float candidate3Strength = maximaValues[2];
//
//        float candidate1Frequency = maximaFrequencies[0];
//        float candidate2Frequency = maximaFrequencies[1];
//        float candidate3Frequency = maximaFrequencies[2];
//        
        maximaFrequencies[maximaCount] = 0;
        
        if (maximaValues[maximaCount]>strongestCandidateAmplitude) { //also check if found strongest candidate
            strongestCandidateAmplitude = maximaValues[maximaCount];
            strongestCandidateIndex = maximaCount;
        }
        maximaCount = maximaCount + 1;
//        carrierFreq = maximaFrequencies[strongestCandidateIndex];

        //VITERBI SEARCH ///

        float switchCost;
        float lowestCost = 100000.0;
        //float *transitionCost = (float *) malloc(maximumNumberOfCandidatesPerFrame*maximumNumberOfCandidatesPerFrame*sizeof(float));
        float transitionCost;

        for (int k=0;k<MIN(_pitchCandidatesCount,maximumNumberOfCandidatesPerFrame);k++){ // k is index for previous frame pitch candidates
            for (int n=0;n<MIN(maximaCount,maximumNumberOfCandidatesPerFrame);n++){ // n is index for current frame pitch candidates
                if((_pitchCandidatesFrequencies[k])*(maximaFrequencies[n]) == 0 ){
                    switchCost = _voicedUnvoicedCost;  // if either candidate is unvoiced then special case cost
                }else{
                    switchCost = _octaveJumpCost*fabsf(log2f(_pitchCandidatesFrequencies[k]/maximaFrequencies[n]));
                }
               transitionCost = switchCost - maximaValues[n] - _pitchCandidatesStrengths[k];
//                transitionCost[(n-1)*maximumNumberOfCandidatesPerFrame+k]
                if (transitionCost<lowestCost) {
                    lowestCost = transitionCost;
                    strongestCandidateIndex  = k;
                }      
            }
        }
        // Viterbi search is done --> update previous pitch and pass pitches for the next frame
        self.previousPitch = _pitchCandidatesFrequencies[strongestCandidateIndex];
        carrierFreq = self.previousPitch;
        
        for (int m=0; m<maximumNumberOfCandidatesPerFrame; m++){
            _pitchCandidatesFrequencies[m] =maximaFrequencies[m];
            _pitchCandidatesStrengths[m] = maximaValues[m];
            _pitchCandidatesCount = maximaCount;
        }
 

    } // end of (localAbsolutePeak != 0))
    
//        // to-do see if better to use vDSP_conv(x, 1, x, 1, result, 1, 2*len_X-1, len_X);

    self.currentPitch = carrierFreq;
    return self.previousPitch;
    
}


@end
