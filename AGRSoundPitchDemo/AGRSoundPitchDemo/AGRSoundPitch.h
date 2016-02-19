//
//  Created by Alex Radnaev on 2/16/15.
//  Copyright (c) 2015 Radnaev.com. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>
#import <math.h>

@interface AGRSoundPitch : NSObject

-(instancetype)initWithBufferSize:(UInt32)bufferSize;
-(float)processData:(float*)data WithBufferSize:(UInt32)bufferSize WithSamplingRate:(Float64)samplingRate WithGlobalAbsolutePeak:(float*)globalAbsolutePeak WithFFTplotPointer:(float*)FFTdataHolder WithACPlotPointer:(float*)ACPlotHolder;

@property (nonatomic,assign) BOOL FFTisSetup;
@property (nonatomic,assign) float* pitchCandidatesFrequencies;
@property (nonatomic,assign) float* pitchCandidatesStrengths;
@property (nonatomic,assign) int pitchCandidatesCount;
@property (nonatomic,assign) float volume;
@property (nonatomic,assign) float currentPitch;
@property (nonatomic,assign) float previousPitch;

@property (nonatomic,assign) float maximumPitchToDetect;
@property (nonatomic,assign) float minimumPitch;
@property (nonatomic,assign) float minimumVolume;
@property (nonatomic,assign) float maximumVolume;

@property (nonatomic,assign) float voicingThreshold;
@property (nonatomic,assign) float silenceThreshold;
@property (nonatomic,assign) float octaveCost;
@property (nonatomic,assign) float voicedUnvoicedCost;
@property (nonatomic,assign) float octaveJumpCost;
@property (nonatomic,assign) float hearingThreshold;
@property (nonatomic,assign) float hearingThresholdDB;

@end
