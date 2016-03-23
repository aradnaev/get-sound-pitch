//
//  ViewController.m
//  AGRSoundPitchDemo
//
//  Created by Alex Radnaev on 2/19/16.
//  Copyright © 2016 Radnaev.com. All rights reserved.
//
#import "AGRSoundPitch.h"
#import "ViewController.h"
#import "AEAudioController.h"
#import "AEUtilities.h"
#import "AEFloatConverter.h"
#import "AEBlockAudioReceiver.h"

@import AVFoundation;

@interface ViewController ()
@property (nonatomic, strong) AEAudioController *audioController;
@property (weak, nonatomic) IBOutlet UILabel *pitchLabel;
@property (nonatomic, assign) int longBufferSize;
@property (nonatomic, assign) float sampleRate;
@property (nonatomic, strong) AGRSoundPitch *soundPitch;
@end

@implementation ViewController{
    int accumulatedBufferLength;
    float *longBuffer;
    float *globalAbsolutePeak;


}


-(void)updateVoiceBuffer:(float**)buffer withBuffersize:(UInt32)bufferSize {
    if(bufferSize == 0){
        NSLog(@"processReceivedData: nothing to process - bufferSize is zero");
        return;
    }
    if(bufferSize > self.longBufferSize){
        NSLog(@"processReceivedData error: received buffer of %u is too long (longer than allocated buffer size of %u)",bufferSize,self.longBufferSize);
        return;
    }
    if(self.longBufferSize % bufferSize !=0 ){
        NSLog(@"processReceivedData error: received buffer is subinteger of long buffer if size %u",self.longBufferSize);
        return;
    }
    
    //TODO: implement circular buffer http://atastypixel.com/blog/a-simple-fast-circular-buffer-implementation-for-audio-processing/
    
    if(accumulatedBufferLength<self.longBufferSize) { // long buffer is not filled yet
        for(int i=0; i<bufferSize; i++){
            longBuffer[i+accumulatedBufferLength] = buffer[0][i];
        }
        accumulatedBufferLength = accumulatedBufferLength + bufferSize;
    }
    
    if(bufferSize == self.longBufferSize){
        accumulatedBufferLength = self.longBufferSize;
    }
    
    if(accumulatedBufferLength==self.longBufferSize)
    {
        for(int i=0; i<self.longBufferSize-bufferSize; i++){
            longBuffer[i] = longBuffer[i+bufferSize]; //move all data to the left by amount of new buffer
        }
        for(int i=0; i<bufferSize; i++){
            longBuffer[i+self.longBufferSize-bufferSize] = buffer[0][i];
        }
    }
    
    
}


- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    
    _longBufferSize = 1024;
    globalAbsolutePeak = (float*)malloc(sizeof(float));
    longBuffer = (float*)malloc(sizeof(float)*self.longBufferSize);
    *globalAbsolutePeak = 0.0;

    
    AudioStreamBasicDescription exerciseInputAudioDescription = AEAudioStreamBasicDescriptionMake(AEAudioStreamBasicDescriptionSampleTypeFloat32 , NO, 1, 44100);
    _sampleRate = exerciseInputAudioDescription.mSampleRate; //TODO: set sample rate
    
    

    _soundPitch = [[AGRSoundPitch alloc] init];
    
    [self.soundPitch initWithBufferSize:self.longBufferSize];
    
    self.audioController = [[AEAudioController alloc]
                            initWithAudioDescription:exerciseInputAudioDescription
                            inputEnabled:YES]; // don't forget to autorelease if you don't use ARC!
    
    [self.audioController stop];
    
    AEFloatConverter *exerciseAudioFloatConverter = [[AEFloatConverter new] initWithSourceFormat:exerciseInputAudioDescription];
    
    
    __weak ViewController *weakSelf = self;
    //id<AEAudioReceiver> receiver = [[AGRAudioReceiver alloc] init]
    id<AEAudioReceiver> receiver = [AEBlockAudioReceiver audioReceiverWithBlock:
                                    ^(void *source,
                                      const AudioTimeStamp *time,
                                      UInt32 frames,
                                      AudioBufferList *audio) {
                                        
                                        int nChannels = self.audioController.numberOfInputChannels;
                                        // NSLog(@"Received audio with %u frames, number of channels %u",(unsigned int)frames,nChannels);
                                        static float **scratchBuffer = NULL; // static will make sure this variable is alive even when the block is complete for next block reusal
                                        //                                      static const UInt32 scratchBufferFrames = 512;
                                        //                                    static float scratchBuffer[2][scratchBufferFrames]; //this allocates memory
                                        //TODO: if number of frames is variable, create static UINT32 scratchBufferFrames and monitor if need to reallocate memory if the frames number has changed
                                        if(scratchBuffer == NULL){
                                            scratchBuffer = (float**)malloc(sizeof(float*) * nChannels);
                                            assert(scratchBuffer);
                                            for (int c = 0; c < nChannels; c++) {
                                                scratchBuffer[c] = (float *)malloc(sizeof(float)*frames);
                                                assert(scratchBuffer[c]);
                                            }
                                        }
                                        AEFloatConverterToFloat(exerciseAudioFloatConverter, audio, scratchBuffer, frames);
                                        [self updateVoiceBuffer:scratchBuffer withBuffersize:frames];
                                        
                                        if(accumulatedBufferLength >= self.longBufferSize)
                                        {
                                            
                                            [weakSelf.soundPitch processData:longBuffer WithBufferSize:weakSelf.longBufferSize WithSamplingRate:self.sampleRate WithGlobalAbsolutePeak:globalAbsolutePeak WithFFTplotPointer:nil WithACPlotPointer:nil];
                                            
                                            dispatch_async(dispatch_get_main_queue(), ^{
                                                [weakSelf.pitchLabel setText:[NSString stringWithFormat:@"%u Hz",(int)round(weakSelf.soundPitch.previousPitch)]];
                                                
                                                
                                            }
                                                           );

                                        }
                                    }];
    [self.audioController addInputReceiver:receiver];
    NSError *error;
    [self.audioController start:&error];
    if (error) {
        NSLog(@"error: %@",[error localizedDescription]);
    }
    
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
