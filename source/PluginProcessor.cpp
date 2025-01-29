#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "leaf-filters.h"
#include "leaf-oscillators.h"
BirlOneHoleAudioProcessor::BirlOneHoleAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
    : AudioProcessor (BusesProperties()
    #if ! JucePlugin_IsMidiEffect
        #if ! JucePlugin_IsSynth
                          .withInput  ("Input",  AudioChannelSet::stereo(), true)
        #endif
                          .withOutput ("Output", AudioChannelSet::stereo(), true)
    #endif
              )
      , tree (*this, nullptr, "PARAMETERS", {
                                              std::make_unique<AudioParameterInt> ("tuning", "Tuning", 0, 4, 0),
                                          })

#endif
{
    birl.setNumToneHoles(NUM_OF_KEYS);
    birl.setShaperDrive(0.0);
    birl.setNoiseBPQ(0.5);
    birl.setNoiseGain(0.5);
    birl.setNoiseBPCutoff(200.0);

}

BirlOneHoleAudioProcessor::~BirlOneHoleAudioProcessor()
{
    //body
}

//==============================================================================
const String BirlOneHoleAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool BirlOneHoleAudioProcessor::acceptsMidi() const
{
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}

bool BirlOneHoleAudioProcessor::producesMidi() const
{
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}

bool BirlOneHoleAudioProcessor::isMidiEffect() const
{
#if JucePlugin_IsMidiEffect
    return true;
#else
    return false;
#endif
}

double BirlOneHoleAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int BirlOneHoleAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
        // so this should be at least 1, even if you're not really implementing programs.
}

int BirlOneHoleAudioProcessor::getCurrentProgram()
{
    return 0;
}

void BirlOneHoleAudioProcessor::setCurrentProgram (int index)
{
}

const String BirlOneHoleAudioProcessor::getProgramName (int index)
{
    return {};
}

void BirlOneHoleAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void BirlOneHoleAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    LEAF_init(&leaf, sampleRate, myMemory, 500000, []() {return (float)rand() / RAND_MAX; });
    tSVF_init(&birl.inLoopLP, SVFTypeLowpass, 3000.0f, 0.7f,&leaf);
    tCycle_init(&sin,&leaf);
    tADSR_init(&breathPosEnv, 500, 500, 1000, 1000,&leaf);

    /* These are for the X-Y keys */
    tExpSmooth_init(&XY1_XEnv, 0.0f, 0.02f,&leaf); // 0.02 = 5 ms
    tExpSmooth_init(&XY1_YEnv, 0.0f, 0.02f,&leaf);
    tExpSmooth_init(&XY2_XEnv, 0.0f, 0.02f, &leaf);
    tExpSmooth_init(&XY2_YEnv, 0.0f, 0.02f,&leaf);

    for (int i = 0; i < NUM_OF_KEYS; ++i){
        tExpSmooth_init(&keysEnv[i], 0.0f, 0.01f, &leaf);
    }
}

void BirlOneHoleAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool BirlOneHoleAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
    #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
        && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

            // This checks if the input layout matches the output layout
        #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
        #endif

    return true;
    #endif
}
#endif
void BirlOneHoleAudioProcessor::synthesize() {
    //    birl.setToneHoleIndex(0);
    //    birl.setToneHole(0.0-tExpSmooth_tick(&keysEnv[0]));
    //    birl.setToneHoleIndex(1);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[1]));
    //    birl.setToneHoleIndex(2);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[2]));
    //    birl.setToneHoleIndex(3);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[3]));
    //    birl.setToneHoleIndex(4);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[4]));
    //    birl.setToneHoleIndex(5);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[5]));
    //    birl.setToneHoleIndex(6);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[6]));
    //    birl.setToneHoleIndex(7);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[7]));
    //    birl.setToneHoleIndex(8);
    //    birl.setToneHole(1.0-tExpSmooth_tick(&keysEnv[8]));
}

void BirlOneHoleAudioProcessor::parameterSmooth(){
    tExpSmooth_setDest(XY1_XEnv, 1.0);
    tExpSmooth_setDest(XY1_YEnv, 1.0);
    tExpSmooth_setDest(XY2_XEnv, 1.0);
    tExpSmooth_setDest(XY2_YEnv, 1.0);

    //    for (int i = 0; i < NUM_OF_KEYS; ++i) {
    //        tExpSmooth_setDest(&keysEnv[i], keysGradient[i]);
    //    }
}

void BirlOneHoleAudioProcessor::changeTuning(){
    int selectedTuning = *tree.getRawParameterValue("tuning");
    switch (selectedTuning) {
        case 0:
            ::setTuning(JUST_INTONATION);
            break;
        case 1:
            ::setTuning(EQUAL_TEMPERED);
            break;
        case 3:
            ::setTuning(MEANTONE);
            break;
        case 4:
            ::setTuning(HIGHLAND_BAGPIPE);
            break;
        default:
            break;
    }
}
void BirlOneHoleAudioProcessor::scaleValues() {
    for (int tonehole = 0; tonehole < NUM_OF_KEYS; ++tonehole) {
        if (keys_[tonehole] > maxKeyArg[tonehole]) {
            maxKeyArg[tonehole] = keys_[tonehole];
        }
        if (keys_[tonehole] < minKeyArg[tonehole]) {
            minKeyArg[tonehole] = keys_[tonehole];
        }
        keysGradient[tonehole] = (double)(keys_[tonehole] - minKeyArg[tonehole]) / (double)(maxKeyArg[tonehole] - minKeyArg[tonehole]);
    }
}
void BirlOneHoleAudioProcessor::buttonListen() {
    // listen for buttons
}
void BirlOneHoleAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    buffer.clear();
    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    Array<float*> channelData;
    for (int channel = 0; channel < totalNumInputChannels; ++channel)
    {
        channelData.add(buffer.getWritePointer (channel));
    }

    for (int sampIndex = 0; sampIndex < buffer.getNumSamples(); ++sampIndex) {
        //        scaleValues();
        parameterSmooth();
        //        synthesize();
        birl.setBreathPressure(gain);
        audioBufferTick = birl.birlTick();

        for (int j = 0; j < totalNumInputChannels; ++j) {
            channelData.getUnchecked(j)[sampIndex] = audioBufferTick;
        }
    }
}

bool BirlOneHoleAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* BirlOneHoleAudioProcessor::createEditor()
{
    return new BirlOneHoleAudioProcessorEditor (*this);
}

void BirlOneHoleAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void BirlOneHoleAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new BirlOneHoleAudioProcessor();
}



