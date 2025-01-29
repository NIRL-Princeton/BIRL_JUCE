#include "PluginProcessor.h"
#include "PluginEditor.h"
#include "leaf.h"
#include "melatonin_audio_sparklines/melatonin_audio_sparklines.h"
BirlAudioProcessor::BirlAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       ), parameters(*this, nullptr, "PARAMETERS", {
         std::make_unique<juce::AudioParameterFloat>("gain", "Gain", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("fundamental", "Fundamental", NormalisableRange<float>(20.0f, 5000.0f), 100.0f),
         std::make_unique<juce::AudioParameterFloat>("num_toneholes", "Num_Toneholes", NormalisableRange<float>(0.0f, 9.0f), 9.0f),
         std::make_unique<juce::AudioParameterFloat>("dcblocker1", "DCBlocker1", NormalisableRange<float>(0.0f, 1.0f), 0.995f),
         std::make_unique<juce::AudioParameterFloat>("dcblocker2", "DCBlocker2", NormalisableRange<float>(0.0f, 1.0f), 0.995f),
         std::make_unique<juce::AudioParameterFloat>("biquad_coeff1", "Biquad_Coeff1", NormalisableRange<float>(-1.0f, 1.0f), 0.169301f),
         std::make_unique<juce::AudioParameterFloat>("biquad_coeff2", "Biquad_Coeff2", NormalisableRange<float>(-1.0f, 1.0f), 0.338601f),
         std::make_unique<juce::AudioParameterFloat>("biquad_coeff3", "Biquad_Coeff3", NormalisableRange<float>(-1.0f, 1.0f), 0.169301f),
         std::make_unique<juce::AudioParameterFloat>("biquad_coeff4", "Biquad_Coeff4", NormalisableRange<float>(-1.0f, 1.0f), -0.482013f),
         std::make_unique<juce::AudioParameterFloat>("biquad_coeff5", "Biquad_Coeff5", NormalisableRange<float>(-1.0f, 1.0f), 0.186622f),
         std::make_unique<juce::AudioParameterFloat>("pf1_cutoff", "PF1_Cutoff", NormalisableRange<float>(0.0f, 20000.0f), 2000.0f),
         std::make_unique<juce::AudioParameterFloat>("pf1_q", "PF1_Q", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("pf2_cutoff", "PF2_Cutoff", NormalisableRange<float>(0.0f, 20000.0f), 5000.0f),
         std::make_unique<juce::AudioParameterFloat>("pf2_q", "PF2_Q", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("lp1_cutoff", "LP1_Cutoff", NormalisableRange<float>(0.0f, 20000.0f), 5000.0f),
         std::make_unique<juce::AudioParameterFloat>("lp1_q", "LP1_Q", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("lp2_cutoff", "LP2_Cutoff", NormalisableRange<float>(0.0f, 20000.0f), 5000.0f),
         std::make_unique<juce::AudioParameterFloat>("lp2_q", "LP2_Q", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("shaper_drive", "Shaper_Drive", NormalisableRange<float>(0.0f, 1.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("shaper_mix", "Shaper_Mix", NormalisableRange<float>(0.0f, 1.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("noise_gain", "Noise_Gain", NormalisableRange<float>(0.0f, 1.0), 0.2),
         std::make_unique<juce::AudioParameterFloat>("noise_bp_cutoff", "Noise_BP_Cutoff", NormalisableRange<float>(0.0f, 20000.0f), 16000.0f),
         std::make_unique<juce::AudioParameterFloat>("noise_bp_q", "Noise_BP_Q", NormalisableRange<float>(0.0f, 1.0f), 1.0f),
         
         
         std::make_unique<juce::AudioParameterFloat>("osc1_on", "Osc1_On", NormalisableRange<float>(0.0f, 1.0f), 1.0f), // true
         std::make_unique<juce::AudioParameterFloat>("osc1_waveform", "Osc1_Waveform", NormalisableRange<float>(0.0f, 3.0f), 1.0f),
         std::make_unique<juce::AudioParameterFloat>("osc1_gain", "Osc1_Gain", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("osc1_octave", "Osc1_Octave", NormalisableRange<float>(-3.0f, 3.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("osc1_semi", "Osc1_Semi", NormalisableRange<float>(-12.0f, 12.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("osc1_detune", "Osc1_Detune", NormalisableRange<float>(-3.0f, 3.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("filter1_on", "Filter1_On", NormalisableRange<float>(0.0f, 1.0f), 1.0f), // true
         std::make_unique<juce::AudioParameterFloat>("filter1_type", "Filter1_Type", NormalisableRange<float>(0.0f, 4.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("filter1_tof2", "Filter1_toF2", NormalisableRange<float>(0.0f, 1.0f), 0.0f), // in percentage
         std::make_unique<juce::AudioParameterFloat>("filter1_freq", "Filter1_Freq", NormalisableRange<float>(20.0f, 22000.0f), 22000.0f),
         std::make_unique<juce::AudioParameterFloat>("filter1_q", "Filter1_Q", NormalisableRange<float>(0.0f, 1.0f), 0.7f),
         std::make_unique<juce::AudioParameterFloat>("amp1_on", "Amp1_On", NormalisableRange<float>(0.0f, 1.0f), 1.0f), // true
         std::make_unique<juce::AudioParameterFloat>("amp1_pan", "Amp1_Pan", NormalisableRange<float>(-1.0f, 1.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("amp1_level", "Amp1_Level", NormalisableRange<float>(-40.0f, 6.0f), 0.0f), // this is in dB!!!
         std::make_unique<juce::AudioParameterFloat>("lfo1_on", "LFO1_On", NormalisableRange<float>(0.0f, 1.0f), 1.0f), // true
         std::make_unique<juce::AudioParameterFloat>("lfo1_rate", "LFO1_Rate", NormalisableRange<float>(0.1f, 17.0f), 0.9f), // in Hz
         
         std::make_unique<juce::AudioParameterFloat>("osc2_on", "Osc2_On", NormalisableRange<float>(0.0f, 1.0f), 1.0f), // true
         std::make_unique<juce::AudioParameterFloat>("osc2_waveform", "Osc2_Waveform", NormalisableRange<float>(0.0f, 3.0f), 2.0f),
         std::make_unique<juce::AudioParameterFloat>("osc2_gain", "Osc2_Gain", NormalisableRange<float>(0.0f, 1.0f), 0.5f),
         std::make_unique<juce::AudioParameterFloat>("osc2_octave", "Osc2_Octave", NormalisableRange<float>(-3.0f, 3.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("osc2_semi", "Osc2_Semi", NormalisableRange<float>(-12.0f, 12.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("osc2_detune", "Osc2_Detune", NormalisableRange<float>(-3.0f, 3.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("filter2_on", "Filter2_On", NormalisableRange<float>(0.0f, 1.0f), 0.0f), // false
         std::make_unique<juce::AudioParameterFloat>("filter2_type", "Filter2_Type", NormalisableRange<float>(0.0f, 4.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("filter2_freq", "Filter2_Freq", NormalisableRange<float>(20.0f, 22000.0f), 22000.0f),
         std::make_unique<juce::AudioParameterFloat>("filter2_q", "Filter2_Q", NormalisableRange<float>(0.0f, 1.0f), 0.7f),
         std::make_unique<juce::AudioParameterFloat>("amp2_on", "Amp2_On", NormalisableRange<float>(0.0f, 1.0f), 0.0f), // false
         std::make_unique<juce::AudioParameterFloat>("amp2_pan", "Amp2_Pan", NormalisableRange<float>(-1.0f, 1.0f), 0.0f),
         std::make_unique<juce::AudioParameterFloat>("amp2_level", "Amp2_Level", NormalisableRange<float>(-40.0f, 6.0f), 0.0f), // this is in dB!!!
         std::make_unique<juce::AudioParameterFloat>("lfo2_on", "LFO2_On", NormalisableRange<float>(0.0f, 1.0f), 0.0f), // false
         std::make_unique<juce::AudioParameterFloat>("lfo2_rate", "LFO2_Rate", NormalisableRange<float>(0.1f, 17.0f), 0.9f), // in Hz

         std::make_unique<juce::AudioParameterFloat>("control", "Control", NormalisableRange<float>(0.0f, 2.0f), 1.0f), // default rule-based
         std::make_unique<juce::AudioParameterFloat>("synth", "Synth", NormalisableRange<float>(0.0f, 1.0f), 1.0f) // default synth
                                     })

#endif
{
    rcvr.addListener(this);
    rcvr.connect(1234);
    
    // Control
    controlNumber = 0;
    prevControlNumber = 0;
    
}

BirlAudioProcessor::~BirlAudioProcessor()
{
}

//==============================================================================
void BirlAudioProcessor::oscMessageReceived(const OSCMessage& message){
    if(message.getAddressPattern().toString()=="/birl/keys"){
        int k = 0;
        for (auto* arg = message.begin(); arg != message.end(); ++arg) {
            if (k < NUM_OF_TONEHOLES) {
                birl::toneholes[k] = getArgValue(*arg);
                if (getArgValue(*arg) > birl::maxToneholeArg[k]) {
                    birl::maxToneholeArg[k] = getArgValue(*arg);
                }
                birl::toneholes[k] = (float) getArgValue(*arg) / (birl::maxToneholeArg[k]);
                birl::toneholes[k] = LEAF_clip(0.0, birl::toneholes[k], 1.0);
                ++k;
            }
        }
    }
    if(message.getAddressPattern().toString()=="/birl/buttons") {
        int b = 0;
        for (auto* arg = message.begin(); arg != message.end(); ++arg) {
            if (b < NUM_OF_BUTTONS) {
                if (getArgValue(*arg) == 1)
                    birl::buttons[b] = true;
                else
                    birl::buttons[b] = false;
                ++b;
            }
        }
    }
    if(message.getAddressPattern().toString()=="/birl/breathpos") {
        for (auto* arg = message.begin(); arg != message.end(); ++arg) {
            currentBreathPos = getArgValue(*arg);
            if (currentBreathPos > 0)
                isExhaling = true;
            if(currentBreathPos > maxBreathPos)
                maxBreathPos = currentBreathPos;
            if (isExhaling) {
                birl::breathArray[0] = LEAF_clip(0.0, currentBreathPos / maxBreathPos, 1.0);
                //birl::setGain(LEAF_clip(0.0, currentBreathPos / maxBreathPos, 1.0));
                birl::SFXPhysicalModelSetBreathPressure(birl::breathArray[0]);
                parameters.getParameter("gain")->setValue(birl::breathArray[0]);
            }
        }
    }
    if(message.getAddressPattern().toString()=="/birl/breathneg") {
        for (auto* arg = message.begin(); arg != message.end(); ++arg) {
            currentBreathNeg = getArgValue(*arg);
            if (currentBreathNeg > 0)
                isExhaling = false;
            if(currentBreathNeg > maxBreathNeg)
                maxBreathNeg = currentBreathNeg;
            if (!isExhaling) {
                birl::breathArray[1] = LEAF_clip(0.0, currentBreathNeg / maxBreathNeg, 1.0);
                //birl::setGain(LEAF_clip(0.0, currentBreathNeg / maxBreathNeg, 1.0));
                birl::SFXPhysicalModelSetBreathPressure(birl::breathArray[1]);
                parameters.getParameter("gain")->setValue(birl::breathArray[1]);
            }
        }
    }
    triggerAsyncUpdate();
}

void BirlAudioProcessor::oscBundleReceived(const OSCBundle& bundle) {
    for (auto* element = bundle.begin(); element != bundle.end(); ++element) {
        if (element->isMessage())
            oscMessageReceived(element->getMessage());
        else if (element->isBundle())
            oscBundleReceived(element->getBundle());
    }
    triggerAsyncUpdate();
}

int BirlAudioProcessor::getArgValue(const juce::OSCArgument& arg) {
    if (arg.isInt32()) {
        return arg.getInt32();
    }
    else {
        AlertWindow::showMessageBoxAsync(AlertWindow::WarningIcon,
            "Argument Type Error",
            "Error: Birl input not integer values",
            "OK");
        return 0;
    }
}

float BirlAudioProcessor::getToneholeValue(int index) {
    return birl::toneholes[index];
}
float BirlAudioProcessor::getBreathArrayValue(int index) {
    return birl::breathArray[index];
}

void BirlAudioProcessor::handleAsyncUpdate() {}


const juce::String BirlAudioProcessor::getName() const
{
    return JucePlugin_Name;
}


bool BirlAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool BirlAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool BirlAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double BirlAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int BirlAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int BirlAudioProcessor::getCurrentProgram()
{
    return 0;
}

void BirlAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String BirlAudioProcessor::getProgramName (int index)
{
    return {};
}

void BirlAudioProcessor::changeProgramName (int index, const juce::String& newName)
{
}

//==============================================================================
void BirlAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    LEAF_init(&leaf,sampleRate, birl::medium_memory, 500000, []() {return (float)rand() / RAND_MAX; });
    tMempool_init(&birl::smallPool, birl::small_memory, 80328, &leaf);
    tMempool_init(&birl::largePool, birl::large_memory, 33554432, &leaf);
    birl::initGlobalSFXObjects(leaf);
    
    if (controlNumber == 0) {
        controlNumber = 1;
        loading = true;
    }

}

void BirlAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool BirlAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    juce::ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
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

void BirlAudioProcessor::processBlock (juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    if (loading)
    {
        for (int i = 0; i < birl::ButtonNil; ++i) {
            for (int j = 0; j < birl::ActionNil; ++j) {
                birl::buttonActionsSFX[i][j] = false;
                birl::buttonActionsUI[i][j] = false;
            }
        }
        switch (prevControlNumber) {
            case 1:
                birl::SFXPhysicalModelPMFree();
                break;
            case 2:
                birl::SFXRuleBasedPMFree();
                break;
            case 3:
                birl::SFXRuleBasedSynthFree();
                break;
            case 4:
                birl::SFXNeuralNetPMFree();
                break;
            case 5:
                birl::SFXNeuralNetSynthFree();
                break;
            default:
                break;
        }
        
        switch (controlNumber) {
            case 1:
                birl::SFXPhysicalModelPMAlloc(leaf);
                break;
            case 2:
                birl::SFXRuleBasedPMAlloc();
                break;
            case 3:
                birl::SFXRuleBasedSynthAlloc();
                break;
            case 4:
                birl::SFXNeuralNetPMAlloc();
                break;
            case 5:
                birl::SFXNeuralNetSynthAlloc();
                break;
            default:
                break;
        }
        
        loading = false;
    }
    
    if (!loading)
    {
        switch (controlNumber) {
            case 1:
                birl::SFXPhysicalModelPMFrame(buffer);
                break;
            case 2:
                birl::SFXRuleBasedPMFrame();
                break;
            case 3:
                birl::SFXRuleBasedSynthFrame();
                break;
            case 4:
                birl::SFXNeuralNetPMFrame();
                break;
            case 5:
                birl::SFXNeuralNetSynthFrame();
                break;
            default:
                break;
        }
    }
        
 
    /*
     MidiMessage m;
     for (MidiMessageMetadata metadata : midiMessages) {
        m = metadata.getMessage();
        handleIncomingMidiMessage(nullptr, m);
     }
     */
    

    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());
    
    birl::controlKnobValues[0][0] = parameters.getParameter("gain")->getValue();
    birl::controlKnobValues[0][1] = parameters.getParameter("fundamental")->getValue();
    birl::controlKnobValues[0][2] = parameters.getParameter("num_toneholes")->getValue();
    birl::controlKnobValues[0][3] = parameters.getParameter("dcblocker1")->getValue();
    birl::controlKnobValues[0][4] = parameters.getParameter("dcblocker2")->getValue();
    birl::controlKnobValues[0][5] = parameters.getParameter("biquad_coeff1")->getValue();
    birl::controlKnobValues[0][6] = parameters.getParameter("biquad_coeff2")->getValue();
    birl::controlKnobValues[0][7] = parameters.getParameter("biquad_coeff3")->getValue();
    birl::controlKnobValues[0][8] = parameters.getParameter("biquad_coeff4")->getValue();
    birl::controlKnobValues[0][9] = parameters.getParameter("biquad_coeff5")->getValue();
    birl::controlKnobValues[0][10] = parameters.getParameter("pf1_cutoff")->getValue();
    birl::controlKnobValues[0][11] = parameters.getParameter("pf1_q")->getValue();
    birl::controlKnobValues[0][12] = parameters.getParameter("pf2_cutoff")->getValue();
    birl::controlKnobValues[0][13] = parameters.getParameter("pf2_q")->getValue();
    birl::controlKnobValues[0][14] = parameters.getParameter("lp1_cutoff")->getValue();
    birl::controlKnobValues[0][15] = parameters.getParameter("lp1_q")->getValue();
    birl::controlKnobValues[0][16] = parameters.getParameter("lp2_cutoff")->getValue();
    birl::controlKnobValues[0][17] = parameters.getParameter("lp2_q")->getValue();
    birl::controlKnobValues[0][18] = parameters.getParameter("shaper_drive")->getValue();
    birl::controlKnobValues[0][19] = parameters.getParameter("shaper_mix")->getValue();
    birl::controlKnobValues[0][20] = parameters.getParameter("noise_gain")->getValue();
    birl::controlKnobValues[0][21] = parameters.getParameter("noise_bp_cutoff")->getValue();
    birl::controlKnobValues[0][22] = parameters.getParameter("noise_bp_q")->getValue();

    birl::controlKnobValues[2][0] = parameters.getParameter("osc1_on")->getValue();
    birl::controlKnobValues[2][1] = parameters.getParameter("osc1_waveform")->getValue();
    birl::controlKnobValues[2][2] = parameters.getParameter("osc1_gain")->getValue();
    birl::controlKnobValues[2][3] = parameters.getParameter("osc1_octave")->getValue();
    birl::controlKnobValues[2][4] = parameters.getParameter("osc1_semi")->getValue();
    birl::controlKnobValues[2][5] = parameters.getParameter("osc1_detune")->getValue();
    birl::controlKnobValues[2][6] = parameters.getParameter("filter1_on")->getValue();
    birl::controlKnobValues[2][7] = parameters.getParameter("filter1_type")->getValue();
    birl::controlKnobValues[2][8] = parameters.getParameter("filter1_tof2")->getValue();
    birl::controlKnobValues[2][9] = parameters.getParameter("filter1_freq")->getValue();
    birl::controlKnobValues[2][10] = parameters.getParameter("filter1_q")->getValue();
    birl::controlKnobValues[2][11] = parameters.getParameter("amp1_on")->getValue();
    birl::controlKnobValues[2][12] = parameters.getParameter("amp1_pan")->getValue();
    birl::controlKnobValues[2][13] = parameters.getParameter("amp1_level")->getValue();
    birl::controlKnobValues[2][14] = parameters.getParameter("lfo1_on")->getValue();
    birl::controlKnobValues[2][15] = parameters.getParameter("lfo1_rate")->getValue();
    
    birl::controlKnobValues[2][16] = parameters.getParameter("osc2_on")->getValue();
    birl::controlKnobValues[2][17] = parameters.getParameter("osc2_waveform")->getValue();
    birl::controlKnobValues[2][18] = parameters.getParameter("osc2_gain")->getValue();
    birl::controlKnobValues[2][19] = parameters.getParameter("osc2_octave")->getValue();
    birl::controlKnobValues[2][20] = parameters.getParameter("osc2_semi")->getValue();
    birl::controlKnobValues[2][21] = parameters.getParameter("osc2_detune")->getValue();
    birl::controlKnobValues[2][22] = parameters.getParameter("filter2_on")->getValue();
    birl::controlKnobValues[2][23] = parameters.getParameter("filter2_type")->getValue();
    birl::controlKnobValues[2][24] = parameters.getParameter("filter2_freq")->getValue();
    birl::controlKnobValues[2][25] = parameters.getParameter("filter2_q")->getValue();
    birl::controlKnobValues[2][26] = parameters.getParameter("amp2_on")->getValue();
    birl::controlKnobValues[2][27] = parameters.getParameter("amp2_pan")->getValue();
    birl::controlKnobValues[2][28] = parameters.getParameter("amp2_level")->getValue();
    birl::controlKnobValues[2][29] = parameters.getParameter("lfo2_on")->getValue();
    birl::controlKnobValues[2][30] = parameters.getParameter("lfo2_rate")->getValue();
    
    
    // for (int channel = 0; channel < buffer.getNumChannels(); channel++) {
        // output data
        float* leftChannel = buffer.getWritePointer(0);
        float* rightChannel;
        if (buffer.getNumChannels() == 1) {
            rightChannel = buffer.getWritePointer(0);
        }
        else {
            rightChannel = buffer.getWritePointer(1);
        }
        
        for (int i = 0; i < buffer.getNumSamples(); i++)
        {
            audio[0] = leftChannel[i];
            audio[1] = rightChannel[i];

            if (!loading)
                switch (controlNumber) {
                    case 1:
                        birl::SFXPhysicalModelPMTick(audio);
                        break;
                    case 2:
                        birl::SFXRuleBasedPMTick(audio);
                        break;
                    case 3:
                        birl::SFXRuleBasedSynthTick(audio);
                        break;
                    case 4:
                        birl::SFXNeuralNetPMTick(audio);
                        break;
                    case 5:
                         birl::SFXNeuralNetSynthTick(audio);
                        break;
                    default:
                        break;
                }
            
           // audio[0] = LEAF_interpolation_linear(leftChannel[i], audio[0], interpVal);
            buffer.setSample(0, i, audio[0]);
            if (buffer.getNumChannels() > 1) {
                //audio[1] = LEAF_interpolation_linear(rightChannel[i], audio[i], interpVal);
                buffer.setSample(1, i, audio[1]);

            }
    //        for (int j = 0; j < totalNumInputChannels; ++j) {
    //            channelData.getUnchecked(j)[i] = audioBufferTick;
    //    }
    //
    //    for (int channel = 0; channel < totalNumInputChannels; ++channel)
    //    {
    //        auto* channelData = buffer.getWritePointer (channel);
        }
    melatonin::printSparkline(buffer);

        
    // }
    // melatonin::printSparkline(buffer);
}

//==============================================================================
bool BirlAudioProcessor::hasEditor() const
{
    return true;
}

juce::AudioProcessorEditor* BirlAudioProcessor::createEditor()
{
    return new BirlAudioProcessorEditor (*this);
}

//==============================================================================
void BirlAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    std::unique_ptr<XmlElement> xml(new XmlElement("Birl"));
    switch (controlNumber) {
        case 1:
//            *pmpm_gain = birl::controlKnobValues[controlNumber - 1][0];
//            juce::MemoryOutputStream(destData, true).writeFloat(*pmpm_gain);
            xml->setAttribute("Gain", birl::controlKnobValues[controlNumber - 1][0]);
            xml->setAttribute("Fundamental", birl::controlKnobValues[controlNumber - 1][1]);
            xml->setAttribute("NumToneholes", birl::controlKnobValues[controlNumber - 1][2]);
            xml->setAttribute("DCBlocker1", birl::controlKnobValues[controlNumber - 1][3]);
            xml->setAttribute("DCBlocker2", birl::controlKnobValues[controlNumber - 1][4]);
            xml->setAttribute("BiquadCoeff1", birl::controlKnobValues[controlNumber - 1][5]);
            xml->setAttribute("BiquadCoeff2", birl::controlKnobValues[controlNumber - 1][6]);
            xml->setAttribute("BiquadCoeff3", birl::controlKnobValues[controlNumber - 1][7]);
            xml->setAttribute("BiquadCoeff4", birl::controlKnobValues[controlNumber - 1][8]);
            xml->setAttribute("BiquadCoeff5", birl::controlKnobValues[controlNumber - 1][9]);
            xml->setAttribute("PF1Cutoff", birl::controlKnobValues[controlNumber - 1][10]);
            xml->setAttribute("PF1Q", birl::controlKnobValues[controlNumber - 1][11]);
            xml->setAttribute("PF2Cutoff", birl::controlKnobValues[controlNumber - 1][12]);
            xml->setAttribute("PF2Q", birl::controlKnobValues[controlNumber - 1][13]);
            xml->setAttribute("LP1Cutoff", birl::controlKnobValues[controlNumber - 1][14]);
            xml->setAttribute("LP1Q", birl::controlKnobValues[controlNumber - 1][15]);
            xml->setAttribute("LP2Cutoff", birl::controlKnobValues[controlNumber - 1][16]);
            xml->setAttribute("LP2Q", birl::controlKnobValues[controlNumber - 1][17]);
            xml->setAttribute("ShaperDrive", birl::controlKnobValues[controlNumber - 1][18]);
            xml->setAttribute("ShaperMix", birl::controlKnobValues[controlNumber - 1][19]);
            xml->setAttribute("NoiseGain", birl::controlKnobValues[controlNumber - 1][20]);
            xml->setAttribute("NoiseBPCutoff", birl::controlKnobValues[controlNumber - 1][21]);
            xml->setAttribute("NoiseBPQ", birl::controlKnobValues[controlNumber - 1][22]);
            break;
        case 2:
            break;
        case 3:
            xml->setAttribute("Osc1_On", birl::controlKnobValues[controlNumber - 1][0]);
            xml->setAttribute("Osc1_Waveform", birl::controlKnobValues[controlNumber - 1][1]);
            xml->setAttribute("Osc1_Gain", birl::controlKnobValues[controlNumber - 1][2]);
            xml->setAttribute("Osc1_Octave", birl::controlKnobValues[controlNumber - 1][3]);
            xml->setAttribute("Osc1_Semi", birl::controlKnobValues[controlNumber - 1][4]);
            xml->setAttribute("Osc1_Detune", birl::controlKnobValues[controlNumber - 1][5]);
            xml->setAttribute("Filter1_On", birl::controlKnobValues[controlNumber - 1][6]);
            xml->setAttribute("Filter1_Type", birl::controlKnobValues[controlNumber - 1][7]);
            xml->setAttribute("Filter1_toF2", birl::controlKnobValues[controlNumber - 1][8]);
            xml->setAttribute("Filter1_Freq", birl::controlKnobValues[controlNumber - 1][9]);
            xml->setAttribute("Filter1_Res", birl::controlKnobValues[controlNumber - 1][10]);
            xml->setAttribute("Amp1_On", birl::controlKnobValues[controlNumber - 1][11]);
            xml->setAttribute("Amp1_Pan", birl::controlKnobValues[controlNumber - 1][12]);
            xml->setAttribute("Amp1_Level", birl::controlKnobValues[controlNumber - 1][13]);
            xml->setAttribute("LFO1_On", birl::controlKnobValues[controlNumber - 1][14]);
            xml->setAttribute("LFO1_Rate", birl::controlKnobValues[controlNumber - 1][15]);
        
            xml->setAttribute("Osc2_On", birl::controlKnobValues[controlNumber - 1][16]);
            xml->setAttribute("Osc2_Waveform", birl::controlKnobValues[controlNumber - 1][17]);
            xml->setAttribute("Osc2_Gain", birl::controlKnobValues[controlNumber - 1][18]);
            xml->setAttribute("Osc2_Octave", birl::controlKnobValues[controlNumber - 1][19]);
            xml->setAttribute("Osc2_Semi", birl::controlKnobValues[controlNumber - 1][20]);
            xml->setAttribute("Osc2_Detune", birl::controlKnobValues[controlNumber - 1][21]);
            xml->setAttribute("Filter2_On", birl::controlKnobValues[controlNumber - 1][22]);
            xml->setAttribute("Filter2_Type", birl::controlKnobValues[controlNumber - 1][23]);
            xml->setAttribute("Filter2_Freq", birl::controlKnobValues[controlNumber - 1][24]);
            xml->setAttribute("Filter2_Res", birl::controlKnobValues[controlNumber - 1][25]);
            xml->setAttribute("Amp2_On", birl::controlKnobValues[controlNumber - 1][26]);
            xml->setAttribute("Amp2_Pan", birl::controlKnobValues[controlNumber - 1][27]);
            xml->setAttribute("Amp2_Level", birl::controlKnobValues[controlNumber - 1][28]);
            xml->setAttribute("LFO2_On", birl::controlKnobValues[controlNumber - 1][29]);
            xml->setAttribute("LFO2_Rate", birl::controlKnobValues[controlNumber - 1][30]);
            break;
        case 4:
            break;
        case 5:
            break;
        default:
            break;
    }
    copyXmlToBinary(*xml, destData);
}

void BirlAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    std::unique_ptr<XmlElement> xml(getXmlFromBinary(data, sizeInBytes));
    
    if (xml.get() != nullptr && xml->hasTagName("Birl")) {
        switch (controlNumber) {
            case 1:
                birl::controlKnobValues[controlNumber - 1][0] = xml->getDoubleAttribute("Gain");
                birl::controlKnobValues[controlNumber - 1][1] = xml->getDoubleAttribute("Fundamental");
                birl::controlKnobValues[controlNumber - 1][2] = xml->getDoubleAttribute("NumToneholes");
                birl::controlKnobValues[controlNumber - 1][3] = xml->getDoubleAttribute("DCBlocker1");
                birl::controlKnobValues[controlNumber - 1][4] = xml->getDoubleAttribute("DCBlocker2");
                birl::controlKnobValues[controlNumber - 1][5] = xml->getDoubleAttribute("BiquadCoeff1");
                birl::controlKnobValues[controlNumber - 1][6] = xml->getDoubleAttribute("BiquadCoeff2");
                birl::controlKnobValues[controlNumber - 1][7] = xml->getDoubleAttribute("BiquadCoeff3");
                birl::controlKnobValues[controlNumber - 1][8] = xml->getDoubleAttribute("BiquadCoeff4");
                birl::controlKnobValues[controlNumber - 1][9] = xml->getDoubleAttribute("BiquadCoeff5");
                birl::controlKnobValues[controlNumber - 1][10] = xml->getDoubleAttribute("PF1Cutoff");
                birl::controlKnobValues[controlNumber - 1][11] = xml->getDoubleAttribute("PF1Q");
                birl::controlKnobValues[controlNumber - 1][12] = xml->getDoubleAttribute("PF2Cutoff");
                birl::controlKnobValues[controlNumber - 1][13] = xml->getDoubleAttribute("PF2Q");
                birl::controlKnobValues[controlNumber - 1][14] = xml->getDoubleAttribute("LP1Cutoff");
                birl::controlKnobValues[controlNumber - 1][15] = xml->getDoubleAttribute("LP1Q");
                birl::controlKnobValues[controlNumber - 1][16] = xml->getDoubleAttribute("LP2Cutoff");
                birl::controlKnobValues[controlNumber - 1][17] = xml->getDoubleAttribute("LP2Q");
                birl::controlKnobValues[controlNumber - 1][18] = xml->getDoubleAttribute("ShaperDrive");
                birl::controlKnobValues[controlNumber - 1][19] = xml->getDoubleAttribute("ShaperMix");
                birl::controlKnobValues[controlNumber - 1][20] = xml->getDoubleAttribute("NoiseGain");
                birl::controlKnobValues[controlNumber - 1][21] = xml->getDoubleAttribute("NoiseBPCutoff");
                birl::controlKnobValues[controlNumber - 1][22] = xml->getDoubleAttribute("NoiseBPQ");
                break;
            case 2:
                break;
            case 3:
                birl::controlKnobValues[controlNumber - 1][0] = xml->getDoubleAttribute("Osc1_On");
                birl::controlKnobValues[controlNumber - 1][1] = xml->getDoubleAttribute("Osc1_Waveform");
                birl::controlKnobValues[controlNumber - 1][2] = xml->getDoubleAttribute("Osc1_Gain");
                birl::controlKnobValues[controlNumber - 1][3] = xml->getDoubleAttribute("Osc1_Octave");
                birl::controlKnobValues[controlNumber - 1][4] = xml->getDoubleAttribute("Osc1_Semi");
                birl::controlKnobValues[controlNumber - 1][5] = xml->getDoubleAttribute("Osc1_Detune");
                birl::controlKnobValues[controlNumber - 1][6] = xml->getDoubleAttribute("Filter1_On");
                birl::controlKnobValues[controlNumber - 1][7] = xml->getDoubleAttribute("Filter1_Type");
                birl::controlKnobValues[controlNumber - 1][8] = xml->getDoubleAttribute("Filter1_toF2");
                birl::controlKnobValues[controlNumber - 1][9] = xml->getDoubleAttribute("Filter1_Freq");
                birl::controlKnobValues[controlNumber - 1][10] = xml->getDoubleAttribute("Filter1_Res");
                birl::controlKnobValues[controlNumber - 1][11] = xml->getDoubleAttribute("Amp1_On");
                birl::controlKnobValues[controlNumber - 1][12] = xml->getDoubleAttribute("Amp1_Pan");
                birl::controlKnobValues[controlNumber - 1][13] = xml->getDoubleAttribute("Amp1_Level");
                birl::controlKnobValues[controlNumber - 1][14] = xml->getDoubleAttribute("LFO1_On");
                birl::controlKnobValues[controlNumber - 1][15] = xml->getDoubleAttribute("LF01_Rate");
                birl::controlKnobValues[controlNumber - 1][16] = xml->getDoubleAttribute("Osc2_On");
                birl::controlKnobValues[controlNumber - 1][17] = xml->getDoubleAttribute("Osc2_Waveform");
                birl::controlKnobValues[controlNumber - 1][18] = xml->getDoubleAttribute("Osc2_Gain");
                birl::controlKnobValues[controlNumber - 1][19] = xml->getDoubleAttribute("Osc2_Octave");
                birl::controlKnobValues[controlNumber - 1][20] = xml->getDoubleAttribute("Osc2_Semi");
                birl::controlKnobValues[controlNumber - 1][21] = xml->getDoubleAttribute("Osc2_Detune");
                birl::controlKnobValues[controlNumber - 1][22] = xml->getDoubleAttribute("Filter2_On");
                birl::controlKnobValues[controlNumber - 1][23] = xml->getDoubleAttribute("Filter2_Type");
                birl::controlKnobValues[controlNumber - 1][24] = xml->getDoubleAttribute("Filter2_Freq");
                birl::controlKnobValues[controlNumber - 1][25] = xml->getDoubleAttribute("Filter2_Res");
                birl::controlKnobValues[controlNumber - 1][26] = xml->getDoubleAttribute("Amp2_On");
                birl::controlKnobValues[controlNumber - 1][27] = xml->getDoubleAttribute("Amp2_Pan");
                birl::controlKnobValues[controlNumber - 1][28] = xml->getDoubleAttribute("Amp2_Level");
                birl::controlKnobValues[controlNumber - 1][29] = xml->getDoubleAttribute("LFO2_On");
                birl::controlKnobValues[controlNumber - 1][30] = xml->getDoubleAttribute("LF02_Rate");
                break;
            case 4:
                break;
            case 5:
                break;
            default:
                break;
        }
        interpVal = xml->getDoubleAttribute("DryWet");
        birl::buttonActionsSFX[0][birl::ActionPress] = xml->getDoubleAttribute("ButtonOctaveUp");
        birl::buttonActionsSFX[1][birl::ActionPress] = xml->getDoubleAttribute("ButtonOctaveDown");
        birl::buttonActionsSFX[2][birl::ActionPress] = xml->getDoubleAttribute("ButtonPrevTuning");
        birl::buttonActionsSFX[3][birl::ActionPress] = xml->getDoubleAttribute("ButtonNextTuning");
        birl::buttonActionsSFX[4][birl::ActionPress] = xml->getDoubleAttribute("ButtonA");
        birl::buttonActionsSFX[5][birl::ActionPress] = xml->getDoubleAttribute("ButtonB");
        birl::buttonActionsSFX[6][birl::ActionPress] = xml->getDoubleAttribute("ButtonC");
        birl::buttonActionsSFX[7][birl::ActionPress] = xml->getDoubleAttribute("ButtonD");
        birl::buttonActionsSFX[8][birl::ActionPress] = xml->getDoubleAttribute("ButtonE");
    }
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new BirlAudioProcessor();
}

void BirlAudioProcessor::handleIncomingMidiMessage(MidiInput* source, const MidiMessage& message) {
    const ScopedValueSetter<bool> scopedInputFlag(isAddingFromMidiInput, true);
    if (message.isNoteOn()) {
        birl::noteOn(message.getNoteNumber(), message.getVelocity());
    } else {
        birl::noteOff(message.getNoteNumber(), message.getVelocity());
    }
    keyboardState.processNextMidiEvent(message);
}

void BirlAudioProcessor::handleNoteOn(MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) {
    if (!isAddingFromMidiInput) {
        auto m = MidiMessage::noteOn(midiChannel, midiNoteNumber, velocity);
        m.setTimeStamp(Time::getMillisecondCounterHiRes() * 0.001);
    }
}

void BirlAudioProcessor::handleNoteOff(MidiKeyboardState*, int midiChannel, int midiNoteNumber, float velocity) {
    if (!isAddingFromMidiInput) {
        auto m = MidiMessage::noteOff(midiChannel, midiNoteNumber);
        m.setTimeStamp(Time::getMillisecondCounterHiRes() * 0.001);
    }
}
