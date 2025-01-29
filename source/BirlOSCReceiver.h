#ifndef BIRLOSCRECEIVER_H
#define BIRLOSCRECEIVER_H

#include <juce_osc/juce_osc.h>
#include "PluginProcessor.h"
#include "PhysicalModelEditor.h"

class BirlOSCReceiver:
            private OSCReceiver::Listener<OSCReceiver::MessageLoopCallback>,
            private AsyncUpdater,
            public PhysicalModelEditor
{
public:
    BirlOSCReceiver(BirlOneHoleAudioProcessor&);
    ~BirlOSCReceiver();
    void oscMessageReceived (const OSCMessage&) override;
    void oscBundleReceived (const OSCBundle&) override;
    
    int getArgValue(const OSCArgument&);
    void handleAsyncUpdate() override;
    void connectButtonClicked();

    void handleConnectError();
    void handleArgTypeError();

protected:
    BirlOneHoleAudioProcessor& processor;
    OSCReceiver rcvr;
    
    int port;
    double absBreath_;
    double maxBreathPos_{1};
    double maxBreathNeg_{1};
    double currentBreathPos_;
    double currentBreathNeg_;
    bool isExhaling_;
};

#endif
