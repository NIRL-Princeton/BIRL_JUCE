//#ifndef RBI_H
//#define RBI_H
//
//#include <juce_core/juce_core.h>
//#include "PluginProcessor.h"
//using namespace juce;
//#define BIRL_PORT 1234
//#define OPENHOLE_MIDI 53
//#define CLIP_RATIO 0.7
//
//class RuleBasedInterface:
//            public Component,
//            private OSCReceiver::Listener<OSCReceiver::MessageLoopCallback>,
//            private AsyncUpdater,
//            public Slider::Listener
//{
//public:
//    RuleBasedInterface(BirlAudioProcessor&);
//    ~RuleBasedInterface();
//    void paint (Graphics& g) override;
//    void resized() override;
//    //=======================
//    void populateKeySliders();
//    void populateButtonButtons();
//    //=========================
//    void oscMessageReceived (const OSCMessage& message) override;
//    void oscBundleReceived (const OSCBundle& bundle) override;
//    //==========================
//    int getArgValue(const OSCArgument& arg);
//    void setMaxBreathPosValue(double);
//    void setMaxBreathNegValue(double);
//    void updateMaxBreathLabel(Label*, double);
//    //==========================
//    void handleAsyncUpdate() override;
//    //==========================
//    void connect();
//    void disconnect();
//    bool isConnected() const;
//    void connectButtonClicked();
//    void handleArgTypeError();
//    //===========================
//    void adjustMidi();
//    void calculateFrequency(double);
//    void octaveButtonPressed();
//    void transpose();
//    int getScaleType();
//    //===========================
//    void sliderValueChanged(Slider* slider) override;
//private:
//
//        Slider monday;
//
//        BirlAudioProcessor& processor;
//        OSCReceiver oscReceiver{};
//        int currentPortNumber{-1};
//        TextButton connectButton { "Connect" };
//
//
//        Slider      keySliders[9];
//        TextButton  buttonButtons[9];
//
//        Slider  breathpos{};
//        Slider  breathneg{};
//        Label   maxbreathposlabel{};
//        Label   maxbreathneglabel{};
//
//        double  absbreath{};
//        double  maxbreathpos{1.0};
//        double  maxbreathneg{1.0};
//
//        double  currentbreathpos{};
//        double  currentbreathneg{};
//
//        bool    isExhaling{true};
//
//        float   keys[9];
//        int     maxKeyArg[9]{1, 1, 1, 1, 1, 1, 1, 1, 1};
//
//        float   buttons[9];
//        bool    upIsLocked{true};
//        bool    downIsLocked{true};
//
//        double  octaveTransposition{0.0};
//        double  currentMidi{OPENHOLE_MIDI};
//
//        int tempScaleType{3};
//
//        double normalscale[13]
//                = {0.0, 1.117313, 2.039101, 3.156414, 3.86314, 4.980453, 5.902237,
//                    7.019547, 8.136864, 8.84359, 9.960899, 10.88269, 12.0};
//        double otonalscale[13]
//                = {0.0, 1.049553, 2.039101, 2.97128, 3.863136, 4.707809, 5.513176,
//                    7.019547, 8.405273, 9.058647, 9.688255, 10.88269, 12.0};
//        double utonalscale[13]
//                = {0.0, 1.1131, 2.311745, 3.15641, 3.594727, 4.980453, 6.486824,
//                    7.292191, 8.136864, 9.024847, 9.960899, 10.950447, 12.0};
//        double equaltemperedscale[13]
//                = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
//                    8.0, 9.0, 10.0, 11.0, 12.0};
//        double meantonescale[13]
//                = {0.0, 0.76, 1.932, 3.103, 3.863, 5.034, 5.795, 6.966, 7.726,
//                    8.897, 10.068, 10.829, 12.0};
//
//
//};
//
//#endif /*RBI_H*/
