
#ifndef MAINCOMPONENT_H_INCLUDED
#define MAINCOMPONENT_H_INCLUDED
#include "../JuceLibraryCode/JuceHeader.h"
#include "AudioSettings.h"
#include "Synth.h"
// =====================================================================================================
class MainContentComponent : public Component,
                             public Slider::Listener

{
private:
	MidiKeyboardState keyboardState;
	MidiKeyboardComponent keyboardComponent;
	SharedResourcePointer<AudioSettings> sharedAudioSettings;
	AudioDeviceManager& inputDeviceManager =    sharedAudioSettings->getInputDeviceManager();
	AudioDeviceManager& outputDeviceManager =   sharedAudioSettings->getOutputDeviceManager();
	PercSynthPlayer synthPlayer;
    Slider dampingSlider,stiffnessSlider,thickSlider, volumeSlider;
    Label dampingLabel,stiffnessLabel,thickLabel,volumeLabel;
    ComboBox materials;

public:
    MainContentComponent()
		: synthPlayer(keyboardState, inputDeviceManager, outputDeviceManager)
        , keyboardComponent(keyboardState, MidiKeyboardComponent::horizontalKeyboard)
	{
		setSize(900, 300);
		//setOpaque(true);
        
        addAndMakeVisible (keyboardComponent);
        keyboardComponent.setAvailableRange(21, 108);
        
        addAndMakeVisible (volumeSlider);
        addAndMakeVisible (volumeLabel);
        volumeLabel.setText ("MASTER VOLUME", dontSendNotification);
        volumeLabel.attachToComponent (&volumeSlider, false);
        volumeSlider.setRange (0.0, 100.0);
        volumeSlider.setValue(50.0);
        volumeSlider.setComponentID ("volume");
        volumeSlider.setSliderStyle(juce::Slider::LinearVertical);
        volumeSlider.addListener (&synthPlayer);
        
        addAndMakeVisible( materials);
        materials.addItem("Bright metallic",1);
        materials.addItem("Dark woody",2);
        materials.setText("Select plate 'material'");
        materials.addListener(&synthPlayer);
        materials.setComponentID("materials");
        
        addAndMakeVisible (dampingSlider);
        addAndMakeVisible (dampingLabel);
        dampingLabel.setText ("Damping", dontSendNotification);
        dampingLabel.attachToComponent (&dampingSlider, true);
        dampingSlider.setRange (0.0, 1.0);
        dampingSlider.setComponentID ("damping");
        dampingSlider.addListener (&synthPlayer);

        addAndMakeVisible (stiffnessSlider);
        addAndMakeVisible (stiffnessLabel);
        stiffnessLabel.setText ("Stiffness", dontSendNotification);
        stiffnessLabel.attachToComponent (&stiffnessSlider, true);
        stiffnessSlider.setRange (0.0, 1.0);
        stiffnessSlider.setComponentID ("stiffness");  // to recognise your slider in the callback
        stiffnessSlider.addListener (&synthPlayer);
        
        addAndMakeVisible (thickSlider);
        addAndMakeVisible (thickLabel);
        thickLabel.setText ("Thickness", dontSendNotification);
        thickLabel.attachToComponent (&thickSlider, true);
        thickSlider.setRange (0.05, 0.1);
        thickSlider.setComponentID ("thickness");  // to recognise your slider in the callback
        thickSlider.addListener (&synthPlayer);
    }
    
    void resized() override
    {   const int sliderLeft = 120;
        volumeSlider.setBounds(750,40,getWidth()/3 - sliderLeft - 50,120);
        dampingSlider.setBounds (sliderLeft, 40, getWidth()/3 - sliderLeft - 10, 40);
        stiffnessSlider.setBounds (sliderLeft, 60, getWidth()/3 - sliderLeft - 10, 60);
        thickSlider.setBounds (4*sliderLeft, 40, getWidth()/3 + sliderLeft - 300, 40);
        materials.setBounds (4*sliderLeft, 90, 1.5*sliderLeft, 50);
        
        int w = (int)keyboardComponent.getKeyWidth() * (7 * 7 + 3), h = 80;
        keyboardComponent.setSize(w, h);
        keyboardComponent.setCentrePosition(getWidth() / 2, getHeight() - getHeight() / 6);
    }
    
    void sliderValueChanged (Slider *slider) override {}

	void paint(Graphics& g) override
    {   g.fillAll(Colours::slategrey);
        g.setFont(18.0f);
        g.drawText("Press on the keys to strike a plate!", 50, 90, 400, 120, Justification::centred, true);
        g.setFont(11.0f);
        g.drawText("lower octaves~bigger plates // higher octaves~smaller plates", 0, 110, 500, 120, Justification::centred, true);
	}

private:

	JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainContentComponent)
};


// This is called by the app startup code to create the main component instance
Component* createMainContentComponent()  {
	return new MainContentComponent(); 
}


#endif  // MAINCOMPONENT_H_INCLUDED
