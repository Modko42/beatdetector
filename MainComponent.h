
#include <vector>
#include <chrono>
#include <thread>
#include <JuceHeader.h>

using namespace juce;
using namespace std;


#pragma once

function<void()> drumhit;

void wait(int mseconds){
    this_thread::sleep_for(std::chrono::milliseconds(mseconds));
}
//==============================================================================
class SimpleThumbnailComponent : public juce::Component,
                                 private juce::ChangeListener
{
public:
    SimpleThumbnailComponent (int sourceSamplesPerThumbnailSample,
                              juce::AudioFormatManager& formatManager,
                              juce::AudioThumbnailCache& cache,
                              AudioTransportSource& transportSourceToUse)
       : thumbnail (sourceSamplesPerThumbnailSample, formatManager, cache),
            transportSource (transportSourceToUse)
    {
        thumbnail.addChangeListener (this);
    }

    void setFile (const juce::File& file)
    {
        thumbnail.setSource (new juce::FileInputSource (file));
    }

    void paint (juce::Graphics& g) override
    {
        if (thumbnail.getNumChannels() == 0)
            paintIfNoFileLoaded (g);
        else
            paintIfFileLoaded (g);
    }

    void paintIfNoFileLoaded (juce::Graphics& g)
    {
        g.fillAll (juce::Colours::black);
        g.setColour (juce::Colours::darkgrey);
        g.drawFittedText ("No File Loaded", getLocalBounds(), juce::Justification::centred, 1);
    }

    void paintIfFileLoaded (juce::Graphics& g)
    {
        g.fillAll (juce::Colours::black);
        
        g.setColour (juce::Colour(0xff0ffA600));
        
        counter = transportSource.getCurrentPosition() / (double)interval_time;
        thumbnail.drawChannels (g, getLocalBounds(), 0 + counter*interval_time, interval_time + counter*interval_time, 1.0f);
    }

    void changeListenerCallback (juce::ChangeBroadcaster* source) override
    {
        if (source == &thumbnail)
            thumbnailChanged();
    }

    void setInterval(int interval_){
        interval_time = interval_;
    }
    
private:
    void thumbnailChanged()
    {
        repaint();
    }

    juce::AudioTransportSource& transportSource;
    juce::AudioThumbnail thumbnail;
    int interval_time = 5;
    int counter = 0;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SimpleThumbnailComponent)
};

//------------------------------------------------------------------------------

class SimplePositionOverlay : public juce::Component,
                              private juce::Timer
{
public:
    SimplePositionOverlay (juce::AudioTransportSource& transportSourceToUse)
       : transportSource (transportSourceToUse)
    {
        startTimer (40);
    }

    void paint (juce::Graphics& g) override
    {
        auto duration = (float) transportSource.getLengthInSeconds();

        if (duration > 0.0)
        {
            auto audioPosition = (float) transportSource.getCurrentPosition();

            //auto drawPosition = (audioPosition / duration) * (float) getWidth();
            float drawPosition = ((float)((int)(audioPosition*100.0f) % (interval_time*100))) /
                (float)interval_time / 100.0f * getWidth();
            g.setColour (juce::Colour(0xffC900ff));
            g.drawLine (drawPosition, 0.0f, drawPosition, (float) getHeight(), 2.0f);
        }
    }

    void mouseDown (const juce::MouseEvent& event) override
    {
        auto duration = transportSource.getLengthInSeconds();

        if (duration > 0.0)
        {
            auto clickPosition = event.position.x;
            auto audioPosition = (clickPosition / (float) getWidth()) * duration;

            transportSource.setPosition (audioPosition);
        }
    }
    
    void setInterval(int interval_){
        interval_time = interval_;
    }
    

private:
    void timerCallback() override
    {
        repaint();
    }

    int interval_time = 5;
    int counter = 0;
    juce::AudioTransportSource& transportSource;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SimplePositionOverlay)
};


class SimplePeakOverlay : public juce::Component,private KeyListener,private Timer
{
public:
    SimplePeakOverlay (juce::AudioTransportSource& transportSourceToUse)
       : transportSource (transportSourceToUse)
    {
        setWantsKeyboardFocus(true);
        addKeyListener(this);
        startTimer(40);
    }

    
    
    void paint (juce::Graphics& g) override
    {

        auto duration = (float) transportSource.getLengthInSeconds();
        if(duration > 0.0)
        {
            if(transportSource.getCurrentPosition() == 0)
                counter = 0;
            else
                counter = transportSource.getCurrentPosition() / (double)interval_time;
            int sampleRate = numberofSamples / duration;
            g.setColour (juce::Colour(0xff000fC4));
            for (float i : calculated_peaks) {
                if(i>counter*sampleRate*interval_time && i<(counter+1)*sampleRate*interval_time && numberofSamples>i){
                    float one_interval = (float)sampleRate*interval_time;
                    float percentage_i = (i-counter*sampleRate*interval_time) / one_interval;
                    g.drawLine (percentage_i*getWidth(), 0.0f, percentage_i*getWidth(), (float)getHeight()/1.3f, 1.0f);
                }
            }
            
            
            for (float i : drum_hits) {
                g.setColour (getColorforHit(i, calculated_peaks));
                if(i>counter*sampleRate*interval_time && i<(counter+1)*sampleRate*interval_time && numberofSamples>i){
                    float one_interval = (float)sampleRate*interval_time;
                    float percentage_i = (i-counter*sampleRate*interval_time) / one_interval;
                    g.drawLine (percentage_i*getWidth(), 0.0f, percentage_i*getWidth(), (float)getHeight()/1.3f, 1.0f);
                }
            }
            
            g.setColour (juce::Colour(0xff000fC4));
            g.drawLine(0.0f,(float)getHeight()/1.2f,0.0f,(float)getHeight(),7.0f);
            g.drawLine(getWidth(),(float)getHeight()/1.2f,getWidth(),(float)getHeight(),7.0f);
            g.drawLine(getWidth()/4,(float)getHeight()/1.2f,getWidth()/4,(float)getHeight(),5.0f);
            g.drawLine(getWidth()/4*3,(float)getHeight()/1.2f,getWidth()/4*3,(float)getHeight(),5.0f);
            g.drawLine(getWidth()/2,(float)getHeight()/1.2f,getWidth()/2,(float)getHeight(),5.0f);
            
            
            g.setColour (juce::Colour(0xff6d6d6d));
            for(int i=1;i<=15;i++)
                if(i%4!=0)g.drawLine(getWidth()/16*i,(float)getHeight()/1.2f,getWidth()/16*i,(float)getHeight(),2.5f);
            
            vector<float> quarter = getQuarterVector();
            g.setColour (juce::Colour(0xff00C40f));
            for(float f : quarter){
                g.drawLine(getWidth()/100*f,(float)getHeight()/1.2f,getWidth()/100*f,(float)getHeight(),2.0f);
            }
        
        }
    }
    
    Colour getColorforHit(float drum_hit,vector<float> peaks){
        int min = 10000;
        for (float f : peaks) {
            if(abs(drum_hit-f)<min)
                min = abs(drum_hit-f);
        }
        int sampleRate = numberofSamples / (float) transportSource.getLengthInSeconds();
        int one_ms_insample = sampleRate / 1000;
        
        double error_in_ms = min / one_ms_insample;
        
        if(error_in_ms>60)return Colour(0xffE91010);
        if(error_in_ms>40)return Colour(0xffF19409);
        if(error_in_ms>20)return Colour(0xffF1F009);
        return Colour(0xff63F109);
        
    }
    
    
    
    vector<float> getQuarterVector(){
        vector<float> percentages;
        if(!transportSource.isPlaying())return percentages;
        float currentLastHit_asQuarter;
        float duration = (float) transportSource.getLengthInSeconds();
        float audioPosition = (float) transportSource.getCurrentPosition();
        float currentSample = (audioPosition / duration) * (float) numberofSamples;
        int first_hit_index = 0;
        if(currentFirstHit_asQuarter == 0){
            if(calculated_peaks.size() == 0)return percentages;
            currentFirstHit_asQuarter = calculated_peaks[0];
            currentLastHit_asQuarter = calculated_peaks[4];
        }
        else{
            auto it = find(calculated_peaks.begin(),calculated_peaks.end(),currentFirstHit_asQuarter);
            first_hit_index = distance(calculated_peaks.begin(),it);
            if(calculated_peaks[first_hit_index+4]<currentSample){
                currentFirstHit_asQuarter = calculated_peaks[first_hit_index+4];
                first_hit_index += 4;
                currentLastHit_asQuarter = calculated_peaks[first_hit_index+4];
            }else{
                currentLastHit_asQuarter = calculated_peaks[first_hit_index+4];
            }
        }
        float onepercent = (currentLastHit_asQuarter - currentFirstHit_asQuarter) / 100;
        
        for(float f : drum_hits){
            if(f>currentFirstHit_asQuarter && f<currentLastHit_asQuarter){
                percentages.push_back((f-currentFirstHit_asQuarter) / onepercent);
            }
        }
        
        for(float f: percentages){
            if(f>100)f = 100;
            if(f<0)f = 0;
        }
        return percentages;
    }
    
    bool keyPressed(const KeyPress &k, Component *c) override {
        //drumhit();
        return true;
    }
  
    void updateVectors(vector<float> peaks_,vector<float> drum_hits_){
        calculated_peaks = peaks_;
        drum_hits = drum_hits_;
        for(int z = 0;z<5;z++){
        calculated_peaks.push_back(calculated_peaks[calculated_peaks.size()-1]+
                                    (calculated_peaks[calculated_peaks.size()-1]-
                                     calculated_peaks[calculated_peaks.size()-2]));
        }
        
    }
    
    void updatedrumHits(vector<float> drum_hits_){
        drum_hits = drum_hits_;
        
    }
    void updatePeaks(vector<float> peaks_){
        calculated_peaks = peaks_;
        for(int z = 0;z<5;z++){
        calculated_peaks.push_back(calculated_peaks[calculated_peaks.size()-1]+
                                    (calculated_peaks[calculated_peaks.size()-1]-
                                     calculated_peaks[calculated_peaks.size()-2]));
        }
    }
    
    void update_base_Peaks(vector<float> peaks_){
        peaks = peaks_;
    }
    void setNumberofSamples(int n){
        numberofSamples = n;
        currentFirstHit_asQuarter = 0;
    }
    
    void mouseDown (const juce::MouseEvent& event) override
    {
        auto duration = transportSource.getLengthInSeconds();

        if (duration > 0.0)
        {
            auto clickPosition = event.position.x;
            auto audioPosition = (clickPosition / (float) getWidth()) * duration;

            transportSource.setPosition (audioPosition);
        }
    }
    
    
    void setInterval(int interval_){
        interval_time = interval_;
    }
    

private:

    void timerCallback() override
    {
        repaint();
    }
    
    juce::AudioTransportSource& transportSource;
    vector<float> calculated_peaks;
    vector<float> peaks; //Just for screenshots and testing
    vector<float> drum_hits;
    double startingState;
    float currentFirstHit_asQuarter = 0;
    int numberofSamples;
    int interval_time = 5;
    int counter = 0;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (SimplePeakOverlay)
};

class Logic{
    AudioTransportSource& transportSource;
    int sampleRate;
    float adjustedBPM = 0;
    float startingState = 0;
    vector<float> data;
    float micLatency_mseconds = 0;
    float latencyfromenvelope_insample = 1000;
    chrono::time_point<chrono::high_resolution_clock> startTimeV;
    vector<float> previous_inputs;//Stores last 4 inputs from previous buffer
    vector<float> previous_outputs;//Same for outputs
    float mic_sensitivity = 0.05;
    
public:
    vector<float> drum_hits;
    vector<float> calculated_peaks;
    vector<float> base_peaks; // For testing
    bool calibrating_out = false;
    bool calibrating_in = false;
    
    Logic (AudioTransportSource& transportSourceToUse)
    : transportSource (transportSourceToUse){;}
    double getAdjustedBPM(){
        return adjustedBPM;
    }
    void setSamplerate(int rate){
        sampleRate = rate;
    }
    
    void setAdjustedBPM(float multiplier){
        adjustedBPM*=multiplier;
    }

    void initialize(vector<float> data_){
        data = data_;
        vector<float> peaks = calculate_peaks(data);
        base_peaks = peaks;
        vector<double> envelope = calculate_envelope(data,false);
        calculated_peaks = updatePeaks(peaks);
        cout<<"Calculated Peak size : "<<calculated_peaks.size()<<endl;
        calculateBPM_andstartSmaple(peaks,envelope);
        cout<<"Calculated Peak size : "<<calculated_peaks.size()<<endl;
        drum_hits.clear();
    }
    
    vector<double> calculate_envelope(vector<float> data,bool use_previous_data){
        vector<double> data_abs;
        vector<double> envelope;
        
        double a0,a1,a2,a3,a4,b0,b1,b2,b3,b4;
          
          a0 = 1;
          a1 = -3.99255385408073;
          a2 =  5.9776892722649;
          a3 = -3.97771692191996;
          a4 =  0.992581503801471;
          
          b0 =  0.041053549 / 10000000000;
          b1 =  0.164214197 / 10000000000;
          b2 =  0.246321296 / 10000000000;
          b3 =  0.164214197 / 10000000000;
          b4 =  0.041053549 / 10000000000;
        
        
        for (float f : data) {
            data_abs.push_back(abs(f));
            envelope.push_back(0);
        }
            
        if(use_previous_data){
            if(!previous_inputs.empty() && !previous_outputs.empty()){
                envelope[0] = (-1) * a1 * previous_outputs[3] - a2 * previous_outputs[2] - a3 * previous_outputs[1] - a4 * previous_outputs[0] +
                  b0 * data_abs[0] + b1 * previous_inputs[3] +  b2 * previous_inputs[2] +  b3 * previous_inputs[1] +  b4 * previous_inputs[0];
                
                envelope[1] = (-1) * a1 * envelope[0] - a2 * previous_outputs[3] - a3 * previous_outputs[2] - a4 * previous_outputs[1] +
                b0 * data_abs[1] +   b1 * data_abs[0] + b2 * previous_inputs[3] +  b3 * previous_inputs[2] +  b4 * previous_inputs[1];
                
                envelope[2] = (-1) * a1 * envelope[1] - a2 * envelope[0] - a3 * previous_outputs[3] - a4 * previous_outputs[2] +
                b0 * data_abs[2] +   b1 * data_abs[1] + b2 * data_abs[0] + b3 * previous_inputs[3] +  b4 * previous_inputs[2];
                
                envelope[3] = (-1) * a1 * envelope[2] - a2 * envelope[1] - a3 * envelope[0] - a4 * previous_outputs[3] +
                b0 * data_abs[3] +   b1 * data_abs[2] + b2 * data_abs[1] + b3 * data_abs[0] + b4 * previous_inputs[3];
            }
        }
        for(int i=4;i<data_abs.size();i++){
            envelope[i] = (-1) * a1 * envelope[i-1] - a2 * envelope[i-2] - a3 * envelope[i-3] - a4 * envelope[i-4] +
            b0 * data_abs[i] + b1 * data_abs[i-1] + b2 * data_abs[i-2] + b3 * data_abs[i-3] + b4 * data_abs[i-4];
        }
    
        if(use_previous_data)
            for(int i=(int)data_abs.size()-4;i<data_abs.size();i++){
                previous_outputs.push_back(envelope[i]);
                previous_inputs.push_back(data_abs[i]);
            }
        
        
        
        return envelope;
    }
    
    vector<float> calculate_peaks(vector<float> data){
        
        vector<float> differ;
        vector<double> envelope;
        
        envelope = calculate_envelope(data,false);
        
        
        for(int i=0;i<envelope.size()-1;i++){
            differ.push_back(envelope[i+1]-envelope[i]);
        }
        
        vector<double> abs_differ;
        
        for (double d : differ) abs_differ.push_back(abs(d));
        
        double mean = accumulate(abs_differ.begin(),abs_differ.end(),0.0) / differ.size();
        
        vector<float> peaks;
        
        for(int i=1;i<differ.size()-1;i++){
            if(differ[i-1]<differ[i] && differ[i]>differ[i+1] && differ[i]>mean*4)
                peaks.push_back(i);
        }
        return peaks;
    }
    
    void calculateBPM_andstartSmaple(vector<float> peaks,vector<double> envelope){
        vector<int> peak_distances;
        
        int possibleBPM = 0;
        int maxcounter = 0;
        
        if(!peaks.empty()){
            for (int i=0;i<peaks.size();i++) {
                int counter = 0;
                for(int j=i+1;j<peaks.size();j++){
                    if(counter++>5)break;
                    peak_distances.push_back(abs(peaks[j]-peaks[i]));
                }
            }
         
            for (int i=0; i<peak_distances.size(); i++) {
                float sum = peak_distances[i];
                float counter = 1;
                for (int j=i+1; j<peak_distances.size(); j++) {
                    if(abs(peak_distances[i]-peak_distances[j])<1000){
                        sum+=peak_distances[j];
                        counter++;
                        peak_distances.erase(peak_distances.begin()+j);
                        j = i+1;
                    }
                }
                if(counter>10){
                    if(counter>maxcounter){
                        possibleBPM = round(sampleRate/(sum/counter)*60);
                        while(possibleBPM<70)possibleBPM*=2;
                        maxcounter = counter;
                    }
                }
                peak_distances.erase(peak_distances.begin()+i);
                i = 0;
            }
        }
        peaks.push_back(5000);
        if(possibleBPM == 0 || possibleBPM > 200){
            possibleBPM = 60;
            peaks.clear();
            peaks.push_back(5000);
        }
        double possiblePeakDistance;
        double max = 0;
        startingState = 0;
        
        for (float k=-20; k<20; k+=0.25) {
            possiblePeakDistance = sampleRate/((possibleBPM+k)/60);
            for (int i = 0; i < 350; i++) {
                double sum = 0;
                
                for (int j = 0; j < 30; j++) {
                    if((peaks[0]+i*100)+j*round(possiblePeakDistance)< envelope.size())
                       sum += envelope[(peaks[0]+i*100)+j*round(possiblePeakDistance)];
                }
                if(sum > max){
                    max = sum;
                    startingState = peaks[0]+i*100;
                    adjustedBPM = possibleBPM + k;
                }
            }
        }
        
        if(adjustedBPM<80)adjustedBPM*=2;
        
        startingState -= latencyfromenvelope_insample;
        calculated_peaks = updatePeaks(peaks);
    }
    
    void insertBuffer(vector<float> buffer){
        float audioPosition = (float) transportSource.getCurrentPosition();
        
        int samplePosition = audioPosition*sampleRate;

        vector<double> envelope_of_buffer = calculate_envelope(buffer, true);
        
        auto it = max_element(envelope_of_buffer.begin(),envelope_of_buffer.end());
        auto i = distance(envelope_of_buffer.begin(),it);
        
        float sampleHit = samplePosition - (envelope_of_buffer.size()-i);
        
        for(int i=0;i<*max_element(envelope_of_buffer.begin(),envelope_of_buffer.end())*1000;i++)
            cout<<"-";
        cout<<endl;
        
        if(*it<mic_sensitivity)return;
        
        if(!drum_hits.empty()){
            if(abs(drum_hits[drum_hits.size()-1]-sampleHit)>7000) {
                drum_hits.push_back(sampleHit - micLatency_mseconds / 1000 * sampleRate);
                drumhit();
            
            }
        }else{
            drum_hits.push_back(sampleHit - micLatency_mseconds / 1000 * sampleRate);
            drumhit();
        }
    }
    
    
    void updatePeaks(){
        vector<float> new_peaks;
        
        double possiblePeakDistance = sampleRate/((adjustedBPM)/60);
        
        int j=0;
        while(startingState+j*round(possiblePeakDistance)<data.size()-1) {
            new_peaks.push_back(startingState+j*round(possiblePeakDistance));
            j++;
        }
        calculated_peaks = new_peaks;
    }
    
    vector<float> updatePeaks(vector<float> peaks){
        cout<<"Peak size : "<<peaks.size()<<endl;
        vector<float> adjusted_peaks;
        
        double possiblePeakDistance = sampleRate/((adjustedBPM)/60);
        cout<<"Peak distance : "<<possiblePeakDistance<<endl;
        
        int j=0;
        while(startingState+j*round(possiblePeakDistance)<data.size()-1) {
            adjusted_peaks.push_back(startingState+j*round(possiblePeakDistance));
            j++;
        }
        cout<<"Adjusted Peak size : "<<adjusted_peaks.size()<<endl;
        return adjusted_peaks;
    }
    
    void calibrateMic(){

        calibrating_out = true;
    }
    
    void startTime(){
        startTimeV = chrono::high_resolution_clock::now();
    }
    
    void stopTime(vector<float> buffer){
        auto stopTimeV = chrono::high_resolution_clock::now();
        chrono::duration<double> diff = stopTimeV - startTimeV;
        auto first_sound = upper_bound(buffer.begin(),buffer.end(),0.3);
        double position =(double)(first_sound-buffer.begin());
        double delay_from_buffer = (buffer.size()-position)/(double)44100;
        micLatency_mseconds = (diff.count()-delay_from_buffer)*1000;
        cout<<"Input latency is: "<<micLatency_mseconds<<"ms"<<endl;
    }
    void setSensitivity(float value){
        mic_sensitivity = value;
    }
};




//------------------------------------------------------------------------------

class MainContentComponent   : public juce::AudioAppComponent,
                               public juce::ChangeListener,
                            public Slider::Listener
{
public:
    MainContentComponent()
      : state (Stopped),
        other_window("Device Manager",Colour(0xff494949),0,true),
        audioSetupComp(deviceManager,0,256,0,256,false,false,false,true),
        thumbnailCache (5),
        thumbnailComp (512, formatManager, thumbnailCache,transportSource),
        positionOverlay (transportSource),
        peakOverlay(transportSource),
        logic(transportSource)
    {
        addAndMakeVisible (&openButton);
        openButton.setButtonText ("Open...");
        openButton.onClick = [this] { openButtonClicked(); };

        addAndMakeVisible (&playButton);
        playButton.setButtonText ("Play");
        playButton.onClick = [this] {playButtonClicked(); };
        playButton.setColour (juce::TextButton::buttonColourId, juce::Colours::green);
        playButton.setEnabled (false);

        addAndMakeVisible (&stopButton);
        stopButton.setButtonText ("Stop");
        stopButton.onClick = [this] { stopButtonClicked(); };
        stopButton.setColour (juce::TextButton::buttonColourId, juce::Colours::red);
        stopButton.setEnabled (false);

        addAndMakeVisible (&thumbnailComp);
        addAndMakeVisible (&positionOverlay);
        addAndMakeVisible(&peakOverlay);

        addAndMakeVisible(&BPMlabel);
        addAndMakeVisible (&plusBPM);
        plusBPM.setButtonText("x2 BPM");
        plusBPM.onClick = [this] {logic.setAdjustedBPM(2.0f);logic.updatePeaks();peakOverlay.updatePeaks(logic.calculated_peaks);BPMlabel.setText(to_string(logic.getAdjustedBPM())+" BPM",dontSendNotification);};
        addAndMakeVisible (&minusBPM);
        minusBPM.setButtonText("/2 BPM");
        minusBPM.onClick = [this] {logic.setAdjustedBPM(0.5f);logic.updatePeaks();peakOverlay.updatePeaks(logic.calculated_peaks);BPMlabel.setText(to_string(logic.getAdjustedBPM())+" BPM",dontSendNotification);};
        
        addAndMakeVisible(mic_sensitivity_slider);
        mic_sensitivity_slider.setRange(0.0,10.0);
        mic_sensitivity_slider.setSkewFactorFromMidPoint(0.1);
        mic_sensitivity_slider.setValue(0.3,dontSendNotification);
        mic_sensitivity_slider.addListener(this);
        mic_sensitivity_slider.setNumDecimalPlacesToDisplay(3);
        
        setSize (600, 400);

        formatManager.registerBasicFormats();
        transportSource.addChangeListener (this);

        setAudioChannels (2, 2);
        
        drumhit = [this](){
            peakOverlay.updatedrumHits(logic.drum_hits);
        };
        
        addAndMakeVisible(calibrateLatencybutton);
        calibrateLatencybutton.setButtonText("Calibrate Mic");
        calibrateLatencybutton.onClick = [this]{
            logic.calibrateMic();
        };
        
        other_window.addToDesktop();
        other_window.setContentOwned(&audioSetupComp,true);
        deviceManager.addChangeListener(this);
        audioSetupComp.setItemHeight(25);
        
        addAndMakeVisible (audiosetupbutton);
        audiosetupbutton.setButtonText ("Device Setup");
        audiosetupbutton.onClick = [this] { audiosetupbuttonClicked(); };

        addAndMakeVisible(thumbnail_size_insec);
        thumbnail_size_insec.setRange(4,30,1);
        thumbnail_size_insec.setValue(5,dontSendNotification);
        thumbnail_size_insec.setTextValueSuffix(" sec");
        thumbnail_size_insec.addListener(this);
    }

    ~MainContentComponent() override
    {
        shutdownAudio();
    }

    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override
    {
        transportSource.prepareToPlay (samplesPerBlockExpected, sampleRate);
    }

    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override
    {
        if(logic.calibrating_out){
            for (auto channel = 0; channel < bufferToFill.buffer->getNumChannels(); ++channel)
                    {
                        // Get a pointer to the start sample in the buffer for this audio output channel
                        auto* buffer = bufferToFill.buffer->getWritePointer (channel, bufferToFill.startSample);
             
                        // Fill the required number of samples with noise between -0.125 and +0.125
                        for (auto sample = 0; sample < bufferToFill.numSamples; ++sample)
                            buffer[sample] = 1.0;
                    }
            logic.calibrating_out = false;
            logic.calibrating_in = true;
            logic.startTime();
            return;
        }
        
        auto* device = deviceManager.getCurrentAudioDevice();
        auto activeInputChannels  = device->getActiveInputChannels();
        auto activeOutputChannels = device->getActiveOutputChannels();
        auto maxInputChannels  = activeInputChannels .getHighestBit() + 1;
        auto maxOutputChannels = activeOutputChannels.getHighestBit() + 1;
        for (auto channel = 0; channel < maxOutputChannels; ++channel){
                    {
                        auto actualInputChannel = channel % maxInputChannels;
         
                        if (! activeInputChannels[channel])
                        {
                            bufferToFill.buffer->clear (channel, bufferToFill.startSample, bufferToFill.numSamples);
                        }
                        else
                        {
                            vector<float> current_buffer;
                            float const* pointer = bufferToFill.buffer->getReadPointer (actualInputChannel);
                            double max = 0;
                            for(int i=0;i<bufferToFill.buffer->getNumSamples();i++){
                                current_buffer.push_back(pointer[i]);
                                if(max<pointer[i]) max = pointer[i];
                            }
                            if(max>mic_sensitivity_slider.getValue() && logic.calibrating_in){
                                logic.calibrating_in = false;
                                logic.stopTime(current_buffer);
                                
                            }
                            if(file_open){
                                logic.insertBuffer(current_buffer);
                            }
                        }
                    }
        }
        
        if (readerSource.get() == nullptr)
            bufferToFill.clearActiveBufferRegion();
        else
            transportSource.getNextAudioBlock (bufferToFill);
    }

    void releaseResources() override
    {
        transportSource.releaseResources();
    }

    void resized() override
    {
        openButton.setBounds (10, 10, 100, 30);
        playButton.setBounds (10, 50, 100, 30);
        stopButton.setBounds (10, 90, 100, 30);
        
        audiosetupbutton.setBounds(140,10,100,30);
        calibrateLatencybutton.setBounds(140,50,100,30);
        
        BPMlabel.setBounds(380,10,100,30);
        plusBPM.setBounds(270,10,80,30);
        minusBPM.setBounds(270,50,80,30);
        
        mic_sensitivity_slider.setBounds(140,90,200,30);
        thumbnail_size_insec.setBounds(380,90,200,30);
        
        juce::Rectangle<int> thumbnailBounds (10, 150, getWidth() - 20, (getHeight() - 170)/1.3);
        thumbnailComp.setBounds (thumbnailBounds);
        positionOverlay.setBounds (thumbnailBounds);
        juce::Rectangle<int> thumbnailBounds2 (10, 150, getWidth() - 20, (getHeight() - 170));
        peakOverlay.setBounds(thumbnailBounds2);
        
        audioSetupComp.setBounds(10,10,500,200);
    }
    
    

    void changeListenerCallback (juce::ChangeBroadcaster* source) override
    {
        if (source == &transportSource)
            transportSourceChanged();
    }
    
    void sliderValueChanged (juce::Slider* slider) override
        {
            thumbnailComp.setInterval(thumbnail_size_insec.getValue());
            positionOverlay.setInterval(thumbnail_size_insec.getValue());
            peakOverlay.setInterval(thumbnail_size_insec.getValue());
            logic.setSensitivity(mic_sensitivity_slider.getValue());
        }

private:
    enum TransportState
    {
        Stopped,
        Starting,
        Playing,
        Stopping
    };

    void changeState (TransportState newState)
    {
        if (state != newState)
        {
            state = newState;

            switch (state)
            {
                case Stopped:
                    stopButton.setEnabled (false);
                    playButton.setEnabled (true);
                    transportSource.setPosition (0.0);
                    break;

                case Starting:
                    playButton.setEnabled (false);
                    transportSource.start();
                    file_open = true;
                    break;

                case Playing:
                    stopButton.setEnabled (true);
                    break;

                case Stopping:
                    transportSource.stop();
                    file_open = false;
                    break;

                default:
                    jassertfalse;
                    break;
            }
        }
    }

    void transportSourceChanged()
    {
        if (transportSource.isPlaying())
            changeState (Playing);
        else
            changeState (Stopped);
    }

    void openButtonClicked()
    {
        
        
        
        juce::FileChooser chooser ("Select a Wave file to play...",
                                   {},
                                   "*.wav;*.mp3");

        if (chooser.browseForFileToOpen())
        {
            juce::File file = chooser.getResult();

            if (AudioFormatReader * reader = formatManager.createReaderFor (file))
            {
                AudioFormatReaderSource* AudioReader = new juce::AudioFormatReaderSource (reader, true);
                std::unique_ptr<juce::AudioFormatReaderSource> newSource (AudioReader);
                cout<<"Sample rate : "<<reader->sampleRate<<endl;
                AudioSampleBuffer *audiobuffer = new AudioSampleBuffer();
                audiobuffer->setSize((int) reader->numChannels, (int) reader->lengthInSamples);
                reader->read(audiobuffer,0,(int) reader->lengthInSamples,0,true,true);
                float const * data = audiobuffer->getReadPointer(0);
                std::vector<float> data_vec;
                for(int i=0;i<audiobuffer->getNumSamples();i++){
                    data_vec.push_back(data[i]);
                }
                
                delete audiobuffer;
                logic.setSamplerate(reader->sampleRate);
                logic.initialize(data_vec);
                BPMlabel.setText(to_string(logic.getAdjustedBPM())+" BPM",dontSendNotification);
                
                peakOverlay.setNumberofSamples(data_vec.size());
                peakOverlay.updateVectors(logic.calculated_peaks,logic.drum_hits);
                peakOverlay.update_base_Peaks(logic.base_peaks);
                
                transportSource.setSource (newSource.get(), 0, nullptr, reader->sampleRate);
                playButton.setEnabled (true);
                thumbnailComp.setFile (file);
                readerSource.reset (newSource.release());
            }
        }
    }

    void playButtonClicked()
    {
        changeState (Starting);
    }

    void stopButtonClicked()
    {
        changeState (Stopping);
    }
    
    void audiosetupbuttonClicked(){
        other_window.setVisible(!other_window.isVisible());
    }

    //==========================================================================
    TextButton openButton;
    TextButton playButton;
    TextButton stopButton;
    TextButton plusBPM;
    TextButton minusBPM;
    TextButton calibrateLatencybutton;
    TextButton audiosetupbutton;
    Label BPMlabel;
    Slider thumbnail_size_insec;
    
    DocumentWindow other_window;
    AudioDeviceSelectorComponent audioSetupComp;
    AudioFormatManager formatManager;
    std::unique_ptr<juce::AudioFormatReaderSource> readerSource;
    AudioTransportSource transportSource;
    TransportState state;
    AudioThumbnailCache thumbnailCache;
    SimpleThumbnailComponent thumbnailComp;
    SimplePositionOverlay positionOverlay;
    SimplePeakOverlay peakOverlay;
    Logic logic;
    Slider mic_sensitivity_slider;
    bool file_open = false;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainContentComponent)
};
