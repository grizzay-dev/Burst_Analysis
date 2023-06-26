#ifndef burst_h
#define burst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TH1.h"
#include "TAxis.h"
#include "TGraph.h"
#include <cmath>
#include <vector>
#include <string>
#include <numeric>

using namespace std;

class Burst {
    //Method of burst detection adapated from: https://www.sciencedirect.com/science/article/pii/S1002007108003432
    public : 
        //Spike train
        TH1C* spike_train;          //{nullptr};
        TH1* isi_histogram {nullptr};         //isi (CMA Method)
        TH1* isi_cs_histogram;      //cumulative sum (CMA Method)
        TH1* isi_cma_histogram;     //cumulative average (CMA Method)
        
        TH1C* spike_train_V2;
        TH1F* V;

        //spike train properties
        int                     total_recording_time;
        int                     injection_start;
        int                     injection_end;

        double                  seconds_per_bin;


        //Burst processing properties
        int                     id;
        int                     n_bins;
        int                     n_spikes{ 0 };
        int                     isi_length;
        int                     ln_length;

        //My Method
        double                  k{ 1 };
        double                  isi_sd{ 0 };
        double                  isi_mean{ 0 };
        double                  isi_cutoff{ 2 }; //set to 2 for standard processing
        double                  isi_threshold{ 0 };
        double                  ln_mean{ 0 };

        vector<double>          spikes_x;
        vector<double>          isi_sequence;
        vector<double>          ln_sequence;
        
        //CMA Method
        double                  isi_skewness{ 0 };
        double                  isi_skewness_manual{ 0 };
        double                  isi_argmax{ 0 };

        // To calculate isi_threshold_primary and isi_threshold_secondary
        double                  alpha_primary{ 0 };
        double                  alpha_neighbour{ 0 };
        double                  cma_primary{ 0 };
        double                  cma_neighbour{ 0 };

        double                  isi_threshold_baseline{ 0 }; //cma_argmax_bin * seconds_per_bin
        double                  isi_threshold_primary{ 0 };
        double                  isi_threshold_neighbour{ 0 };

        int                     cma_argmax_bin{ 0 };
        int                     cma_n_bursts{ 0 };


        vector<int>             isi_cumulative_sum; //(CMA Method)
        vector<double>          isi_cumulative_avg; //(CMA Method)
        vector<double>          cma_neighbour_locations;
        vector<vector<int>>     cma_burst_locations;

        //Burst metrics
        bool                    bursting{ false };
        int                     burst_type{ 0 };
        int                     n_bursts{ 0 };

        vector<vector<int>>     burst_locations;

        double                  burst_total_duration{ 0 };
        double                  burst_average_duration{ 0 };
        double                  burst_frequency{ 0 };
        double                  ML{ 0 };

        //Constructor
        Burst(int entry_id, TH1C* entry, TH1F* voltage, double constant = 1){
            spike_train = entry;
            V = voltage;
            id = entry_id;
            k = constant;

            SpkDetection(V);
            
            //set seconds per bin using bin width
            seconds_per_bin = spike_train_V2->GetBinWidth(0);
            n_bins = spike_train_V2->GetNbinsX();

            //get x-axis
            TAxis* t = spike_train_V2->GetXaxis();
            injection_start = t->GetXmin();
            injection_end = t->GetXmax();
        };

        virtual void    PrintMetrics();
        virtual void    Spikes();
        virtual void    ISI();
        virtual void    ISIHistogram();  //CMA and skewness
        virtual void    CalculateSkewness();
        virtual void    CalculateAlphaValues(); //CMA
        virtual double  FindClosestValue(double value, vector<double> vect, int start_bin);
        virtual void    CMA(); //CMA
        virtual void    LN();
        virtual void    DetectBurst();
        virtual void    CMADetectBurst();
        virtual void    ProcessNeuron();
        virtual void    ProcessTree();
        virtual void    TreeToCSV();
        virtual void    SpkDetection(TH1F* V);
        virtual void    CleanSpikes(double window);
};

void PrintSequence(vector<double> sequence, string name = "Sequence:"){
    int length = sequence.size();
    for(int i = 0; i < length; i++){
        if(i == 0){
            printf("\n%-40s(%i) ", name.c_str(), length);
        }
        cout << sequence.at(i);
        if(i != length -1) cout << ", ";
    }
    cout << " - END";
}

double CumulativeAverage(int start, int end, vector<double> vect){
    double sum = 0.0;
    double delta, n, avg;

    for (int x = start; x < end; x++) {
        delta = vect[x + 1] - vect[x];
        //delta = vect.at(x + 1) - vect.at(x);
        sum += delta;
    }
    n = (double)end - (double)start;
    avg = sum / n;
    return avg;
}

double Burst::FindClosestValue(double value, vector<double> vect, int start_bin){
    double delta, prev_delta, near_delta, near_val;
    int i_near;
    near_delta = value;
    prev_delta = value + 1;
    
    for(int i = start_bin; i < vect.size(); i++){
        
        delta = abs(vect.at(i) - value);
        if(delta > prev_delta){
            //break;
        }
        if(delta < near_delta){
            i_near = i;
            near_delta = delta;
            near_val = i_near * seconds_per_bin;
        }

        prev_delta = delta;
    }
    return near_val;
}

void Burst::CalculateAlphaValues() {
    if (isi_skewness_manual < 1) {
        alpha_primary       = 1;
        alpha_neighbour     = 0.5;
        
    }
    else if(isi_skewness_manual < 4) {
        alpha_primary       = 0.7;
        alpha_neighbour     = 0.5;

    }
    else if(isi_skewness_manual < 9) {
        alpha_primary       = 0.5;
        alpha_neighbour     = 0.3;
    }
    else {
        alpha_primary       = 0.3;
        alpha_neighbour     = 0.1;
    }
    cma_primary = isi_argmax * alpha_primary;
    cma_neighbour = isi_argmax * alpha_neighbour;

    isi_threshold_primary = FindClosestValue(cma_primary, isi_cumulative_avg, cma_argmax_bin);
    isi_threshold_neighbour = FindClosestValue(cma_neighbour, isi_cumulative_avg, cma_argmax_bin);
}

void Burst::ISIHistogram() {
    // Only continue if atleast 2 spikes.
    if (n_spikes < 2) return;

    //Create histogram:
    double x_max = *max_element(isi_sequence.begin(), isi_sequence.end());
    //x_max += (x_max * 0.1);  //add 10% space to end of x-axis
    double isi_bins = x_max / seconds_per_bin;
    if(isi_histogram != nullptr) {
        isi_histogram->Reset();
        }
        else {
            isi_histogram = new TH1I("isi_histogram", "ISI Histogram", isi_bins, 0, x_max);
        } 

    //Fill histogram
    for (double isi : isi_sequence) {
        isi_histogram->Fill(isi);
    }

    //Save skewness
    isi_skewness = isi_histogram->GetSkewness();

}

void Burst::CMA() {

    int isi_bins;
    isi_bins = isi_histogram->GetNbinsX();

    int current_sum = 0;
    for (int i = 0; i < isi_bins; i++) {
        current_sum += isi_histogram->GetBinContent(i);
        isi_cumulative_sum.push_back(current_sum);
        isi_cumulative_avg.push_back(isi_cumulative_sum.at(i) / ((double)i+1));
        
        if (isi_argmax < isi_cumulative_avg.at(i)) {
            isi_argmax = isi_cumulative_avg.at(i);
            cma_argmax_bin = i;
        }
    }

    isi_threshold_baseline = cma_argmax_bin * seconds_per_bin;
}

void Burst::PrintMetrics(){

    printf("\n\nBURST METRICS");
    printf("\n%-40s%-20i", "ID", id);
    printf("\n%-40s%-20i", "Burst Type", burst_type);
    printf("\n%-40s%-20i", "# of Bursts", n_bursts);
    printf("\n%-40s%-20f", "Total Burst Duration", burst_total_duration);
    printf("\n%-40s%-20f", "Avg burst duration", burst_average_duration);
    printf("\n%-40s%-20f", "Burst frequency", burst_frequency);
    printf("\n%-40s%-20f", "ML", ML);
    printf("\n%-40s%-20f", "ISI Threshold", isi_threshold);
    printf("\n%-40s%-20f", "ISI ArgMax", isi_argmax);
    printf("\n%-40s%-20i", "ISI ArgMax Bin", cma_argmax_bin);
    printf("\n%-40s%-20f", "(CMA) Primary Alpha", cma_primary);
    printf("\n%-40s%-20f", "(CMA) Neighbour Alpha", cma_neighbour);
    printf("\n%-40s%-20f", "(CMA) Baseline Threshold", isi_threshold_baseline);
    printf("\n%-40s%-20f", "(CMA) Primary Threshold", isi_threshold_primary);
    printf("\n%-40s%-20f", "(CMA) Neighbour Threshold", isi_threshold_neighbour);
    printf("\n%-40s%-20f", "ISI Skewness", isi_skewness);
    printf("\n");
    PrintSequence(spikes_x, "Spike Train");
    PrintSequence(isi_sequence, "ISI Seq");
    PrintSequence(cma_neighbour_locations, "Neighbour Spikes");
    
    if (n_bursts > 0) {
        int x = 1;
        printf("\nBurst Locations:");
        for (vector<int> i : burst_locations) {
            printf("\n%41s%-5i%.2f%6s   %.2f", "#", x, spikes_x.at(i.at(0)), "--->", spikes_x.at(i.at(1)));
            x++;
        }
            
    }

    if (cma_n_bursts > 0){
        int x = 1;
        printf("\n(CMA)Burst Locations:");
        for (vector<int> i : cma_burst_locations) {
            printf("\n%41s%-5i%.2f%6s   %.2f", "#", x, spikes_x.at(i.at(0)), "--->", spikes_x.at(i.at(1)));
            x++;
        }
    }
}

void Burst::Spikes() {

    //calculate n_spikes
    for (int bin = 0; bin < n_bins; bin++) {
        double bin_val = spike_train_V2->GetBinContent(bin);
        if (bin_val > 0) {
            n_spikes = n_spikes + 1;
            spikes_x.push_back((bin * seconds_per_bin) + injection_start);
        }
    }
}

void Burst::ISI(){
    if(n_spikes==0) return;
    isi_length = n_spikes - 1;

    for(int x = 0; x < isi_length; x++){
        isi_sequence.push_back(spikes_x.at(x + 1) - spikes_x.at(x));
    }
    
    double sum = accumulate(isi_sequence.begin(), isi_sequence.end(), 0.0);
    isi_mean = sum / double(isi_sequence.size());

    for (double x : isi_sequence) {
        isi_sd += pow(x - isi_mean, 2);
    }

    isi_sd = sqrt(isi_sd/double(isi_sequence.size()));
    isi_threshold = isi_sd / isi_mean;
    
}

void Burst::LN(){
    for(int x = 0; x < isi_length; x++){
        if(isi_sequence.at(x) < isi_mean) {
            ln_sequence.push_back(isi_sequence.at(x));
        }
    }
    double ln_sum = accumulate(ln_sequence.begin(), ln_sequence.end(), 0.0);
    ln_mean = (double)ln_sum / (double)(ln_sequence.size());
}

void Burst::CMADetectBurst(){
    int start = 0;
    int end = 0;
    int min_isi_required = 2; //equivalent to 3 spikes
    int n_isi = 0;
    bool detected = false;

    while(end < isi_length){
        if(isi_sequence.at(end) < isi_threshold_neighbour){
            n_isi += 1;
            if(isi_sequence.at(end) > isi_threshold_primary){
                cma_neighbour_locations.push_back(spikes_x.at(end));
            }
            end += 1;

            if (end == isi_length-1){
                if(n_isi >= min_isi_required){
                    // ORIGINAL vector<int> x1_x2 = {start, (end+1)};
                    // MOSTLY WORKING vector<int> x1_x2 = {start, (end)};
                    vector<int> x1_x2 = {start, (end)};
                    cma_burst_locations.push_back(x1_x2);
                    cma_n_bursts += 1;
                }
            }
        } else {
            if(n_isi >= min_isi_required){
                vector<int> x1_x2 = {start, (end)};
                cma_burst_locations.push_back(x1_x2);
                cma_n_bursts += 1;
            }
            n_isi = 0;
            start = end +1;
            end = start + 1;
        }
    }
}

void Burst::DetectBurst(){
    int start = 0;
    int end =  2;
    bool detected = false;

    while (end < isi_length + 2){
        double burst_average = 1.0;
        burst_average = CumulativeAverage(start, end, spikes_x);

        if((burst_average < (ln_mean * k)) && (burst_average > 0)){
            detected = true;
            end += 1;
            bursting = true;
            burst_type = 1;
        } else {
            if(detected){
                n_bursts += 1;
                burst_total_duration += spikes_x.at(end - 1) - spikes_x.at(start);
                vector<int> x1_x2 = { start, (end - 1) };
                burst_locations.push_back(x1_x2);

                //reset for next sequence
                start = end;
                end = start + 2;
                detected = false;
                
            } else {
                start = end - 1;
                end = start + 2;
                detected = false;
            }
        }
    }
    
    //Set remaining burst metrics
    ML = ln_mean;
    if(n_bursts > 0){
        burst_average_duration = burst_total_duration / n_bursts;
        burst_frequency = double(n_bursts) / (injection_end - injection_start); //current injection from 400ms to 100ms 
    }
}

void Burst::ProcessNeuron(){
    //method: 0 = standard, 1 = cma
    cout << "\rProcessing entry: " << id;
    
    Spikes();
    CleanSpikes(5);
    spikes_x.clear();
    n_spikes = 0;
    Spikes();

    

    if(n_spikes > 0){
        ISI();
        ISIHistogram();
        CMA();
        CalculateSkewness();
        CalculateAlphaValues();
        if(isi_threshold >= isi_cutoff){
            LN();
            DetectBurst();  
            CMADetectBurst();
        } 
    }
    
    
}

void Burst::CalculateSkewness(){
    double n = isi_sequence.size();
    double x_bar;
    double std_dev;
    double n_component = n/((n-1)*(n-2));
    double sum_component = 0;
    
    if(isi_sequence.empty()) {
        x_bar = 0;
        std_dev = 0;
    } else {
        x_bar = accumulate(isi_sequence.begin(), isi_sequence.end(), 0.0) / isi_sequence.size();

        for (int i; i < isi_sequence.size(); i++){
            std_dev += pow(isi_sequence.at(i)-x_bar, 2.0);
        }
        std_dev = std_dev/n;
        std_dev = sqrt(std_dev);
        
    }
    for (int i; i < isi_sequence.size(); i++){

        sum_component += (pow(((isi_sequence.at(i)-x_bar)/std_dev), 3.0));
    }
    isi_skewness_manual = n_component * sum_component;
}

void Burst::ProcessTree() {

}

void Burst::TreeToCSV() {

}

void Burst::SpkDetection(TH1F *V){
    struct Frame {
        float   x_1{ 0. };
        float   y_1{ 0. };
        float   x_2{ 0. };
        float   y_2{ 0. };
        float   m_1{ 0. };
        float   m_2{ 0. };
        int     spk_state { 0 }; // 0 = DEFAULT; 1 = UP; 2 = DOWN; 3 = DOWN_SPIKE;
        float   vm_max { 0. };
        float   vm_min { 0. };
        float   refractory_threshold { 1. };
        float   vm_peak_x { 1. }; //starting bin
        float   vm_peak_y { 0. };
        float   vm_peak_delta { 0. };

    } Spk_Frame;
    
    //Detect spikes
    int n_bins_Vm = V->GetNbinsX();
    double s_per_bin_Vm = V->GetBinWidth(0);
    Spk_Frame.vm_max = V->GetMaximum();
    Spk_Frame.vm_min = V->GetMinimum();
    Spk_Frame.refractory_threshold = 0.1 * (Spk_Frame.vm_max - Spk_Frame.vm_min);

    // TH1C* spk_hist_2;
    delete spike_train_V2;
    
    spike_train_V2 = new TH1C("spk_hist_2", "Spike Detection v2", 500 + n_bins_Vm * s_per_bin_Vm, 0, 500 + n_bins_Vm * s_per_bin_Vm); // check parameters to fit size
    TAxis* hist_axis_X = V->GetXaxis();
    spike_train_V2->GetXaxis()->SetRange(hist_axis_X->GetXmin(), hist_axis_X->GetXmax());

    Spk_Frame.vm_peak_y = V->GetBinContent(1); //initialize starting voltage

    for (int i = 1; i < n_bins_Vm; i++) {

        // Calculate gradient between current points (m_2)
        Spk_Frame.x_1 = i;
        Spk_Frame.y_1 = V->GetBinContent(i);
        Spk_Frame.x_2 = i + 1;
        Spk_Frame.y_2 = V->GetBinContent(i+1);
        Spk_Frame.m_2 = (Spk_Frame.y_2 - Spk_Frame.y_1)/(Spk_Frame.x_2 - Spk_Frame.x_1);

        // Check if gradient (m_2) is positive
        if(Spk_Frame.m_2 > Spk_Frame.m_1){

            // Set new peak of current spike
            Spk_Frame.vm_peak_x = Spk_Frame.x_2;
            Spk_Frame.vm_peak_y = Spk_Frame.y_2;

            Spk_Frame.spk_state = 1;

        } else { // else if (Spk_Frame.m_2 < Spk_Frame.m_1)
            
            //printf("\n declining... peak: %f to %f", Spk_Frame.vm_peak_y, Spk_Frame.y_2);
            // Check if spike has fell enough from the peak to reset
            Spk_Frame.vm_peak_delta = Spk_Frame.vm_peak_y - Spk_Frame.y_2;

            if (Spk_Frame.vm_peak_delta > Spk_Frame.refractory_threshold && Spk_Frame.spk_state != 3) {
                // Record spike
                //printf("\n [PEAK - %f, %f [delta - %f]", Spk_Frame.vm_peak_x, Spk_Frame.vm_peak_y, Spk_Frame.vm_peak_delta);
                
                // Record peak at previous vm_peak_x
                spike_train_V2->Fill(500 + Spk_Frame.vm_peak_x * s_per_bin_Vm);

                // Reset peak max and spike state
                Spk_Frame.vm_peak_x = Spk_Frame.x_2;
                Spk_Frame.vm_peak_y = Spk_Frame.y_2;
                Spk_Frame.spk_state = 3;
            }
        }
    } 


    // Set axis for spike_hist_v2 to match the rest
    
}

// Iterate through each spike in the train and check if it is valid.
void Burst::CleanSpikes(double window){    



    for (double spike : spikes_x){

        // Frame for scanning before and after spikes
        struct Frame {
            float   start_x { 0 };
            float   end_x { 0 };
            float   x_1{ 0. };
            float   y_1{ 0. };
            float   x_2{ 0. };
            float   y_2{ 0. };
            float   m{ 0. };
            float   m_max { 0. };

            Frame(int spike_x, TH1F* v_hist, double delta){
                
                // set start x & y
                x_1 = spike_x - delta;
                 if(x_1 < 1) {
                    x_1 = 1;
                }
                y_1 = v_hist->GetBinContent(v_hist->GetXaxis()->FindBin(x_1));

                // set end x & y
                x_2 = x_1 + (v_hist->GetBinWidth(0));
                if(x_2 > v_hist->GetNbinsX()) {
                    x_2 = v_hist->GetNbinsX();
                }
                y_2 =  v_hist->GetBinContent(v_hist->GetXaxis()->FindBin(x_2));

                start_x = x_1;
                end_x = start_x + (delta * 2);
            }
        };


        //cout << "\nSpike: " << spike;
        //cout << "\n    BinContent: " << spike_train_V2->GetBinContent(spike);
        //cout << "\n    BinsWidth: " << V->GetBinWidth(0) << " NBins: " << V->GetNbinsX();


        Frame Spk_Frame = Frame(spike, V, window);

        while (Spk_Frame.x_1 < Spk_Frame.end_x){

            //calculate m
            Spk_Frame.m = (Spk_Frame.y_2 - Spk_Frame.y_1) / (Spk_Frame.x_2 - Spk_Frame.x_1);
            if (Spk_Frame.m > Spk_Frame.m_max) {Spk_Frame.m_max = Spk_Frame.m;}

            //print
            //if(spike == 1350) {printf("\n    1 - (%f, %f) - (%f, %f) - m: %f, m_max: %f", Spk_Frame.x_1, Spk_Frame.y_1, Spk_Frame.x_2, Spk_Frame.y_2, Spk_Frame.m, Spk_Frame.m_max);}
            
            // Move frame 1 bin
            Spk_Frame.x_1   +=  V->GetBinWidth(0);
            Spk_Frame.y_1   =   V->GetBinContent(V->GetXaxis()->FindBin(Spk_Frame.x_1));

            Spk_Frame.x_2   =  Spk_Frame.x_1 + (V->GetBinWidth(0));
            Spk_Frame.y_2   =   V->GetBinContent(V->GetXaxis()->FindBin(Spk_Frame.x_2));
        }
    
        if(Spk_Frame.m_max < 20) {
            spike_train_V2->SetBinContent(spike, 0);
        }
    
    }
}


#endif

// working on 139 to get isi neighbour working
