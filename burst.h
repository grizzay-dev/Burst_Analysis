#ifndef burst_h
#define burst_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TH1.h"
#include "TH1.h"
#include "TGraph.h"
#include <cmath>
#include <vector>
#include <string>

using namespace std;

class Burst {
    //Method of burst detection adapated from: https://www.sciencedirect.com/science/article/pii/S1002007108003432
    public : 
        //Spike train
        TH1C* spike_train;// {nullptr};
        TH1* isi_histogram;

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

        double                  k{ 1 };
        double                  isi_sd{ 0 };
        double                  isi_mean{ 0 };
        double                  isi_cutoff{ 2 };
        double                  isi_threshold{ 0 };
        double                  ln_mean{ 0 };
        double                  isi_skewness{ 0 };

        vector<double>          spikes_x;
        vector<double>          isi_sequence;
        vector<double>          ln_sequence;

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
        Burst(int entry_id, TH1C* entry, double constant = 1){
            spike_train = entry;
            id = entry_id;
            k = constant;

            //set spike_train parameters
            seconds_per_bin = spike_train->GetBinWidth(0);
            n_bins = spike_train->GetNbinsX();

            TAxis* x_axis = spike_train->GetXaxis();
            injection_start = x_axis->GetXmin();
            injection_end = x_axis->GetXmax();
            total_recording_time = injection_end - injection_start;
        };

        virtual void    PrintMetrics();
        virtual void    ISIHistogram();
        virtual void    Spikes();
        virtual void    ISI();
        virtual void    LN();
        virtual void    DetectBurst();
        virtual void    ProcessNeuron();
};

void PrintSequence(vector<double> sequence, int length, string name = "Sequence:"){
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

void Burst::ISIHistogram() {
    // Only continue if atleast 2 spikes.
    if (n_spikes < 2) return;

    //Create histogram:
    double x_max = *max_element(isi_sequence.begin(), isi_sequence.end());
    x_max += (x_max * 0.1);  //add 10% space to end of x-axis
    double isi_bins = x_max / seconds_per_bin;

    isi_histogram = new TH1I("isi_histogram", "ISI Histogram", isi_bins, 0, x_max);
    
    //Fill histogram
    for (double isi : isi_sequence) {
        isi_histogram->Fill(isi);
    }
    
    //Save skewness
    isi_skewness = isi_histogram->GetSkewness();

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
    printf("\n");
    PrintSequence(spikes_x, n_spikes, "Spike Train");
    PrintSequence(isi_sequence, isi_length, "ISI Seq");
    
    if (n_bursts > 0) {
        int x = 1;
        printf("\nBurst Locations:");
        for (vector<int> i : burst_locations) {
            printf("\n%41s%-5i%.2f%6s   %.2f", "#", x, spikes_x.at(i.at(0)), "--->", spikes_x.at(i.at(1)));
            x++;
        }
            
    }
}

void Burst::Spikes() {
    for (int bin = 0; bin <= n_bins; bin++) {
        double bin_value = spike_train->GetBinContent(bin);
        if (bin_value > 0) {
            n_spikes += 1;
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
    cout << "\nProcessing entry: " << id;
    Spikes();
    
    if(n_spikes > 0){
        ISI();
        ISIHistogram();
        if(isi_threshold >= isi_cutoff){
            LN();
            DetectBurst();  
        } 
    }
}

#endif