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

class Burst {
    //Method of burst detection adapated from: https://www.sciencedirect.com/science/article/pii/S1002007108003432
    public : 
        //Spike train
        TH1C        *spike_train;

        //Burst processing properties
        int         id;
        int         n_bins;
        int         n_spikes                = 0;
        int         isi_length;
        int         ln_length;

        double      k;
        double*     spikes_x;
        double*     isi_sequence;
        double*     ln_sequence;
        double      isi_sd                  = 0;
        double      isi_mean                = 0;
        double      isi_threshold           = 0;
        double      ln_mean                 = 0;

        //Burst metrics
        bool        bursting                = false;
        int         burst_type              = 0;
        int         n_bursts                = 0;
        int**       burst_locations;

        double      burst_total_duration    = 0;
        double      burst_average_duration  = 0;
        double      burst_frequency         = 0;
        double      ML                      = 0;

        //Constructor
        Burst(int entry_id, TH1C* entry, double constant = 1){
            spike_train = entry;
            id = entry_id;
            k = constant;
        };

        virtual void    PrintMetrics();
        virtual void    Spikes();
        virtual void    ISI();
        virtual void    LN();
        virtual void    DetectBurst();
        virtual void    BurstLocations();
        virtual void    ProcessNeuron();
};

void PrintSequence(double sequence[], int length, string name = "Sequence:"){
    for(int i = 0; i < length; i++){
        if(i == 0){
            cout << name << ": ";
        }
        cout << sequence[i];
        if(i != length -1) cout << ", ";
    }
    cout << " - END\n";
}

double CumulativeAverage(int start, int end, double data[])
{
  double sum = 0.0, delta, avg;
  for(int x = start; x < end; x++){
    delta = data[x+1] - data[x];
    sum += delta;
  }
  avg = sum / (end-start);
  return avg;
}

void Burst::PrintMetrics(){

    cout << "\nBurst Metrics:\n";
    cout << "Entery (id): "                 << id                       << ".\n";
    cout << "Burst type: "                  << burst_type               << ".\n";
    cout << "# of bursts: "                 << n_bursts                 << ".\n";
    if(n_bursts > 0){
        for(int z = 0; z < n_bursts; z++){
            cout << "    #" << (z+1) << ": " << burst_locations[z][0] << " ---> " << burst_locations[z][1] << "\n";
        }
    }
    cout << "Total burst duration: "        << burst_total_duration     << ".\n";
    cout << "Average burst duration: "      << burst_average_duration   << ".\n";
    cout << "Burst frequency: "             << burst_frequency          << ".\n";
    cout << "ML: "                          << ML                       << ".\n";
    cout << "ISI threshold: "               << isi_threshold            << ".\n";
    PrintSequence(spikes_x, n_spikes, "Spike Train");
    PrintSequence(isi_sequence, isi_length, "ISI Seq");
    
}

void Burst::Spikes(){
    n_bins = spike_train->GetNbinsX();
    //Set spike count to 0
    n_spikes = 0;
    for(int bin = 0; bin < n_bins; bin++){
        int spike_value;
        spike_value = spike_train->GetBinContent(bin);
        if(spike_value == 1) n_spikes ++;
    }   
    
    if(n_spikes > 0){
        //assign array length
        spikes_x = (double*) malloc(sizeof(double)*n_spikes);

        //fill spikes_x with spike x position
        int counter = 0;
        for(int bin = 0; bin < n_bins; bin++){
            int spike_value;
            spike_value = spike_train->GetBinContent(bin);
            if(spike_value==1){
                spikes_x[counter] = (bin * 0.25) + 400; //0.25s per bin, recordings start at 400ms.
                counter++;
            }
        }
    }
}

void Burst::ISI(){
    if(n_spikes==0) return;

    isi_length = n_spikes - 1;
    
    //resize isi_sequence
    isi_sequence = (double*) malloc(sizeof(double)*isi_length);
    for(int x = 0; x < isi_length; x++){
        isi_sequence[x] = spikes_x[x+1] - spikes_x[x];
    }

    double sum = 0.0;
    for(int x = 0; x < isi_length; x++){
        sum += isi_sequence[x];
    }

    isi_mean = sum / isi_length;
    for(int x = 0; x < isi_length; ++x){
        isi_sd += pow(isi_sequence[x] - isi_mean, 2);
    }
    isi_sd = sqrt(isi_sd/double(isi_length));
    isi_threshold = isi_sd / isi_mean;
}

void Burst::LN(){
    ln_length = 0;
    for(int x = 0; x < isi_length; x++){
        if(isi_sequence[x] < isi_mean) ln_length++;
    }

    //resize ln_sequence
    ln_sequence = (double*) malloc(sizeof(double)*ln_length);
    int ln_counter = 0;
    for(int x = 0; x < isi_length; x++){
        if(isi_sequence[x] < isi_mean){
            ln_sequence[ln_counter] = isi_sequence[x];
            ln_counter++;
        }
    }

    double ln_sum = 0;
    for(int x = 0; x < ln_length; x++){
        ln_sum += ln_sequence[x];
    }

    ln_mean = (double)ln_sum / (double)ln_length;
}

void Burst::DetectBurst(){
    
    int start = 0;
    int end =  2;
    bool detected = false;

    while (end < isi_length + 2){
        double burst_average = CumulativeAverage(start, end, spikes_x);
        burst_average = CumulativeAverage(start, end, spikes_x);
        if((burst_average < (ln_mean * k)) && (burst_average > 0)){
            detected = true;
            end += 1;
            bursting = true;
            burst_type = 1;
        } else {
            if(detected){
                n_bursts += 1;
                burst_total_duration += (spikes_x[end-1] - spikes_x[start]);

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
        burst_frequency = double(n_bursts) / 600.0; //current injection from 400ms to 600ms 
    }
}

void Burst::BurstLocations(){
    if(n_bursts > 0){
        //burst_locations = (int**) malloc((n_bursts * 2) * sizeof(int)); THIS DOESNT WORK SEE BELOW
        //NEED TO CHECK THIS malloc() for validity
        burst_locations = (int**) malloc (sizeof(int*) * n_bursts * 2);
        for(int i = 0; i < n_bursts; i++){
            burst_locations[i] = (int*) malloc(sizeof(int*) * 2);
        }

        //reusing burst detect code not optomised...
        int start = 0;
        int end =  2;
        int burst_counter = 0;
        bool detected = false;

        while (end < isi_length + 2){
            double burst_average = CumulativeAverage(start, end, spikes_x);

            if((burst_average < (ln_mean * k)) && (burst_average > 0)){
                detected = true;
                end += 1;
            } else {
                if(detected){
                    //save locations
                    burst_locations[burst_counter][0] = start;
                    burst_locations[burst_counter][1] =  end-1;
                    burst_counter++;
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
    }
}

void Burst::ProcessNeuron(){
    cout << "\nProcessing entry: " << id;
    
    Spikes();
    
    if(n_spikes > 0){
        ISI();
        if(isi_threshold >= 2){
            LN();
            DetectBurst();
            BurstLocations();    
        } 
    }
}

#endif