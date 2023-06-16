#define analysis_cxx
#include "analysis.h"
#include "burst.h"

#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <fstream>
#include <cmath>

    void analysis::ProcessSingle(Long64_t entry, double k = 1) {

    if (fChain == 0) return;
    fChain->GetEntry(entry);

    //process spike train
    Burst b(entry, spike, Vm_1nA, k);

    b.ProcessNeuron();

    //Draw spike train
    TCanvas* c1 = new TCanvas("canvas", "multipads", 1600, 800);
    c1->Divide(1, 3, 0, 0);

    //https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html
    // 
    //TOP - isi_histogram
    c1->cd(1);
    gStyle->SetOptStat("s");
    c1->SetLogy();
    b.isi_histogram->GetYaxis()->SetTickLength(-0.01);
    //b.isi_histogram->Fit("gaus", "WW", "");
    gPad->SetBottomMargin(1);
    b.isi_histogram->Draw();


    //CENTRE - spike_train
    c1->cd(2);
    double bin_max = b.spike_train_V2->GetMaximum();

    b.spike_train_V2->GetYaxis()->SetRangeUser(0, (2 * 1.5));
    b.spike_train_V2->GetYaxis()->SetTickLength(-0.01);
    gPad->SetBottomMargin(1);
    b.spike_train_V2->Draw();


    //BOTTOM - Vm_1nA
    c1->cd(3);
    Vm_1nA->GetYaxis()->SetTickLength(-0.01);
    Vm_1nA->Draw();


    if (b.n_spikes > 0) {
        b.PrintMetrics();
        TLine* line;
        
        //Draw lines
        for (int i = 0; i < b.n_bursts; i++) {

            int x1 = b.burst_locations.at(i).at(0);
            int x2 = b.burst_locations.at(i).at(1);

            line = new TLine(b.spikes_x.at(x1), (2 * 1.25), b.spikes_x.at(x2), (2 * 1.25));
            
            line->SetLineWidth(10);
            line->SetLineColorAlpha(kRed, 0.1);

            c1->cd(2);
            line->Draw();
            

        }

    if (b.cma_n_bursts > 0){
        TLine* cma_line;
        for (int i = 0; i < b.cma_n_bursts; i++){
            int cma_x1 = b.cma_burst_locations.at(i).at(0);
            int cma_x2 = b.cma_burst_locations.at(i).at(1);

            cma_line = new TLine(b.spikes_x.at(cma_x1), (1.5), b.spikes_x.at(cma_x2), (1.5));

            cma_line->SetLineWidth(10);
            cma_line->SetLineColorAlpha(kBlue, 0.1);
            c1->cd(2);
            cma_line->Draw();
        }

        for (int i = 0; i < b.cma_neighbour_locations.size(); i++){
            TMarker* marker;
            marker = new TMarker(b.cma_neighbour_locations.at(i), 2, 3);
            c1->cd(2);
            marker->Draw();
        }
    }

        printf("\n\n Procesing complete");
        c1->Print("burst.png");
    }
    else {
        printf("\nNo spikes detected");
    }

}

void analysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analysis.C
//      root> analysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   //Create root file for burst data tree
    TFile* output = new TFile("burst_p4_1nA.root", "recreate");
    //Parameters
    int     id,
            burst,
            nBursts,
            cmaNBurst,
            nSpikes;
    double  k,
            burstFreq,
            totalBurstDur,
            avgBurstDur,
            ML,
            isiThreshold,
            isiSkewness;
    
    //Create tree
    TTree* tree = new TTree("burst_data", "Burst Data");
    //Add branches
    tree->Branch("id", &id, "id/I");
    tree->Branch("burst_type", &burst, "burst_type/I");
    tree->Branch("n_bursts", &nBursts, "n_bursts/I");
    tree->Branch("cma_n_bursts", &cmaNBurst, "cma_n_bursts/I");
    tree->Branch("burst_total_duration", &totalBurstDur, "burst_total_duration/D");
    tree->Branch("burst_frequency", &burstFreq, "burst_frequency/D");
    tree->Branch("burst_average_duration", &avgBurstDur, "burst_average_duration/D");
    tree->Branch("ML", &ML, "ML/D");
    tree->Branch("isi_threshold", &isiThreshold, "isi_threshold/D");
    tree->Branch("isi_skewness", &isiSkewness, "isi_skewness/D");
    tree->Branch("k", &k, "k/D");
    tree->Branch("n_spikes", &nSpikes, "n_spikes/I");


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //Process spike train
      //TH1C* hist = (TH1C*)spike->Clone();
      Burst x(jentry, spike, Vm_1nA);

      x.ProcessNeuron();

      //Save to burst tree
      id = x.id;
      k = x.k;
      burst = x.burst_type;
      nBursts = x.n_bursts;
      cmaNBurst = x.cma_n_bursts;
      burstFreq = x.burst_frequency;
      totalBurstDur = x.burst_total_duration;
      avgBurstDur = x.burst_average_duration;
      ML = x.ML;
      isiThreshold = x.isi_threshold;
      isiSkewness = x.isi_skewness;
      nSpikes = x.n_spikes;

      tree->Fill();
   }

   //commit tree to output file
    output->Write();
    output->Close();


    TTreeToCSV();

}


void analysis::TTreeToCSV() {
    //TTree
    TFile* F = new TFile("LCMV5_P4_n100000_001.root");
    TTree* T = (TTree*)F->Get("neur_test_large");
    T->AddFriend("burst_data", "burst_p4_1nA.root");

    //Output CSV
    ofstream csv_output;
    csv_output.open("ttree_output_p4_1nA.csv");

    //TTree leafs
    int             id;
    Float_t         K_D_gmax;
    Float_t         KV3_1_gmax;
    Float_t         KV2_FAST_gmax;
    Float_t         KV2_SLOW_gmax;
    Float_t         KV1_4_gmax;
    Float_t         KV4_2_gmax;
    Float_t         K_M_gmax;
    Float_t         SK_gmax;
    Float_t         NA_T_AX_gmax;
    Float_t         NA_T_SD_gmax;
    Float_t         NA_P_gmax;
    Float_t         CA_LVA_gmax;
    Float_t         CA_HVA_gmax;
    Float_t         LEAK_gmax;
    Float_t         IH_gmax;
    int             n_spikes;
    int             burst_type;
    int             cma_n_bursts;
    double          isi_skewness;
    int             cma_burst;
    double          spike_freq;

    

    T->SetBranchAddress("id", &id);
    T->SetBranchAddress("K_D_gmax", &K_D_gmax);
    T->SetBranchAddress("KV3_1_gmax", &KV3_1_gmax);
    T->SetBranchAddress("KV2_FAST_gmax", &KV2_FAST_gmax);
    T->SetBranchAddress("KV2_SLOW_gmax", &KV2_SLOW_gmax);
    T->SetBranchAddress("KV1_4_gmax", &KV1_4_gmax);
    T->SetBranchAddress("KV4_2_gmax", &KV4_2_gmax);
    T->SetBranchAddress("K_M_gmax", &K_M_gmax);
    T->SetBranchAddress("SK_gmax", &SK_gmax);
    T->SetBranchAddress("NA_T_AX_gmax", &NA_T_AX_gmax);
    T->SetBranchAddress("NA_T_SD_gmax", &NA_T_SD_gmax);
    T->SetBranchAddress("NA_P_gmax", &NA_P_gmax);
    T->SetBranchAddress("CA_LVA_gmax", &CA_LVA_gmax);
    T->SetBranchAddress("CA_HVA_gmax", &CA_HVA_gmax);
    T->SetBranchAddress("LEAK_gmax", &LEAK_gmax);
    T->SetBranchAddress("IH_gmax", &IH_gmax);
    T->SetBranchAddress("n_spikes", &n_spikes);
    T->SetBranchAddress("burst_type", &burst_type);
    T->SetBranchAddress("isi_skewness", &isi_skewness);
    T->SetBranchAddress("cma_n_bursts", &cma_n_bursts);


    //SET HEADINGS for outptu
    csv_output << "id,K_D_gmax,KV3_1_gmax,KV2_FAST_gmax,KV1_4_gmax,KV4_2_gmax,K_M_gmax,SK_gmax,NA_T_AX_gmax,NA_T_SD_gmax,NA_P_gmax,CA_LVA_gmax,CA_HVA_gmax,IH_gmax,n_spikes,spike_freq,isi_skewness,cma_n_bursts,cma_burst\n";
    
    //Loop through entries and write to file
    Long64_t nentries = T->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t ientry = 0; ientry < nentries; ientry++)
    {
        T->GetEntry(ientry);

        spike_freq = n_spikes / 5.0;

        if(cma_n_bursts>0){
            cma_burst = 1;
        } else {
            cma_burst = 0;
        }

        //output row to csv
        csv_output << id << "," << K_D_gmax << "," << KV3_1_gmax << "," << KV2_FAST_gmax << "," << KV1_4_gmax << "," << KV4_2_gmax << ",";
        csv_output << K_M_gmax << "," << SK_gmax << "," << NA_T_AX_gmax << "," << NA_T_SD_gmax << "," << NA_P_gmax << "," << CA_LVA_gmax << "," << CA_HVA_gmax << ",";
        csv_output << IH_gmax << "," << n_spikes << "," << spike_freq << "," << isi_skewness << "," << cma_n_bursts << "," << cma_burst;

        //complete line
        csv_output << endl;

    }
    csv_output.close();
}
