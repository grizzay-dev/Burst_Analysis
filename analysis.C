
#define analysis_cxx
#include "analysis.h"
#include "burst.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <fstream>
#include <cmath>

void analysis::ProcessSingle(Long64_t entry, double k){
  
  if (fChain == 0) return;
  fChain->GetEntry(entry);

  //Draw spike train
  TCanvas *c1 = new TCanvas("canvas", "multipads", 1600, 600);
  c1->Divide(1,2,0,0);
  c1->cd(1);
  spike->Draw();
  c1->cd(2);
  Vm_1nA->Draw();
  
  //process spike train
  Burst b(entry, spike, k);
  b.ProcessNeuron();
  b.PrintMetrics();
  TLine *line;
  //Draw lines
  for(int i = 0; i < b.n_bursts; i++){
    int x1 = b.burst_locations[i][0];
    int x2 = b.burst_locations[i][1];

    line = new TLine(b.spikes_x[x1], 0.5, b.spikes_x[x2], 0.5);
    line -> SetLineWidth(10);
    line -> SetLineColorAlpha(kRed, 0.35);

    c1->cd(1);
    line->Draw();
    c1->cd(2);
    line->Draw();
    c1->cd(1);
  }
  c1->Print("burst.png");
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
  TFile *output = new TFile("burst.root", "recreate");
  //Parameters
  int     id,
          burst,
          nBursts;
  double  k,
          burstFreq,
          totalBurstDur,
          avgBurstDur,
          ML,
          isiThreshold;
  //Create tree
  TTree *tree = new TTree("burst_data", "Burst Data");
  //Add branches
  tree->Branch("id", &id, "id/I");
  tree->Branch("burst_type", &burst, "burst_type/I");
  tree->Branch("n_bursts", &nBursts, "n_bursts/I");
  tree->Branch("burst_total_duration", &totalBurstDur, "burst_total_duration/D");
  tree->Branch("burst_frequency", &burstFreq, "burst_frequency/D");
  tree->Branch("burst_average_duration", &avgBurstDur, "burst_average_duration/D");
  tree->Branch("ML", &ML, "ML/D");
  tree->Branch("isi_threshold", &isiThreshold, "isi_threshold/D");
  tree->Branch("k", &k, "k/D");   
  
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //Process spike train
    //TH1C* hist = (TH1C*)spike->Clone();
    Burst x(jentry, spike);
    x.ProcessNeuron(); 
    
    //Save to burst tree
    id              = x.id;
    k               = x.k;
    burst           = x.burst_type;
    nBursts         = x.n_bursts;
    burstFreq       = x.burst_frequency;
    totalBurstDur   = x.burst_total_duration;
    avgBurstDur     = x.burst_average_duration;
    ML              = x.ML;
    isiThreshold    = x.isi_threshold;
    
    tree->Fill();
  }

  //commit tree to output file
  output->Write();
  output->Close();


  
}
