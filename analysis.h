//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 23 10:11:28 2022 by ROOT version 6.26/02
// from TTree neur_test_large/[output] results of neuron test
// found on file: ./neuron_P4.root
//////////////////////////////////////////////////////////

#ifndef analysis_h
#define analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TH1.h"
#include "TH1.h"
#include "TGraph.h"

class analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
   Float_t         spk_width[30];
   Float_t         spk_width_sd[30];
   Float_t         spk_height[30];
   Float_t         spk_height_sd[30];
   Float_t         spk_dVmdt[30];
   Float_t         spk_dVmdt_sd[30];
   Float_t         Vm_min[30];
   Float_t         Vm_min_sd[30];
   Float_t         ISI[30];
   Float_t         ISI_sd[30];
   Float_t         tm2spk[30];
   Int_t           spk_num[30];
   Int_t           simu_state;
   Int_t           FR_0;
   Float_t         FI_rheo;
   Float_t         FI_rheo_fit;
   Float_t         FI_slope;
   Float_t         FI_slope_err;
   TH1F            *Vm_0nA;
   TH1F            *Vm_0p1nA;
   TH1F            *Vm_0p2nA;
   TH1F            *Vm_0p5nA;
   TH1F            *Vm_1nA;
   TH1F            *Vm_soma;
   TH1C            *spike;
   TGraph          *FI;

   // List of branches
   TBranch        *b_K_D_gmax;   //!
   TBranch        *b_KV3_1_gmax;   //!
   TBranch        *b_KV2_FAST_gmax;   //!
   TBranch        *b_KV2_SLOW_gmax;   //!
   TBranch        *b_KV1_4_gmax;   //!
   TBranch        *b_KV4_2_gmax;   //!
   TBranch        *b_K_M_gmax;   //!
   TBranch        *b_SK_gmax;   //!
   TBranch        *b_NA_T_AX_gmax;   //!
   TBranch        *b_NA_T_SD_gmax;   //!
   TBranch        *b_NA_P_gmax;   //!
   TBranch        *b_CA_LVA_gmax;   //!
   TBranch        *b_CA_HVA_gmax;   //!
   TBranch        *b_LEAK_gmax;   //!
   TBranch        *b_IH_gmax;   //!
   TBranch        *b_spk_width;   //!
   TBranch        *b_spk_width_sd;   //!
   TBranch        *b_spk_height;   //!
   TBranch        *b_spk_height_sd;   //!
   TBranch        *b_spk_dVmdt;   //!
   TBranch        *b_spk_dVmdt_sd;   //!
   TBranch        *b_Vm_min;   //!
   TBranch        *b_Vm_min_sd;   //!
   TBranch        *b_ISI;   //!
   TBranch        *b_ISI_sd;   //!
   TBranch        *b_tm2spk;   //!
   TBranch        *b_spk_num;   //!
   TBranch        *b_simu_state;   //!
   TBranch        *b_FR_0;   //!
   TBranch        *b_FI_rheo;   //!
   TBranch        *b_FI_rheo_fit;   //!
   TBranch        *b_FI_slope;   //!
   TBranch        *b_FI_slope_err;   //!
   TBranch        *b_Vm_0nA;   //!
   TBranch        *b_Vm_0p1nA;   //!
   TBranch        *b_Vm_0p2nA;   //!
   TBranch        *b_Vm_0p5nA;   //!
   TBranch        *b_Vm_1nA;   //!
   TBranch        *b_Vm_soma;   //!
   TBranch        *b_spike;   //!
   TBranch        *b_FI;   //!

   analysis(TTree *tree=0);
   virtual ~analysis();
   virtual void     TTreeToCSV();
   virtual void     ProcessSingle(Long64_t entry, double k);
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_cxx
analysis::analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./LCMV5_P4_n100000_001.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./LCMV5_P4_n100000_001.root");
      }
      f->GetObject("neur_test_large",tree);

   }
   Init(tree);
}

analysis::~analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Vm_0nA = 0;
   Vm_0p1nA = 0;
   Vm_0p2nA = 0;
   Vm_0p5nA = 0;
   Vm_1nA = 0;
   Vm_soma = 0;
   spike = 0;
   FI = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("K_D_gmax", &K_D_gmax, &b_K_D_gmax);
   fChain->SetBranchAddress("KV3_1_gmax", &KV3_1_gmax, &b_KV3_1_gmax);
   fChain->SetBranchAddress("KV2_FAST_gmax", &KV2_FAST_gmax, &b_KV2_FAST_gmax);
   fChain->SetBranchAddress("KV2_SLOW_gmax", &KV2_SLOW_gmax, &b_KV2_SLOW_gmax);
   fChain->SetBranchAddress("KV1_4_gmax", &KV1_4_gmax, &b_KV1_4_gmax);
   fChain->SetBranchAddress("KV4_2_gmax", &KV4_2_gmax, &b_KV4_2_gmax);
   fChain->SetBranchAddress("K_M_gmax", &K_M_gmax, &b_K_M_gmax);
   fChain->SetBranchAddress("SK_gmax", &SK_gmax, &b_SK_gmax);
   fChain->SetBranchAddress("NA_T_AX_gmax", &NA_T_AX_gmax, &b_NA_T_AX_gmax);
   fChain->SetBranchAddress("NA_T_SD_gmax", &NA_T_SD_gmax, &b_NA_T_SD_gmax);
   fChain->SetBranchAddress("NA_P_gmax", &NA_P_gmax, &b_NA_P_gmax);
   fChain->SetBranchAddress("CA_LVA_gmax", &CA_LVA_gmax, &b_CA_LVA_gmax);
   fChain->SetBranchAddress("CA_HVA_gmax", &CA_HVA_gmax, &b_CA_HVA_gmax);
   fChain->SetBranchAddress("LEAK_gmax", &LEAK_gmax, &b_LEAK_gmax);
   fChain->SetBranchAddress("IH_gmax", &IH_gmax, &b_IH_gmax);
   fChain->SetBranchAddress("spk_width", spk_width, &b_spk_width);
   fChain->SetBranchAddress("spk_width_sd", spk_width_sd, &b_spk_width_sd);
   fChain->SetBranchAddress("spk_height", spk_height, &b_spk_height);
   fChain->SetBranchAddress("spk_height_sd", spk_height_sd, &b_spk_height_sd);
   fChain->SetBranchAddress("spk_dVmdt", spk_dVmdt, &b_spk_dVmdt);
   fChain->SetBranchAddress("spk_dVmdt_sd", spk_dVmdt_sd, &b_spk_dVmdt_sd);
   fChain->SetBranchAddress("Vm_min", Vm_min, &b_Vm_min);
   fChain->SetBranchAddress("Vm_min_sd", Vm_min_sd, &b_Vm_min_sd);
   fChain->SetBranchAddress("ISI", ISI, &b_ISI);
   fChain->SetBranchAddress("ISI_sd", ISI_sd, &b_ISI_sd);
   fChain->SetBranchAddress("tm2spk", tm2spk, &b_tm2spk);
   fChain->SetBranchAddress("spk_num", spk_num, &b_spk_num);
   fChain->SetBranchAddress("simu_state", &simu_state, &b_simu_state);
   fChain->SetBranchAddress("FR_0", &FR_0, &b_FR_0);
   fChain->SetBranchAddress("FI_rheo", &FI_rheo, &b_FI_rheo);
   fChain->SetBranchAddress("FI_rheo_fit", &FI_rheo_fit, &b_FI_rheo_fit);
   fChain->SetBranchAddress("FI_slope", &FI_slope, &b_FI_slope);
   fChain->SetBranchAddress("FI_slope_err", &FI_slope_err, &b_FI_slope_err);
   fChain->SetBranchAddress("Vm_0nA", &Vm_0nA, &b_Vm_0nA);
   fChain->SetBranchAddress("Vm_0p1nA", &Vm_0p1nA, &b_Vm_0p1nA);
   fChain->SetBranchAddress("Vm_0p2nA", &Vm_0p2nA, &b_Vm_0p2nA);
   fChain->SetBranchAddress("Vm_0p5nA", &Vm_0p5nA, &b_Vm_0p5nA);
   fChain->SetBranchAddress("Vm_1nA", &Vm_1nA, &b_Vm_1nA);
   fChain->SetBranchAddress("Vm_soma", &Vm_soma, &b_Vm_soma);
   fChain->SetBranchAddress("spike", &spike, &b_spike);
   fChain->SetBranchAddress("FI", &FI, &b_FI);
   Notify();
}

Bool_t analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_cxx
