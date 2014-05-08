#include "SSLeptonAnalysis/interface/SSLeptonAnalysis.h"
#include <iostream>

#define SSLeptonAnalysisDEBUG 0

using namespace std;

void SSLeptonAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{}

// ----------------------------------------------------------------------------------------------------
SSLeptonAnalysis::SSLeptonAnalysis()  : 
  name_("SSLeptonAnalysis")
{}

// ----------------------------------------------------------------------------------------------------
SSLeptonAnalysis::~SSLeptonAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::Init(LoopAll& l) {
  
  if (l.typerun != l.kReduce) {
    if (puHist != "" && puHist != "auto" ) {
      if(DEBUG)
        cout << "Opening PU file"<<endl;
      TFile* puFile = TFile::Open( puHist );
      if (puFile) {
        TH1 * target = 0;
	
        if( puTarget != "" ) {
          TFile * puTargetFile = TFile::Open( puTarget );
          assert( puTargetFile != 0 );
          target = (TH1*)puTargetFile->Get("pileup");
          if( target == 0 ) { target = (TH1*)puTargetFile->Get("target_pu"); }
          target->Scale( 1. / target->Integral() );
        }
	
        if( puMap != "" ) {
          loadPuMap(puMap, puFile, target);
        } else {
	  loadPuWeights(0, puFile, target);
        }
        puFile->Close();
      }
      else {
        cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl;
        weights[0].resize(50);
        for (unsigned int i=0; i<weights[0].size(); i++) weights[0][i] = 1.0;
      }
      if(DEBUG)
        cout << "Opening PU file END"<<endl;
    } else if ( puHist == "auto" ) {
      TFile * puTargetFile = TFile::Open( puTarget );
      assert( puTargetFile != 0 );
      puTargetHist = (TH1*)puTargetFile->Get("pileup");
      if( puTargetHist == 0 ) {
        puTargetHist = (TH1*)puTargetFile->Get("target_pu");
      }
      puTargetHist = (TH1*)puTargetHist->Clone();
      puTargetHist->SetDirectory(0);
      puTargetHist->Scale( 1. / puTargetHist->Integral() );
      puTargetFile->Close();
    }
  }
  
  
  if (l.typerun == l.kReduce)
    muCorrector_ = new MuScleFitCorrector(muFitParametersFile);
  
  hltSelection.push_back("HLT_Mu17_Mu8_v*");
  hltSelection.push_back("HLT_IsoMu24_v*");
  hltSelection.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  hltSelection.push_back("HLT_DoubleMu8_Mass8_v*");
  hltSelection.push_back("HLT_DoubleMu14_Mass8_v*");
  hltSelection.push_back("HLT_DoubleEle8_CaloIdT_TrkIdVL_Mass_v*");
  hltSelection.push_back("HLT_DoubleEle14_Mass8_CaloIdT_TrkIdVL_pfMHT40_v*");
  hltSelection.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*");
  hltSelection.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoV_v*");
  hltSelection.push_back("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Mass*_v*");
  hltSelection.push_back("HLT_Mu14_Ele14_Mass8_CaloIdT_CaloIsoVL_TrkIdVL_pfMHT40_v*");
}

int SSLeptonAnalysis::categories(TLorentzVector* p1, TLorentzVector* p2, int mixed) {

  if (mixed == 1) {
    if (fabs(p1->Eta())<1.479 and fabs(p2->Eta()) < 1.479)
      return 0;
    else
      return 1;
  } else if (mixed == 0){
    if (fabs(p1->Eta())<0.8 and fabs(p2->Eta())<0.8)
      return 0;
    else
      return 1;
  } else {
    if (fabs(p1->Eta()) < 1.479 and fabs(p2->Eta())<1.479)
      return 0;
    else if (fabs(p1->Eta()) < 1.479)
      return 1;
    else
      return 2;
  }
}

bool SSLeptonAnalysis::removeZLeptons(LoopAll& l, int l1, bool isMuon) {

  if (isMuon) {
    TLorentzVector* p1 = (TLorentzVector*)(l.mu_glo_p4_corr->At(l1));
    for (unsigned int i=0; i<l.mu_glo_n; i++) {
      if (i == l1)
	continue;
      if (l.mu_glo_charge[l1] != l.mu_glo_charge[i]) {
	TLorentzVector* p2 = (TLorentzVector*)(l.mu_glo_p4_corr->At(i));
	float temp_mass = (*p1+*p2).M();
	if (temp_mass > 91.186-10. or temp_mass < 91.186+10)
	  return true;
      }
    }
    return false;
  } else {
    TLorentzVector* p1 = (TLorentzVector*)(l.el_std_p4_corr->At(l1));
    for (unsigned int i=0; i<l.el_std_n; i++) {
      if (i == l1)
	continue;
      if (l.el_std_charge[l1] != l.el_std_charge[i]) {
	TLorentzVector* p2 = (TLorentzVector*)(l.el_std_p4_corr->At(i));
	float temp_mass = (*p1+*p2).M();
	if (temp_mass > 91.186-10. or temp_mass < 91.186+10)
	  return true;
      }
    }
  }

  return false;
}


// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::Analysis(LoopAll& l, Int_t jentry) {

  //apply pileup reweighting
  float weight = 1.;
  float pu_weight = 1.;
  int cur_type = l.itype[l.current];
  
  if (cur_type != 0) {
    unsigned int n_pu = l.pu_n;
    pu_weight = 1;//getPuWeight(l.pu_n, cur_type, &(l.sampleContainer[l.current_sample_index]), jentry == 1);
    weight = pu_weight * l.sampleContainer[l.current_sample_index].weight();
  }

  std::vector<unsigned int> goodEl, goodMu;
  
  for (int i=0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)(l.el_std_p4_corr->At(i));
    
    if (ElectronMVACuts(l, i) and p4->Et() > 25. and !removeZLeptons(l, i, false))
      goodEl.push_back(i);
  }

  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)(l.mu_glo_p4_corr->At(i));
    if ((muIDTight and l.mu_glo_id_tight[i] == 1 and p4->Pt() > 25.) or
	(!muIDTight and l.mu_glo_id_loose[i] == 1 and p4->Pt() > 25.))
      if (!removeZLeptons(l, i, true))
	goodMu.push_back(i);
  }

  Float_t mass[100];
  Int_t type[100], cat[100];
  Int_t index1[100], index2[100];
  Int_t pairs = 0;
  Float_t sumPt[100];

  if (goodMu.size() != 0) {
    for (unsigned int i=0; i<goodMu.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)(l.mu_glo_p4_corr->At(goodMu[i]));
      for (unsigned int j=i+1; j<goodMu.size(); j++) {
	if (l.mu_glo_charge[goodMu[i]] == l.mu_glo_charge[goodMu[j]]) {
	  TLorentzVector* p2 = (TLorentzVector*)(l.mu_glo_p4_corr->At(goodMu[j]));
	  
	  sumPt[pairs] = p1->Pt() + p2->Pt();
	  mass[pairs] = (*p1+*p2).M();
	  type[pairs] = 0;
	  cat[pairs] = categories(p1, p2, 0);
	  index1[pairs] = goodMu[i];
	  index2[pairs] = goodMu[j];
	  pairs++;
	}
      }
    }
  }

  if (goodEl.size() != 0) {	
    for (unsigned int i=0; i<goodEl.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)(l.el_std_p4_corr->At(goodEl[i]));
      for (unsigned int j=i+1; j<goodEl.size(); j++) {
	if (l.el_std_charge[goodEl[i]] == l.el_std_charge[goodEl[j]]) {
	  TLorentzVector* p2 = (TLorentzVector*)(l.el_std_p4_corr->At(goodEl[j]));
	  sumPt[pairs] = p1->Pt() + p2->Pt();
	  mass[pairs] = (*p1+*p2).M();
	  type[pairs] = 1;
	  cat[pairs] = categories(p1, p2, 1);
	  index1[pairs] = goodEl[i];
	  index2[pairs] = goodEl[j];
	  pairs++;
	}
      }
    }
  }

  if (goodMu.size() != 0 and goodEl.size() != 0) {
    for (unsigned int i=0; i<goodEl.size(); i++) {
      TLorentzVector* p1 = (TLorentzVector*)(l.el_std_p4_corr->At(goodEl[i]));
      for (unsigned int j=0; j<goodMu.size(); j++) {
	if (l.el_std_charge[goodEl[i]] == l.mu_glo_charge[goodMu[j]]) {
	  TLorentzVector* p2 = (TLorentzVector*)(l.mu_glo_p4_corr->At(goodMu[j]));
	  sumPt[pairs] = p1->Pt() + p2->Pt();
	  mass[pairs] = (*p1+*p2).M();
	  type[pairs] = 2;
	  cat[pairs] = categories(p1, p2, 2);
	  index1[pairs] = goodEl[i];
	  index2[pairs] = goodMu[j];
	  pairs++;
	}
      }
    }
  }

  if (pairs > 0) {
    Tree(l, pairs, type, sumPt, mass, cat, index1, index2, weight, pu_weight);
    return true;
  }

  return false;
}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::FillReductionVariables(LoopAll& l, int jentry) {

  for (int i=0; i<l.el_std_n; i++) {
    l.el_std_ch_ctf[i] = l.tk_charge[l.el_std_tkind[i]];
  }

  l.mu_glo_p4_corr->Clear();
  
  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p = (TLorentzVector*)(l.mu_glo_p4->At(i)->Clone());
    if (l.itype[l.current] == 0)
      muCorrector_->applyPtCorrection(*p, l.mu_glo_charge[i]);
    else
      muCorrector_->applyPtSmearing(*p, l.mu_glo_charge[i]);
    
    new((*(l.mu_glo_p4_corr))[i]) TLorentzVector();
    ((TLorentzVector*)l.mu_glo_p4_corr->At(i))->SetXYZM(p->X(), p->Y(), p->Z(), p->M());

    if (l.MuonLooseID2012(i) && l.MuonIsolation2012(i, p->Pt(), false))
      l.mu_glo_id_loose[i] = 1;
    else
      l.mu_glo_id_loose[i] = 0;

    if (l.MuonTightID2012(i) && l.MuonIsolation2012(i, p->Pt(), true))
      l.mu_glo_id_tight[i] = 1;
    else
      l.mu_glo_id_tight[i] = 0;
  }

  l.doJetMatching(*l.jet_algoPF1_p4,*l.genjet_algo1_p4,l.jet_algoPF1_genMatched,l.jet_algoPF1_vbfMatched,l.jet_algoPF1_bgenMatched,l.jet_algoPF1_cgenMatched,l.jet_algoPF1_lgenMatched,l.jet_algoPF1_genPt,l.jet_algoPF1_genDr);
  MetCorrections2012(l);

//  if (l.itype[l.current] < 0) {
//    for (int i=0; i<l.gp_n; i++) {
//      if (l.gp_pdgid[i] == 25 && l.gp_status[i] == 3) {
//        TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(i);
//        *(l.higgs) = *p4;
//        break;
//      }
//    }
//    for (int i=0; i<l.mu_glo_n; i++) {
//      Int_t mc1=-1;
//      Int_t mc2=-1;
//      Int_t pho=-1;
//      
//      l.FindMCLeptons(i, mc1, mc2, pho, 13);
//      // et, eta, phi
//      
//      if (mc1 != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(mc1);
//	l.mc_et[i] = p4->Et();
//	l.mc_phi[i] = p4->Phi();
//	l.mc_eta[i] = p4->Eta();
//      } else {
//	l.mc_et[i] = -999.;
//	l.mc_phi[i] = -999.;
//	l.mc_eta[i] = -999.;
//      }
//      
//      if (pho != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(pho);
//	l.fsr_et[i] = p4->Et();
//	l.fsr_phi[i] = p4->Phi();
//	l.fsr_eta[i] = p4->Eta();
//      } else {
//	l.fsr_et[i] = -999.;
//	l.fsr_phi[i] = -999.;
//	l.fsr_eta[i] = -999.;
//      }
//    }
//    for (int i=0; i<l.el_std_n; i++) {
//      Int_t mc1=-1;
//      Int_t mc2=-1;
//      Int_t pho=-1;
//      
//      l.FindMCLeptons(i, mc1, mc2, pho, 11);
//      // et, eta, phi
//      
//      if (mc1 != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(mc1);
//	l.mc_et[i] = p4->Et();
//	l.mc_phi[i] = p4->Phi();
//	l.mc_eta[i] = p4->Eta();
//      } else {
//	l.mc_et[i] = -999.;
//	l.mc_phi[i] = -999.;
//	l.mc_eta[i] = -999.;
//      }
//      
//      if (pho != -1) {
//	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(pho);
//	l.fsr_et[i] = p4->Et();
//	l.fsr_phi[i] = p4->Phi();
//	l.fsr_eta[i] = p4->Eta();
//      } else {
//	l.fsr_et[i] = -999.;
//	l.fsr_phi[i] = -999.;
//	l.fsr_eta[i] = -999.;
//      }
//    }
//  }    
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::SelectEventsReduction(LoopAll& l, int jentry) {
  
  if (l.itype[l.current] == 0 && !checkEventHLT(l, hltSelection))
    if (!checkEventHLT(l, hltSelection))
      return false;
  
  // Two muons/electrons with pT > 20
  int goodMu = 0;
  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.mu_glo_p4->At(i);
    if (p4->Pt() > 20.)
      goodMu++;
  }
  
  int goodEl = 0;
  for(int i =0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.el_std_p4->At(i);
    if (p4->Et() > 20.)
      goodEl++;
  }

  return ((goodMu+goodEl) > 1);
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::SkimEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool SSLeptonAnalysis::SelectEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) {

  //l.higgs = new TLorentzVector(0,0,0,0);
  l.mu_glo_p4_corr = new TClonesArray("TLorentzVector", 100);
  
  l.Branch_mu_glo_p4_corr(outputTree);
  l.Branch_mu_glo_id_tight(outputTree);
  l.Branch_mu_glo_id_loose(outputTree);
  l.Branch_mc_et(outputTree);
  l.Branch_mc_eta(outputTree);
  l.Branch_mc_phi(outputTree);
  l.Branch_el_std_ch_ctf(outputTree);
  //l.Branch_fsr_et(outputTree);
  //l.Branch_fsr_eta(outputTree);
  //l.Branch_fsr_phi(outputTree);
  //l.Branch_higgs(outputTree);

  l.Branch_jet_algoPF1_genMatched(outputTree);
  l.Branch_jet_algoPF1_bgenMatched(outputTree);
  l.Branch_jet_algoPF1_cgenMatched(outputTree);
  l.Branch_jet_algoPF1_lgenMatched(outputTree);
  l.Branch_jet_algoPF1_vbfMatched(outputTree);
  l.Branch_jet_algoPF1_genPt(outputTree);
  l.Branch_jet_algoPF1_genDr(outputTree);

  l.Branch_shiftMET_pt(outputTree);
  l.Branch_shiftMET_phi(outputTree);
  l.Branch_smearMET_pt(outputTree);
  l.Branch_smearMET_phi(outputTree);
  l.Branch_shiftsmearMET_pt(outputTree);
  l.Branch_shiftsmearMET_phi(outputTree);
  l.Branch_shiftscaleMET_pt(outputTree);
  l.Branch_shiftscaleMET_phi(outputTree);
  l.Branch_shiftMET_eta(outputTree);
  l.Branch_shiftMET_e(outputTree);
  l.Branch_shiftscaleMET_eta(outputTree);
  l.Branch_shiftscaleMET_e(outputTree);
}

// ----------------------------------------------------------------------------------------------------
void SSLeptonAnalysis::ResetAnalysis()
{}

void SSLeptonAnalysis::Tree(LoopAll& l, Int_t pairs, Int_t* type, Float_t* sumPt,  Float_t* mass, Int_t* cat, Int_t* index1, Int_t* index2,
			    Float_t weight, Float_t pu_weight) {


  l.FillTree("run", l.run);
  l.FillTree("lumis", l.lumis);
  l.FillTree("event", (double)l.event);
  l.FillTree("itype", (int)l.itype[l.current]);
  l.FillTree("nvtx", (int)l.vtx_std_n);
  l.FillTree("rho", (float)l.rho_algo1);
  
  l.FillTree("pairs", pairs);
  l.FillTree("type",  type, pairs);
  l.FillTree("sumpt", sumPt, pairs);
  l.FillTree("mass",  mass, pairs);
  l.FillTree("cat",    cat, pairs);

  Float_t et1[100], eta1[100], phi1[100];
  Float_t et2[100], eta2[100], phi2[100];
  int ch1_1[100], ch2_1[100], ch3_1[100];
  int ch1_2[100], ch2_2[100], ch3_2[100];
  Float_t id1[100], id2[100], iso1[100], iso2[100];
  Float_t corrmet[100], corrmetPhi[100];
  for (unsigned int i=0; i<pairs; i++) {
    if (type[i] == 0) {
      id1[i] = 9999;
      id2[i] = 9999;
      iso1[i] = 9999;
      iso2[i] = 9999;
      ch_1[i] = l.mu_glo_charge[index1[i]];
      ch_2[i] = l.mu_glo_charge[index2[i]];
      ch1_1[i] = l.mu_glo_charge[index1[i]];
      ch2_1[i] = l.mu_glo_charge[index1[i]];
      ch3_1[i] = l.mu_glo_charge[index1[i]];
      ch1_2[i] = l.mu_glo_charge[index2[i]];
      ch2_2[i] = l.mu_glo_charge[index2[i]];
      ch3_2[i] = l.mu_glo_charge[index2[i]];
      et1[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index1[i]))->Et();
      eta1[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index1[i]))->Eta();
      phi1[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index1[i]))->Phi();
      et2[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i]))->Et();
      eta2[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i]))->Eta();
      phi2[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i]))->Phi();

      MetCorrections2012_Simple( l, *((TLorentzVector*) l.mu_glo_p4_corr->At(index1[i])), *((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i])));
      corrmet[i] = l.correctedpfMET;
      corrmetPhi[i] = l.correctedpfMET_phi;
      
    } else if (type[i] == 1) {
      id1[i]  = l.el_std_mva_trig[index1[i]];
      id2[i]  = l.el_std_mva_trig[index2[i]];
      iso1[i] = ElectronIsolation(l, index1[i], fabs(((TLorentzVector*)l.el_std_sc->At(index1[i]))->Eta()), ((TLorentzVector*) l.el_std_p4_corr->At(index1[i]))->Pt());
      iso2[i] = ElectronIsolation(l, index2[i], fabs(((TLorentzVector*)l.el_std_sc->At(index2[i]))->Eta()), ((TLorentzVector*) l.el_std_p4_corr->At(index2[i]))->Pt());
      ch1_1[i] = l.el_std_charge[index1[i]];
      ch2_1[i] = l.el_std_ch_gsf[index1[i]];
      ch3_1[i] = l.el_std_ch_scpix[index1[i]];
      ch1_2[i] = l.el_std_charge[index2[i]];
      ch2_2[i] = l.el_std_ch_gsf[index2[i]];
      ch3_2[i] = l.el_std_ch_scpix[index2[i]];
      ch_1[i] = l.el_std_charge[index1[i]];
      ch_2[i] = l.el_std_charge[index2[i]]; 
      et1[i] = ((TLorentzVector*)  l.el_std_p4_corr->At(index1[i]))->Et();
      eta1[i] = ((TLorentzVector*) l.el_std_p4_corr->At(index1[i]))->Eta();
      phi1[i] = ((TLorentzVector*) l.el_std_p4_corr->At(index1[i]))->Phi();
      et2[i] = ((TLorentzVector*)  l.el_std_p4_corr->At(index2[i]))->Et();
      eta2[i] = ((TLorentzVector*) l.el_std_p4_corr->At(index2[i]))->Eta();
      phi2[i] = ((TLorentzVector*) l.el_std_p4_corr->At(index2[i]))->Phi();
      
      MetCorrections2012_Simple( l, *((TLorentzVector*) l.el_std_p4_corr->At(index1[i])), *((TLorentzVector*) l.el_std_p4_corr->At(index2[i])));
      corrmet[i] = l.correctedpfMET;
      corrmetPhi[i] = l.correctedpfMET_phi;
    } else {
      id1[i]  = l.el_std_mva_trig[index1[i]];
      id2[i] = 9999;
      iso1[i] = ElectronIsolation(l, index1[i], fabs(((TLorentzVector*)l.el_std_sc->At(index1[i]))->Eta()), ((TLorentzVector*) l.el_std_p4_corr->At(index1[i]))->Pt());
      iso2[i] = 9999;
      ch_1[i] = l.el_std_charge[index1[i]];
      ch_2[i] = l.mu_glo_charge[index2[i]];
      ch1_1[i] = l.el_std_charge[index1[i]];
      ch2_1[i] = l.el_std_ch_gsf[index1[i]];
      ch3_1[i] = l.el_std_ch_scpix[index1[i]];
      ch1_2[i] = l.mu_glo_charge[index2[i]];
      ch2_2[i] = l.mu_glo_charge[index2[i]];
      ch3_2[i] = l.mu_glo_charge[index2[i]]; 
      et1[i] = ((TLorentzVector*)  l.el_std_p4_corr->At(index1[i]))->Et();
      eta1[i] = ((TLorentzVector*) l.el_std_p4_corr->At(index1[i]))->Eta();
      phi1[i] = ((TLorentzVector*) l.el_std_p4_corr->At(index1[i]))->Phi();
      et2[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i]))->Et();
      eta2[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i]))->Eta();
      phi2[i] = ((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i]))->Phi();

      MetCorrections2012_Simple( l, *((TLorentzVector*) l.el_std_p4_corr->At(index1[i])), *((TLorentzVector*) l.mu_glo_p4_corr->At(index2[i])));
      corrmet[i] = l.correctedpfMET;
      corrmetPhi[i] = l.correctedpfMET_phi;
    }
  }

  l.FillTree("et1", et1, pairs);
  l.FillTree("eta1", eta1, pairs);
  l.FillTree("phi1", phi1, pairs);
  l.FillTree("et2", et2, pairs);
  l.FillTree("eta2", eta2, pairs);
  l.FillTree("phi2", phi2, pairs);

  l.FillTree("ch1_1", ch1_1, pairs);
  l.FillTree("ch2_1", ch2_1, pairs);
  l.FillTree("ch3_1", ch3_1, pairs);
  l.FillTree("ch1_2", ch1_2, pairs);
  l.FillTree("ch2_2", ch2_2, pairs);
  l.FillTree("ch3_2", ch3_2, pairs);
  l.FillTree("ch_1", ch_1, pairs);
  l.FillTree("ch_2", ch_2, pairs);

  l.FillTree("id1",    id1, pairs);
  l.FillTree("iso1",   iso1, pairs);
  l.FillTree("id2",    id2, pairs);
  l.FillTree("iso2",   iso2, pairs);
  l.FillTree("weight", (float)weight);
  l.FillTree("pu_weight", (float)pu_weight);
  l.FillTree("met", (float)l.met_pfmet);
  l.FillTree("metPhi", (float)l.met_phi_pfmet);

  l.FillTree("corrmet", corrmet, pairs);
  l.FillTree("corrmetPhi", corrmetPhi, pairs);


  Bool_t jetid_flags[100];
  for(int ijet=0; ijet<l.jet_algoPF1_n; ++ijet ) {
    jetid_flags[ijet] = PileupJetIdentifier::passJetId((*l.jet_algoPF1_cutbased_wp_level_ext)[ijet][0], PileupJetIdentifier::kLoose);
  }

  int njets30 = l.SelectNHighestPtJets(l.el_std_p4_corr, l.mu_glo_p4_corr, jetid_flags, 30, 2.4, 0.6).size();
  l.FillTree("njets30", (int) njets30);

  std::vector<int> selectedJets = l.SelectNHighestPtJets(l.el_std_p4_corr, l.mu_glo_p4_corr, jetid_flags, 20, 2.4, 0.6);
  int njets = selectedJets.size();
  l.FillTree("njets20", (int) njets);

  Float_t jetet[100], jeteta[100], jetphi[100], jetm[100];
  for (int i=0; i<njets; i++) {
    TLorentzVector* j = (TLorentzVector*)l.jet_algoPF1_p4->At(selectedJets[i]);
    jetet[i] = j->Et();
    jeteta[i] = j->Eta();
    jetphi[i] = j->Phi();
    jetm[i] = j->M();
  }
  
  l.FillTree("jetet", jetet, njets);
  l.FillTree("jeteta", jeteta, njets);
  l.FillTree("jetphi", jetphi, njets);
  l.FillTree("jetM",  jetm, njets);


  Float_t btag[100];
  //Float_t btag2[100];
  if (njets > 0) {
    std::vector<float> jetEnergies;
    for (unsigned int i=0; i<selectedJets.size(); i++)
      jetEnergies.push_back(((TLorentzVector*)l.jet_algoPF1_p4->At(selectedJets[i]))->Et());

    for (unsigned int i1=0; i1<selectedJets.size()-1; ++i1) {
      for (unsigned int i2=i1+1; i2<selectedJets.size(); ++i2) {
	if (jetEnergies[i1] < jetEnergies[i2]) {
	  std::swap(selectedJets[i1], selectedJets[i2]);
	  std::swap(jetEnergies[i1], jetEnergies[i2]);
	}
      }
    }
    
    for (unsigned int i=0; i<njets; i++) 
      btag[i] = l.jet_algoPF1_csvBtag[selectedJets[i]];
  }

  l.FillTree("btag", btag, njets);
}

bool SSLeptonAnalysis::checkEventHLT(LoopAll& l, std::vector<std::string> paths) {
  
  bool result = false;

  std::vector<unsigned short> hltNumbers;
  
  for (unsigned int i=0; i<paths.size(); i++) {
    //std::cout << i << " " << paths[i] << std::endl;
    //std::cout << "_______________________" << std::endl;
    TRegexp e(TString(paths[i].c_str()), true);
    for (unsigned int j=0; j<l.hlt_path_names_HLT->size(); j++) {
      TString str1((*l.hlt_path_names_HLT)[j].c_str());
      //std::cout << (*l.hlt_path_names_HLT)[j] << std::endl;
      if (str1.Contains(e)) {
        //std::cout << (*l.hlt_path_names_HLT)[j] << std::endl;
        hltNumbers.push_back(j);
      }
    }
  }
  
  //system ("sleep 100000000");
  for (int j=0; j< hltNumbers.size(); j++) {
    for (int i=0; i<(*l.hlt_bit).size(); i++) {
      if (hltNumbers[j] == (*l.hlt_bit)[i]) {
        result = true;
        break;
      }
    }
  }

  return result;
}

void SSLeptonAnalysis::MetCorrections2012(LoopAll& l) {
  //shift met (reduction step)
  //in mc smearing should be applied first and then shifting
  //both these are performed at analysis step
  TLorentzVector unpfMET;
  unpfMET.SetPxPyPzE (l.met_pfmet*cos(l.met_phi_pfmet),
		      l.met_pfmet*sin(l.met_phi_pfmet),
		      0,
		      sqrt(l.met_pfmet*cos(l.met_phi_pfmet) * l.met_pfmet*cos(l.met_phi_pfmet)
			   + l.met_pfmet*sin(l.met_phi_pfmet) * l.met_pfmet*sin(l.met_phi_pfmet)));
  
  bool isMC = l.itype[l.current]!=0;
  
  TLorentzVector shiftMET_corr = l.shiftMet(&unpfMET,isMC);
  l.shiftMET_pt = shiftMET_corr.Pt();
  l.shiftMET_phi = shiftMET_corr.Phi();
  l.shiftMET_eta = shiftMET_corr.Eta();
  l.shiftMET_e = shiftMET_corr.Energy();
}


void SSLeptonAnalysis::MetCorrections2012_Simple(LoopAll& l,TLorentzVector lead_p4 ,TLorentzVector sublead_p4)
{
    // mc: scaling and shifting, data: scaling (analysis step)
    TLorentzVector unpfMET;
    unpfMET.SetPxPyPzE (l.met_pfmet*cos(l.met_phi_pfmet),
           l.met_pfmet*sin(l.met_phi_pfmet),
           0,
           sqrt(l.met_pfmet*cos(l.met_phi_pfmet) * l.met_pfmet*cos(l.met_phi_pfmet)
           + l.met_pfmet*sin(l.met_phi_pfmet) * l.met_pfmet*sin(l.met_phi_pfmet)));

     bool isMC = l.itype[l.current]!=0;

     //take shifted met for data
     TLorentzVector shiftedMET;
     double shiftedMETpt = l.shiftMET_pt;
     double shiftedMETe = l.shiftMET_e;
     double shiftedMETeta = l.shiftMET_eta;
     double shiftedMETphi = l.shiftMET_phi;

     shiftedMET.SetPtEtaPhiE(shiftedMETpt,shiftedMETeta,shiftedMETphi,shiftedMETe);
     if (isMC) {
       //smear raw met for mc
       TLorentzVector smearMET_corr = l.correctMet_Simple( lead_p4, sublead_p4 , &unpfMET, true, false);
       l.smearMET_pt = smearMET_corr.Pt();
       l.smearMET_phi = smearMET_corr.Phi();
       //shift smeared met for mc
       TLorentzVector shiftsmearMET_corr = l.shiftMet(&smearMET_corr,isMC);
       l.shiftsmearMET_pt = shiftsmearMET_corr.Pt();
       l.shiftsmearMET_phi = shiftsmearMET_corr.Phi();
       l.correctedpfMET = l.shiftsmearMET_pt;
       l.correctedpfMET_phi = l.shiftsmearMET_phi;
     } else {
       //scale shifted met for data
       TLorentzVector shiftscaleMET_corr = l.correctMet_Simple( lead_p4, sublead_p4 , &shiftedMET, false , true);
       l.shiftscaleMET_pt = shiftscaleMET_corr.Pt();
       l.shiftscaleMET_phi = shiftscaleMET_corr.Phi();
       l.correctedpfMET = l.shiftscaleMET_pt;
       l.correctedpfMET_phi = l.shiftscaleMET_phi;
     }
}

float SSLeptonAnalysis::ElectronIsolation(LoopAll& l, int el_ind, float eta, float pT) {

  if(eta>2.5 || (eta>1.442 && eta<1.566)) 
    return 9999.;

  double Aeff=0.;
  
  if(eta<1.0)                   Aeff=0.135;
  if(eta>=1.0 && eta<1.479) Aeff=0.168;
  if(eta>=1.479 && eta<2.0) Aeff=0.068;
  if(eta>=2.0 && eta<2.2)   Aeff=0.116;
  if(eta>=2.2 && eta<2.3)   Aeff=0.162;
  if(eta>=2.3 && eta<2.4)   Aeff=0.241;
  if(eta>=2.4)                  Aeff=0.23;
  float thisiso=l.el_std_pfiso_charged[el_ind]+std::max(l.el_std_pfiso_neutral[el_ind]+l.el_std_pfiso_photon[el_ind]-l.rho_algo1*Aeff,0.);
    
  return thisiso/pT;
}

bool SSLeptonAnalysis::ElectronMVACuts(LoopAll& l, int el_ind) {

  bool pass=false;
  if(el_ind < 0 || el_ind >= l.el_std_n) 
    return pass;

  if(l.el_std_mva_trig[el_ind] < eleIDCut) 
    return pass;

  TLorentzVector* thisel = (TLorentzVector*) l.el_std_p4_corr->At(el_ind);
  TLorentzVector* thissc = (TLorentzVector*) l.el_std_sc->At(el_ind);
  float thiseta = fabs(thissc->Eta());
  float thispt = thisel->Pt();
  
  float relIso = ElectronIsolation(l, el_ind, thiseta, thispt);
  if (relIso > eleIsoCut)
    return pass;
  
  if(l.el_std_hp_expin[el_ind]>1) 
    return pass;

  //if(l.el_std_conv[el_ind]==0)    return pass;

  pass=true;
  return pass;
}
