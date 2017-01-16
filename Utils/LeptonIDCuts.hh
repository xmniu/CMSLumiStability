#ifndef EWKANA_UTILS_LEPTONIDCUTS_HH
#define EWKANA_UTILS_LEPTONIDCUTS_HH

#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
#include "BaconAna/Utils/interface/TTrigger.hh"
#include "BaconAna/DataFormats/interface/TElectron.hh"
#include "BaconAna/DataFormats/interface/TMuon.hh"
#include <TMath.h>
#include <cassert>
		   
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho=0, const Int_t isoType=1);
Bool_t passAntiMuonID(const baconhep::TMuon *muon, const Double_t rho=0);
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho=0);

Bool_t passEleID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho=0);
Bool_t passEleTightID(const baconhep::TElectron *electron, const Double_t rho);
Bool_t passEleLooseID(const baconhep::TElectron *electron, const Double_t rho=0);
Bool_t passAntiEleID(const baconhep::TElectron *electron, const Double_t rho=0);

Bool_t isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Int_t triType);
Bool_t isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Int_t triType);

Bool_t isEleTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Bool_t isData, Int_t triType);
Bool_t isEleTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Bool_t isData, Int_t triType);

Double_t getEffAreaEl(const Double_t eta);
Double_t getEffAreaMu(const Double_t eta);

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const baconhep::TMuon *muon, const Double_t rho, const Int_t isoType)
{
  if(muon->nTkLayers  < 6)            return kFALSE;
  if(muon->nPixHits   < 1)            return kFALSE;
  if(fabs(muon->d0)   > 0.2)          return kFALSE;
  if(fabs(muon->dz)   > 0.5)          return kFALSE;
  if(muon->muNchi2    > 10)           return kFALSE;
  if(muon->nMatchStn  < 2)            return kFALSE;
  if(muon->nValidHits < 1)            return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kGlobal)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
 
  if(isoType==0){
  //No isolation requirement in selection
  }
  else if(isoType==1){
  //PF-based combined relative, deltaBeta-corrected
  //Cone: 0.4
  //Tight cut (eff 95%): 0.15 
    //cout<<"pfCombIso 15%"<<endl;
    Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
    if(iso > 0.15*(muon->pt)) return kFALSE;
  }
  else if(isoType==2){
  //PF-based combined relative, deltaBeta-corrected
  //Cone 0.4
  //Loose cut (eff 98%): 0.25
    Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
    if(iso > 0.25*(muon->pt)) return kFALSE;
  }
  else if(isoType==3){
  //Tracker-based relative
  //Cone 0.3
  //Tight cut (eff 95%): 0.05
    Double_t iso = muon->trkIso;
    if(iso > 0.05*(muon->pt)) return kFALSE;
  }
  else if(isoType==4){
  //Tracker-based relative
  //Cone 0.3
  //Loose cut (eff 98%): 0.10
    Double_t iso = muon->trkIso;
    if(iso > 0.10*(muon->pt)) return kFALSE;
  }
  else{
    cout<<"Wrong input for isoType..."<<endl;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passAntiMuonID(const baconhep::TMuon *muon, const Double_t rho)
{
  if(muon->nTkLayers  < 6)            return kFALSE;
  if(muon->nPixHits   < 1)            return kFALSE;
  if(fabs(muon->d0)   > 0.2)          return kFALSE;
  if(fabs(muon->dz)   > 0.5)          return kFALSE;
  if(muon->muNchi2    > 10)           return kFALSE;
  if(muon->nMatchStn  < 2)            return kFALSE;
  if(muon->nValidHits < 1)            return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kGlobal)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;
  
  Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
  if(iso < 0.3*(muon->pt)) return kFALSE;

  return kTRUE;
}


//--------------------------------------------------------------------------------------------------
Bool_t passMuonLooseID(const baconhep::TMuon *muon, const Double_t rho)
{
  if(!(muon->typeBits & baconhep::EMuType::kGlobal) && !(muon->typeBits & baconhep::EMuType::kTracker)) return kFALSE;
  if(!(muon->typeBits & baconhep::EMuType::kPFMuon)) return kFALSE;

  Double_t iso = muon->chHadIso + TMath::Max(muon->neuHadIso + muon->gammaIso - 0.5*(muon->puIso),Double_t(0));
  if(iso > 0.15*(muon->pt)) return kFALSE;
  
  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const baconhep::TElectron *electron, const TLorentzVector tag, const Double_t rho)
{ // Medium Electron ID for PU20 bx25

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;

  if(electron->isConv)            return kFALSE;

  Double_t ea = getEffAreaEl(tag.Eta());
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    if(iso >= 0.0766*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 2)                                    return kFALSE;
    if(electron->sieie >= 0.0101)                                     return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.0336)                              return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.0103)                              return kFALSE;
    if(electron->hovere >= 0.0876)                                    return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0174*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0118)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.373)                                  return kFALSE;
  } else {
    if(iso >= 0.0678*(tag.Pt()))                                      return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.0283)                          return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.114)                           return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.00733)                         return kFALSE;
    if(electron->hovere           >= 0.0678)                          return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.0898*(electron->ecalEnergy))   return kFALSE;
    if(fabs(electron->d0) >= 0.0739)                                  return kFALSE;
    if(fabs(electron->dz) >= 0.602)                                  return kFALSE;
  }

  return kTRUE;

}

//--------------------------------------------------------------------------------------------------
Bool_t passEleTightID(const baconhep::TElectron *electron, const Double_t rho)
{ // Phys14 Tight Electron ID for PU20 bx25

  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;

  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;

  if(!(electron->typeBits & baconhep::EEleType::kEcalDriven)) return kFALSE;

  // conversion rejection
  if(electron->isConv) return kFALSE;

  Double_t ea = getEffAreaEl(electron->scEta);
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  // barrel/endcap dependent requirements
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    if(iso >= 0.069537*(electron->pt))                                return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie >= 0.009947)                                   return kFALSE;
    if(fabs(electron->dPhiIn) >= 0.028092)                            return kFALSE;
    if(fabs(electron->dEtaIn) >= 0.006046)                            return kFALSE;
    if(electron->hovere >= 0.045772)                                  return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.020118*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) >= 0.008790)                                return kFALSE;
    if(fabs(electron->dz) >= 0.021226)                                return kFALSE;
  } else {
    // endcap
    if(iso >= 0.078265*(electron->pt))                                return kFALSE;
    if(electron->nMissingHits > 1)                                    return kFALSE;
    if(electron->sieie            >= 0.028237)                        return kFALSE;
    if(fabs(electron->dPhiIn)     >= 0.030159)                        return kFALSE;
    if(fabs(electron->dEtaIn)     >= 0.007057)                        return kFALSE;
    if(electron->hovere           >= 0.067778)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) >= 0.098919*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) >= 0.027984)                                return kFALSE;
    if(fabs(electron->dz) >= 0.133431)                                return kFALSE;
  }

  return kTRUE;

}

//--------------------------------------------------------------------------------------------------
Bool_t passAntiEleID(const baconhep::TElectron *electron, const Double_t rho)
{ // CSA14 Medium working point
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & baconhep::EEleType::kEcalDriven)) return kFALSE;
  
  // conversion rejection
  if(electron->isConv) return kFALSE;
     
  Double_t ea = getEffAreaEl(electron->scEta);
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);

  // barrel/endcap dependent requirements      
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    if(iso > 0.13*(electron->pt))                                 return kFALSE;
    if(electron->nMissingHits > 2)                                return kFALSE;
    if(electron->sieie  	  > 0.011)                        return kFALSE;
    if(fabs(electron->dPhiIn)     > 0.033)                        return kFALSE;
    if(fabs(electron->dEtaIn)     > 0.0011)                       return kFALSE;
    if(electron->hovere 	  > 0.091)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.034*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) > 0.012)                                return kFALSE;
    if(fabs(electron->dz) > 0.39)                                 return kFALSE;
  } else {
    // endcap
    if(iso > 0.13*(electron->pt))                                 return kFALSE;
    if(electron->nMissingHits > 1)                                return kFALSE;
    if(electron->sieie  	  > 0.031)                        return kFALSE;
    if(fabs(electron->dPhiIn)     > 0.047)                        return kFALSE;
    if(fabs(electron->dEtaIn)     > 0.024)                        return kFALSE;
    if(electron->hovere 	  > 0.099)                        return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.086*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0) > 0.016)                                return kFALSE;
    if(fabs(electron->dz) > 0.78)                                 return kFALSE;
  }
  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t passEleLooseID(const baconhep::TElectron *electron, const Double_t rho)
{ // Phys14 Veto working point
  const Double_t ECAL_GAP_LOW  = 1.4442;
  const Double_t ECAL_GAP_HIGH = 1.566;
  
  if((fabs(electron->scEta)>ECAL_GAP_LOW) && (fabs(electron->scEta)<ECAL_GAP_HIGH)) return kFALSE;
  
  if(!(electron->typeBits & baconhep::EEleType::kEcalDriven)) return kFALSE;

  // conversion rejection
  if(electron->isConv) return kFALSE;
       
  Double_t ea = getEffAreaEl(electron->scEta);
  Double_t iso = electron->chHadIso + TMath::Max(electron->neuHadIso + electron->gammaIso - rho*ea, 0.);     

  // barrel/endcap dependent requirements
  if(fabs(electron->scEta)<=ECAL_GAP_LOW) {
    // barrel
    if(iso > 0.158721*(electron->pt))                                return kFALSE;
    if(electron->nMissingHits > 2)                                   return kFALSE;
    if(electron->sieie        > 0.011586)                            return kFALSE;
    if(fabs(electron->dPhiIn) > 0.230374)                            return kFALSE;
    if(fabs(electron->dEtaIn) > 0.013625)                            return kFALSE;
    if(electron->hovere       > 0.181130)                            return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.295751*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0)     > 0.094095)                            return kFALSE;
    if(fabs(electron->dz)     > 0.713070)                            return kFALSE;
  } else {
    // endcap
    if(iso > 0.177032*(electron->pt))                                return kFALSE;
    if(electron->nMissingHits > 3)                                   return kFALSE;
    if(electron->sieie        > 0.031849)                            return kFALSE;
    if(fabs(electron->dPhiIn) > 0.255450)                            return kFALSE;
    if(fabs(electron->dEtaIn) > 0.011932)                            return kFALSE;
    if(electron->hovere       > 0.223870)                            return kFALSE;
    if(fabs(1.0-electron->eoverp) > 0.155501*(electron->ecalEnergy)) return kFALSE;
    if(fabs(electron->d0)     > 0.342293)                            return kFALSE;
    if(fabs(electron->dz)     > 0.953461)                            return kFALSE;
  }

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t isMuonTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Int_t triType) {
  if (triType==0) return triggerMenu.pass("HLT_IsoMu20_v*",hltBits);
  else if (triType==1) return triggerMenu.pass("HLT_Mu20_v*",hltBits);
  else if (triType==2) return triggerMenu.pass("HLT_Mu50_v*",hltBits);
  else if (triType==3) return (triggerMenu.pass("HLT_IsoMu20_v*",hltBits) || triggerMenu.pass("HLT_IsoTkMu20_eta2p1_v*",hltBits)); 
  else if (triType==4) return triggerMenu.pass("HLT_IsoMu22_v*",hltBits);
  else if (triType==5) return triggerMenu.pass("HLT_IsoMu24_v*",hltBits);
  else if (triType==6) return triggerMenu.pass("HLT_IsoMu27_v*",hltBits);
  else return kFALSE; 
}

Bool_t isMuonTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Int_t triType) {
  if (triType==0) {
    if (isL1) return triggerMenu.passObj("HLT_IsoMu20_v*","hltL1sL1SingleMu18",hltMatchBits);
//    else return triggerMenu.passObj("HLT_IsoMu20_v*","hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09",hltMatchBits);
    else return triggerMenu.passObj("HLT_IsoMu20_v*","hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09",hltMatchBits);
  }
  else if (triType==1) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Mu20_v*","hltL3fL1sMu18L1f0L2f10QL3Filtered20Q",hltMatchBits);
  }
  else if (triType==2) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Mu50_v*","hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q",hltMatchBits);
  }
  else if (triType==3) {
    if (isL1) return kFALSE;
    else return (triggerMenu.passObj("HLT_IsoMu20_v*","hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09",hltMatchBits) || triggerMenu.passObj("HLT_IsoTkMu20_eta2p1_v*","hltL3crIsoL1sMu18Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09",hltMatchBits));
  }
  else if (triType==4) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_IsoMu22_v*","hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09",hltMatchBits);
  }
  else if (triType==5) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_IsoMu24_v*","hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09",hltMatchBits);
  }
  else return kFALSE; 
}

Bool_t isEleTrigger(baconhep::TTrigger triggerMenu, TriggerBits hltBits, Bool_t isData, Int_t triType) {
    if (triType==0) return triggerMenu.pass("HLT_Ele23_WPLoose_Gsf_v*",hltBits);
    else if (triType==1) return triggerMenu.pass("HLT_Ele25_eta2p1_WPTight_Gsf_v*",hltBits);
    else if (triType==2) return triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",hltBits);
    else if (triType==3) return triggerMenu.pass("HLT_Ele27_eta2p1_WPLoose_Gsf_*",hltBits);
    else if (triType==4) return triggerMenu.pass("HLT_Ele22_WP75_Gsf_v*",hltBits);
    else if (triType==5) return triggerMenu.pass("HLT_Ele23_WP75_Gsf_v*",hltBits);
    else if (triType==6) return triggerMenu.pass("HLT_Ele22_eta2p1_WP75_Gsf_v*",hltBits);
    else return kFALSE;
}

Bool_t isEleTriggerObj(baconhep::TTrigger triggerMenu, TriggerObjects hltMatchBits, Bool_t isL1, Bool_t isData, Int_t triType) {
  if (triType==0) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele23_WPLoose_Gsf_v*","hltEle23WPLooseGsfTrackIsoFilter",hltMatchBits);
  }
  else if (triType==1) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele25_eta2p1_WPTight_Gsf_v* ","hltEle25erWPTightGsfTrackIsoFilter",hltMatchBits);
  }
  else if (triType==2) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele27_WPTight_Gsf_v*","hltEle27WPTightGsfTrackIsoFilter",hltMatchBits);
  }
  else if (triType==3) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele27_eta2p1_WPLoose_Gsf_*","hltEle27erWPLooseGsfTrackIsoFilter",hltMatchBits);
  }
  else if (triType==4) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele22_WP75_Gsf_v*","hltEle22WP75GsfTrackIsoFilter",hltMatchBits);
  }
  else if (triType==5) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele23_WP75_Gsf_v*","hltEle23WP75GsfTrackIsoFilter",hltMatchBits);
  }
  else if (triType==6) {
    if (isL1) return kFALSE;
    else return triggerMenu.passObj("HLT_Ele22_eta2p1_WP75_Gsf_v*","hltEle22WP75GsfTrackIsoFilter",hltMatchBits);
  }
    else return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Double_t getEffAreaEl(const Double_t eta) {
  if      (fabs(eta) < 1.0) return 0.1752;
  else if (fabs(eta) < 1.479)  return 0.1862;
  else if (fabs(eta) < 2.0)  return 0.1411;
  else if (fabs(eta) < 2.2)  return 0.1534;
  else if (fabs(eta) < 2.3)  return 0.1903;
  else if (fabs(eta) < 2.4)  return 0.2243;
  else if (fabs(eta) < 2.5)  return 0.2687;
  else return -1; // This should never happen, a cut on |eta|<2.5 is applied before this function is used.
}

Double_t getEffAreaMu(const Double_t eta) {
  // not used, so i didn't update? 
  // well, ok this probably should be used....
  if      (fabs(eta) < 0.8) return 0.0913;
  else if (fabs(eta) < 1.3)  return 0.0765;
  else if (fabs(eta) < 2.0)  return 0.0546;
  else if (fabs(eta) < 2.2)  return 0.0728;
  else if (fabs(eta) < 2.5)  return 0.1177;
  else return -1; // This should never happen, a cut on |eta|<2.4 is applied in selection.
}

#endif

