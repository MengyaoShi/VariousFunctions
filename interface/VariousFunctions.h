#ifndef VariousFunctions_interface_VariousFunctions_h
#define VariousFunctions_interface_VariousFunctions_h

#include <vector>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include <string>
#include "TH2F.h"
#include "THStack.h"
#include "TF2.h"
class VariousFunctions { 

  public:
    static reco::GenParticleRef findDaughterInDaughters(const reco::GenParticleRef& , const double, const bool);
    static bool findIfInDaughters(const reco::GenParticleRef& , const double, const bool);
    static int tauDecayMode(const reco::GenParticleRef &);
    static reco::LeafCandidate::LorentzVector sumTauP4(const reco::GenParticleRef&, const int, const bool);
    static double getDiTauDR(const reco::GenParticleRef&, const reco::GenParticleRef&, const bool);
    static double getMuTauDR(const reco::GenParticleRef&, const reco::GenParticleRef&, const bool);
    static double getDiThingDR_1(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static void formatAndDrawCanvasAndHist2D(TCanvas&, TH2F*, 
                                           const Int_t, const Int_t, const Int_t, 
                                           const Color_t, const Size_t, const Style_t,
                                           const char*, const Float_t, const Float_t, const Float_t, 
				           const char*, const Float_t, const Float_t, const Float_t,
  					   const char*, const Float_t, const Float_t, const Float_t);
    static void DrawHstack(TCanvas&, TH1F*, TH1F*);
    static double getHigherPt(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static double getLowerPt(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static reco::GenParticleRef getHigherPtObj(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static reco::GenParticleRef getLowerPtObj(const reco::GenParticleRef&, const reco::GenParticleRef&);
    static double getInvMass(const reco::GenParticleRef&, const reco::GenParticleRef&);
};
#endif
