
#ifndef  HiggsTools_ANALYZERS_H
#define  HiggsTools_ANALYZERS_H

#include <cmath>
#include <vector>
#include <iostream>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"
#include "MCParticle.h"

namespace HiggsTools{
	///build the resonance from 2 particles from an arbitrary list of input ReconstructedPartilces. Keep the closest to the mass given as input
	struct resonanceZBuilder {
		float m_resonance_mass;
		resonanceZBuilder(float arg_resonance_mass);
		ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) ;
	};

  // temporary duplication. When arg_use_MC_Kinematics = true, the true kinematics will be used instead of the track momenta
  struct resonanceZBuilder2 {
          float m_resonance_mass;
          bool m_use_MC_Kinematics;
          resonanceZBuilder2(float arg_resonance_mass, bool arg_use_MC_Kinematics);
          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
                          ROOT::VecOps::RVec<int> recind,
                          ROOT::VecOps::RVec<int> mcind,
                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                          ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) ;
  };

  struct resonanceZBuilderHiggsPairs {
          float m_resonance_mass;
          bool m_use_MC_Kinematics;
          resonanceZBuilderHiggsPairs(float arg_resonance_mass, bool arg_use_MC_Kinematics);
          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
                          ROOT::VecOps::RVec<int> recind,
                          ROOT::VecOps::RVec<int> mcind,
                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                          ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                          ROOT::VecOps::RVec<int> parents,
                          ROOT::VecOps::RVec<int> daugthers);
  };
  // build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV^M
  // technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair^M
  struct resonanceBuilder_mass_recoil {
          float m_resonance_mass;
          float m_recoil_mass;
          float chi2_recoil_frac;
          float ecm;
          bool m_use_MC_Kinematics;
          resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
                          ROOT::VecOps::RVec<int> recind,
                          ROOT::VecOps::RVec<int> mcind,
                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                          ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                          ROOT::VecOps::RVec<int> parents,
                          ROOT::VecOps::RVec<int> daugthers) ;
  };
  
  struct resonanceBuilder_mass_recoil2 {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass_recoil2(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
                    ROOT::VecOps::RVec<int> recind,
                    ROOT::VecOps::RVec<int> mcind,
                    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                    ROOT::VecOps::RVec<int> parents,
                    ROOT::VecOps::RVec<int> daugthers) ;
}; 

  /// muon scale shifts
  struct momentum_scale {
    momentum_scale(float arg_scaleunc);
    float scaleunc = 1.;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };
  /// to be added to your ReconstructedParticle.h :
  
  /// select ReconstructedParticles with a given type 
  struct sel_type {
    sel_type( int arg_pdg, bool arg_chargeconjugate);
    int m_pdg = 13;
    bool m_chargeconjugate = true;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  struct sel_isol {
    sel_isol(float arg_isocut);
    float m_isocut = 9999.;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> particles, ROOT::VecOps::RVec<float> var );
  };
  
  /// boost along x of the ReconstructedParticles 
  struct BoostAngle {
    BoostAngle( float arg_angle );
    float m_angle = -0.015 ;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in );
  };

  struct gen_sel_pdgIDInt {
    gen_sel_pdgIDInt(int arg_pdg, bool arg_chargeconjugate);
    int m_pdg = 13;
    bool m_chargeconjugate = true;
    ROOT::VecOps::RVec<int>  operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in);
  };

  struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
            
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    ROOT::VecOps::RVec<double>  operator() ( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                                             ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rp) ;
  };

  // Same as coneIsolation, but the cone distance is built from the polar angle (theta)
  // instead of pseudorapidity (eta): dR = sqrt(dTheta^2 + dPhi^2). theta is the natural
  // angular variable at an e+e- collider, so this makes the metric consistent (referee
  // request). The eta-based coneIsolation remains for backward compatibility.
  struct coneIsolationTheta {

    coneIsolationTheta(float arg_dr_min, float arg_dr_max);

    double deltaR(double theta1, double phi1, double theta2, double phi2) { return TMath::Sqrt(TMath::Power(theta1-theta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    ROOT::VecOps::RVec<double>  operator() ( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                                             ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rp) ;
  };

  // calculate the number of foward leptons
  struct polarAngleCategorization {
    polarAngleCategorization(float arg_thetaMin, float arg_thetaMax);
    float thetaMin = 0;
    float thetaMax = 5;
    int operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  // perturb the momentum scale with a given constant
  struct lepton_momentum_scale {
    lepton_momentum_scale(float arg_scaleunc);
    float scaleunc = 1.;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  struct BoostFrame_gen {
    BoostFrame_gen( float arg_low_energy );
    float low_energy = 30;
    ROOT::VecOps::RVec<edm4hep::MCParticleData> operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in);
  };

  struct BoostFrame {
    BoostFrame( float arg_low_energy );
    float low_energy = 30;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  };

  struct resonanceBuilder_mass_recoil_boosted {
        float m_resonance_mass;
        float m_recoil_mass;
        float chi2_recoil_frac;
        float ecm;
        bool m_use_MC_Kinematics;
        resonanceBuilder_mass_recoil_boosted(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics);
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
                        ROOT::VecOps::RVec<int> recind,
                        ROOT::VecOps::RVec<int> mcind,
                        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                        ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                        ROOT::VecOps::RVec<int> parents,
                        ROOT::VecOps::RVec<int> daugthers) ;
  }; 

  struct recoilBuilder_boosted {
    recoilBuilder_boosted(float arg_sqrts);
    float m_sqrts = 240.0;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) ;
  };

  /// return the costheta of the input ReconstructedParticles
	ROOT::VecOps::RVec<float> get_cosTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

	/// return the costheta of the input missing momentum 
  ROOT::VecOps::RVec<float> get_cosTheta_miss(ROOT::VecOps::RVec<Float_t>Px, ROOT::VecOps::RVec<Float_t>Py, ROOT::VecOps::RVec<Float_t>Pz, ROOT::VecOps::RVec<Float_t>E);
  ///return muon_quality_check result (at least one muon plus and one muon minus)
	ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  muon_quality_check(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sort_greater(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sort_greater_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  //ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  get_subleading(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float Reweighting_wzp_kkmc(float pT, float m);

  ///get acolinearity
  ROOT::VecOps::RVec<float> acolinearity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  
  ///get acoplanarity
  ROOT::VecOps::RVec<float> acoplanarity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
 
  /// isolation of a RecoPrticle wrt the others
  ROOT::VecOps::RVec<float> Isolation( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> particles, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in );
  
  /// cf Delphes Merger module - used for glbal sum
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> Merger( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) ;
  
  bool from_Higgsdecay(int i, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind);
  bool from_Higgsdecay(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents, ROOT::VecOps::RVec<int> daugther);
  // get the gen p from reco
  std::vector<float> gen_p_from_reco(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);

  ROOT::VecOps::RVec<edm4hep::MCParticleData> get_photons(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);
  ROOT::VecOps::RVec<float> leptonResolution_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons, ROOT::VecOps::RVec<int> recind,
                                               ROOT::VecOps::RVec<int> mcind,
                                               ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                                               ROOT::VecOps::RVec<edm4hep::MCParticleData> mc);
  // return list of pdg from decay of a list of mother particle
  std::vector<int> gen_decay_list(ROOT::VecOps::RVec<int> mcin, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind);
  
  //Higgstrahlungness
  float Higgsstrahlungness(float mll, float mrecoil);
  float deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float deltaRPrime(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float deltaPhi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float deltaTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float deltaEta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float min_deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float min_deltaRPrime(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float min_deltaPhi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float min_deltaTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float min_deltaEta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float max_deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float max_deltaRPrime(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float max_deltaPhi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float max_deltaTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);
  float max_deltaEta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> whizard_zh_select_prompt_leptons(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, 
                                                                                          ROOT::VecOps::RVec<int> recind, 
                                                                                          ROOT::VecOps::RVec<int> mcind, 
                                                                                          ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, 
                                                                                          ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, 
                                                                                          ROOT::VecOps::RVec<int> parents, 
                                                                                          ROOT::VecOps::RVec<int> daugther);

  bool whizard_zh_from_prompt(int i, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind);

  ROOT::VecOps::RVec<edm4hep::MCParticleData> get_gen_pdg(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, int pdgId, bool abs);

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> missing(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, float ecm, float p_cutoff = 0.0);
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> visible(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, float p_cutoff = 0.0);

  float ZHChi2(float mZ, float mH, float chi2_H_frac = 0.5);

  // ------------------------------------------------------------------
  // ZH hadronic (Z->qq) helpers, ported from jeyserma/FCCPhysics
  // (analyses/h_zh/utils_hadronic.h, utils.h and functions/functions.h)
  // ------------------------------------------------------------------

  /// convert clustered jets (px, py, pz, E, m vectors) to ReconstructedParticleData
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> jets2rp(ROOT::VecOps::RVec<float> px, ROOT::VecOps::RVec<float> py, ROOT::VecOps::RVec<float> pz, ROOT::VecOps::RVec<float> e, ROOT::VecOps::RVec<float> m);

  /// jet quality selection (p > 5 GeV), can reduce the jet multiplicity
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> select_jets(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, std::vector<std::vector<int>> constituents, int njets_sel, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco);

  /// pick the best clustering hypothesis (inclusive, N=2, 4, 6) via a chi2 on m(qq), p(qq), recoil
  int best_clustering_idx(ROOT::VecOps::RVec<float> mz, ROOT::VecOps::RVec<float> pz, ROOT::VecOps::RVec<float> mrec, ROOT::VecOps::RVec<int> njets, ROOT::VecOps::RVec<int> njets_target, int ecm = 240);

  /// build the Z(qq) resonance from the jet pair minimizing a chi2 on mass, recoil and momentum
  struct resonanceBuilder_mass_recoil_hadronic {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    resonanceBuilder_mass_recoil_hadronic(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs);
  };

  /// pair 4 jets into two W candidates minimizing the (mW, mW) chi2
  ROOT::VecOps::RVec<TLorentzVector> pair_WW_N4(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// MC truth: inclusive WW event where both Ws decay to the same lepton flavor (e or mu, incl. via tau)
  bool is_ww_leptonic(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind);

  /// MC truth: H->ZZ event with both Zs decaying invisibly
  bool is_hzz_invisible(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> ind);

  /// filter reco particles based on a property within (min, max); abs = compare |prop|
  struct sel_range {
    sel_range(float arg_min, float arg_max, bool arg_abs = false);
    float m_min = 0.;
    float m_max = 1.;
    bool m_abs = false;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<float> prop);
  };

  /// 3D acolinearity acos(-v1.v2/|v1||v2|) between the two leading particles (FCCPhysics definition, scalar)
  float acolinearity_scalar(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// acoplanarity pi - |dphi| between the two leading particles (scalar)
  float acoplanarity_scalar(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in);

  /// missing-energy four-vector from the visible reco particles
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> missingEnergy(float ecm, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, float p_cutoff = 0.0);

  /// |cos(theta)| of the missing-energy vector
  float get_cosTheta_miss(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> met);
}
#endif
