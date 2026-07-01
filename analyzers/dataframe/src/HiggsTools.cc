#include "FCCAnalyses/HiggsTools.h"
#include "FCCAnalyses/ReconstructedParticle.h"
#include <algorithm>

using namespace HiggsTools;

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  HiggsTools::muon_quality_check(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
	ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  //at least one muon + and one muon - in each event
  int n_muon_plus = 0;
	int n_muon_minus = 0;
	int n = in.size();
	for (int i = 0; i < n; ++i) {
		if (in[i].charge == 1.0){
			++n_muon_plus;
		}
		else if (in[i].charge == -1.0){
			++n_muon_minus;
		}
	}
	if (n_muon_plus >= 1 && n_muon_minus >= 1){
		result = in;
	}
	return result;
}

ROOT::VecOps::RVec<float> HiggsTools::get_cosTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
   ROOT::VecOps::RVec<float> result;
	 for (auto & p: in) {
		 TLorentzVector tlv;
		 tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
		 result.push_back(tlv.CosTheta());
	 }
	 return result;
}

ROOT::VecOps::RVec<float> HiggsTools::get_cosTheta_miss(ROOT::VecOps::RVec<Float_t>Px, ROOT::VecOps::RVec<Float_t>Py, ROOT::VecOps::RVec<Float_t>Pz, ROOT::VecOps::RVec<Float_t>E) {
  ROOT::VecOps::RVec<float> result;
  for (int i =0; i < Px.size(); ++i) {
		TLorentzVector tlv;
		tlv.SetPxPyPzE(Px.at(i), Py.at(i), Pz.at(i), E.at(i));
    result.push_back(tlv.CosTheta());
  }
  return result;
} 


HiggsTools::resonanceZBuilder::resonanceZBuilder(float arg_resonance_mass) {m_resonance_mass = arg_resonance_mass;}
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::resonanceZBuilder::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs) {

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  int n = legs.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      //set initial charge == 0
      reso.charge = 0;
      TLorentzVector reso_lv;
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
                                //prevent +2 and -2 charged Z 
            if (reso.charge == legs[i].charge) continue;
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;
            leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&] (edm4hep::ReconstructedParticleData i ,edm4hep::ReconstructedParticleData j) { return (abs( m_resonance_mass -i.mass)<abs(m_resonance_mass-j.mass)); };
                std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator first = result.begin();
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator last = result.begin() + 1;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}



HiggsTools::resonanceZBuilder2::resonanceZBuilder2(float arg_resonance_mass, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_use_MC_Kinematics = arg_use_MC_Kinematics;}
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::resonanceZBuilder2::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
				ROOT::VecOps::RVec<int> recind ,
				ROOT::VecOps::RVec<int> mcind ,
				ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco ,
				ROOT::VecOps::RVec<edm4hep::MCParticleData> mc )   {

  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  int n = legs.size();
  if (n >1) {
    ROOT::VecOps::RVec<bool> v(n);
    std::fill(v.end() - 2, v.end(), true);
    do {
      edm4hep::ReconstructedParticleData reso;
      //set initial charge == 0
      reso.charge = 0;
      TLorentzVector reso_lv; 
      for (int i = 0; i < n; ++i) {
          if (v[i]) {
    				//prevent +2 and -2 charged Z 
            if (reso.charge == legs[i].charge) continue;
            reso.charge += legs[i].charge;
            TLorentzVector leg_lv;

		// Ideal detector resolution: use the kinematics of the MC particle instead
		if ( m_use_MC_Kinematics) {

		     // ugly: particles_begin is not filled in RecoParticle.
		     // hence: either need to keep trace of the index of the legs into the RecoParticle collection,
		     // or (as done below) use the track index to map the leg to the MC particle :-(

		     int track_index = legs[i].tracks_begin ;   // index in the Track array
		     int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index( track_index, recind, mcind, reco );
		     if ( mc_index >= 0 && mc_index < mc.size() ) {
			 int pdgID = mc.at( mc_index).PDG;
		         leg_lv.SetXYZM(mc.at(mc_index ).momentum.x, mc.at(mc_index ).momentum.y, mc.at(mc_index ).momentum.z, mc.at(mc_index ).mass );
		     }
		}

		else {   //use the kinematics of the reco'ed particle
		     leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
		}

            reso_lv += leg_lv;
          }
      }
      reso.momentum.x = reso_lv.Px();
      reso.momentum.y = reso_lv.Py();
      reso.momentum.z = reso_lv.Pz();
      reso.mass = reso_lv.M();
      result.emplace_back(reso);
    } while (std::next_permutation(v.begin(), v.end()));
  }
  if (result.size() > 1) {
    auto resonancesort = [&] (edm4hep::ReconstructedParticleData i ,edm4hep::ReconstructedParticleData j) { return (abs( m_resonance_mass -i.mass)<abs(m_resonance_mass-j.mass)); };
		std::sort(result.begin(), result.end(), resonancesort);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator first = result.begin();
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator last = result.begin() + 1;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> onlyBestReso(first, last);
    return onlyBestReso;
  } else {
    return result;
  }
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  HiggsTools::sort_greater(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
	ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  //at least one muon + and one muon - in each event
	
  int n = in.size();
  if (n == 0 ){
    return result;
  }
  ROOT::VecOps::RVec<float> pT = FCCAnalyses::ReconstructedParticle::get_pt(in);
  std::vector< std::pair <float, edm4hep::ReconstructedParticleData> > vect;
  for (int i = 0; i<n ; ++i){
    vect.push_back(std::make_pair(pT.at(i),in.at(i)));
  }
  std::stable_sort(vect.begin(), vect.end(),
                 [](const auto& a, const auto& b){return a.first > b.first;});
  //std::sort(vect.begin(), vect.end());
  for (int i = 0; i<n ; ++i){
    result.push_back(vect.at(i).second);
  } 
  return result;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  HiggsTools::sort_greater_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
	ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  //at least one muon + and one muon - in each event
	
  int n = in.size();
  if (n == 0 ){
    return result;
  }

  ROOT::VecOps::RVec<float> p = FCCAnalyses::ReconstructedParticle::get_p(in);
  std::vector< std::pair <float, edm4hep::ReconstructedParticleData> > vect;
  for (int i = 0; i<n ; ++i){
    vect.push_back(std::make_pair(p.at(i),in.at(i)));
  }
  std::stable_sort(vect.begin(), vect.end(),
                 [](const auto& a, const auto& b){return a.first > b.first;});
  //std::sort(vect.begin(), vect.end());
  for (int i = 0; i<n ; ++i){
    result.push_back(vect.at(i).second);
  } 
  return result;
}

float HiggsTools::Reweighting_wzp_kkmc(float pT, float m) {
  float scale;
  if (m > 220.){
    if ( pT > 0. && pT <= 1.){scale = 1.0032322;}
    else if ( pT > 1. && pT <= 2.){scale = 0.95560480;}
    else if ( pT > 2. && pT <= 3.){scale = 0.94986398;}
    else if ( pT > 3. && pT <= 4.){scale = 0.95134129;}
    else if ( pT > 4. && pT <= 5.){scale = 0.94456404;}
    else if ( pT > 5. && pT <= 6.){scale = 0.94726464;}
    else if ( pT > 6. && pT <= 7.){scale = 0.94101542;}
    else if ( pT > 7. && pT <= 8.){scale = 0.91753618;}
    else if ( pT > 8. && pT <= 9.){scale = 0.91804730;}
    else if ( pT > 9. && pT <= 10.){scale = 0.92097238;}
    else if ( pT > 10. && pT <= 11.){scale = 0.91521958;}
    else if ( pT > 11. && pT <= 12.){scale = 0.94550474;}
    else if ( pT > 12. && pT <= 13.){scale = 0.91403417;}
    else if ( pT > 13. && pT <= 14.){scale = 0.87701843;}
    else if ( pT > 14. && pT <= 15.){scale = 0.89537075;}
    else if ( pT > 15. && pT <= 16.){scale = 0.90811988;}
    else if ( pT > 16. && pT <= 17.){scale = 0.90657018;}
    else if ( pT > 17. && pT <= 18.){scale = 0.93739754;}
    else if ( pT > 18. && pT <= 19.){scale = 0.98795371;}
    else if ( pT > 19. && pT <= 20.){scale = 2.6656045;}
    else {scale = 1.;}
  }
  else {scale = 1.;}

  return scale;
}

ROOT::VecOps::RVec<float> HiggsTools::acolinearity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<float> result;
  if(in.size() < 2) return result;

  TLorentzVector p1;
  p1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);

  TLorentzVector p2;
  p2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);

  float acol = abs(p1.Theta() - p2.Theta());

  result.push_back(acol);
  return result;
}

ROOT::VecOps::RVec<float> HiggsTools::acoplanarity(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in){
  ROOT::VecOps::RVec<float> result;
  if(in.size() < 2) return result;

  TLorentzVector p1;
  p1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);

  TLorentzVector p2;
  p2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);

  float acop = abs(p1.Phi() - p2.Phi());
  if(acop > M_PI) acop = 2 * M_PI - acop;
  acop = M_PI - acop;

  result.push_back(acop);
  return result;
}


// perturb the scale of the particles
HiggsTools::momentum_scale::momentum_scale(float arg_scaleunc) : scaleunc(arg_scaleunc) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  HiggsTools::momentum_scale::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    /*
     TLorentzVector lv;
     lv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
     lv *= (1. + scaleunc);
     p.momentum.x = lv.Px();
     p.momentum.y = lv.Py();
     p.momentum.z = lv.Pz();
     //p.energy = lv.E();
    */
    p.momentum.x = p.momentum.x*(1. + scaleunc);
    p.momentum.y = p.momentum.y*(1. + scaleunc);
    p.momentum.z = p.momentum.z*(1. + scaleunc);
    result.emplace_back(p);
  }
  return result;
}

/// to be added to your ReconstructedParticle.cc
sel_type::sel_type( int arg_pdg, bool arg_chargeconjugate) : m_pdg(arg_pdg), m_chargeconjugate( arg_chargeconjugate )  {};

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  HiggsTools::sel_type::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  result.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if ( m_chargeconjugate ) {
        if ( std::abs( p.type ) == std::abs( m_pdg)  ) result.emplace_back(p);
    }
    else {
        if ( p.type == m_pdg ) result.emplace_back(p);
    }
  }
  return result;
}

ROOT::VecOps::RVec<float> Isolation( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> particles, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) {

  ROOT::VecOps::RVec<float> result;

  float DeltaRMax = 0.5 ;
  float PTMin = 0.5 ;

  for (size_t i = 0; i < particles.size(); ++i) {

    auto & par1 = in[i];
    float pt1 = sqrt( pow(par1.momentum.x, 2) + pow(par1.momentum.y, 2) ) ;

    TVector3 p1(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z );
    float sum = 0;

    for (size_t j = 0; j < in.size(); ++j) {
      auto & par2 = in[j];
      if ( par1.energy == par2.energy && par1.momentum.x == par2.momentum.x && par1.momentum.y == par2.momentum.y && par1.momentum.z == par2.momentum.z ) continue;

      TVector3 p2( in.at(j).momentum.x, in.at(j).momentum.y, in.at(j).momentum.z );
      if ( sqrt( pow(par2.momentum.x, 2) + pow(par2.momentum.y, 2) ) < PTMin ) continue;
      float delta_ij = fabs( p1.DeltaR( p2 ) );
      if ( delta_ij < DeltaRMax)  sum +=  sqrt( pow(par2.momentum.x, 2) + pow(par2.momentum.y, 2) ) ;

    }

    float isolation = sum / pt1 ;
    result.push_back( isolation );
  }

  return result;
}

sel_isol::sel_isol( float arg_isocut ) : m_isocut (arg_isocut) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::sel_isol::operator() (  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> particles, ROOT::VecOps::RVec<float> var ) { 
  
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  for (size_t i=0; i < particles.size(); ++i) {
    auto & p = particles[i];
    if ( var[i] < m_isocut) result.emplace_back( p );
  }

  return result;
}

BoostAngle::BoostAngle( float arg_angle ) : m_angle ( arg_angle ) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::BoostAngle::operator() ( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  float ta =  tan( m_angle ) ;
  for ( size_t i=0; i < in.size(); ++i) {
    auto & p = in[i];
    edm4hep::ReconstructedParticleData newp = p;
    float e = p.energy ;
    float px = p.momentum.x;
    float e_prime = e * sqrt( 1 + ta*ta ) + px * ta ;
    float px_prime = px * sqrt( 1 + ta*ta ) + e * ta ;
    newp.momentum.x = px_prime;
    newp.energy = e_prime;
    result.push_back( newp );
  }
  return result;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> Merger( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) {
  // cf Delphes Merger module. Returns a vector with one single RecoParticle corresponding to the global
  // sum (cf the "MissingET" object in the Delphes files )
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  edm4hep::ReconstructedParticleData sum;
  float px=0;
  float py=0;
  float pz=0;
  float energy=0;

  for (size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    px = px + p.momentum.x;
    py = py + p.momentum.y;
    pz = pz + p.momentum.z;
    energy = energy + p.energy ;
  }

  sum.momentum.x = px;
  sum.momentum.y = py;
  sum.momentum.z = pz;
  sum.energy = energy;
  sum.charge = 0;
  sum.mass = 0;

  result.push_back( sum );
  return result;
}


HiggsTools::resonanceZBuilderHiggsPairs::resonanceZBuilderHiggsPairs(float arg_resonance_mass, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::resonanceZBuilderHiggsPairs::resonanceZBuilderHiggsPairs::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
				ROOT::VecOps::RVec<int> recind,
				ROOT::VecOps::RVec<int> mcind,
				ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                ROOT::VecOps::RVec<int> parents,
                ROOT::VecOps::RVec<int> daugthers)   {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    int n = legs.size();
    std::vector<bool> hDecay;
    
    //cout << "*** BUILD RESO ***" << endl;
    //cout << "Number of leptons: " << n << endl;
  
    if (n >1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            edm4hep::ReconstructedParticleData reso;
            //set initial charge == 0
            reso.charge = 0;
            TLorentzVector reso_lv;
            bool muonFromHiggsDecay = false;
            bool oneMuonFromHiggsDecay = false;
            for (int i = 0; i < n; ++i) {
                if (v[i]) {

                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    // Ideal detector resolution: use the kinematics of the MC particle instead
                    if ( m_use_MC_Kinematics) {

                         // ugly: particles_begin is not filled in RecoParticle.
                         // hence: either need to keep trace of the index of the legs into the RecoParticle collection,
                         // or (as done below) use the track index to map the leg to the MC particle :-(

                         int track_index = legs[i].tracks_begin ;   // index in the Track array
                         int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index( track_index, recind, mcind, reco );
                         if ( mc_index >= 0 && mc_index < mc.size() ) {
                         int pdgID = mc.at( mc_index).PDG;
                             leg_lv.SetXYZM(mc.at(mc_index ).momentum.x, mc.at(mc_index ).momentum.y, mc.at(mc_index ).momentum.z, mc.at(mc_index ).mass );
                         }
                    }

                    else {   //use the kinematics of the reco'ed particle
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }
                    

                
                    
                    // find the Higgs MC particle
                    ROOT::VecOps::RVec<int> higgsParticle = gen_sel_pdgIDInt(25, false)(mc);
                    if(higgsParticle.size() > 0) {
                   
                    
                        //std::vector<int> tmp = gen_decay_list(higgsParticle, mc, daugthers);
                        //if(std::abs(tmp.at(0)) == 23) {
                        
                        int track_index = legs[i].tracks_begin ;   // index in the Track array
                        int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco); // MC index of the muon
                        
                        muonFromHiggsDecay = from_Higgsdecay(mc_index, mc, parents);
                        if(muonFromHiggsDecay) oneMuonFromHiggsDecay = true;
                        //if(muonFromHiggsDecay) std::cout << muonFromHiggsDecay << std::endl;
                        
                        //}
                    }
                    reso_lv += leg_lv;
                }
            }
      
      
            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            //if(oneMuonFromHiggsDecay == true) continue; // neglect wrong paired muons
            
            hDecay.push_back(oneMuonFromHiggsDecay);
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);

        } while (std::next_permutation(v.begin(), v.end()));
    }
  
    if (result.size() > 1) {
  
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> bestReso;
  
        // sort on mZ
        //auto resonancesort = [&] (edm4hep::ReconstructedParticleData i ,edm4hep::ReconstructedParticleData j) { return (abs( m_resonance_mass -i.mass)<abs(m_resonance_mass-j.mass)); };
		//std::sort(result.begin(), result.end(), resonancesort);
        //ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator first = result.begin();
        //ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>::const_iterator last = result.begin() + 1;
        //ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> onlyBestReso(first, last);
        //cout << " ->" << onlyBestReso.size() << " " << result.size() << endl;
        //return onlyBestReso;
        
        // sort on recoil
        
        //cout << "*** PAIR SELECTOR ***" << endl;
  
        
        int idx_min = -1;
        float d_min = 9e9;
        //cout << "-------------" << endl;
        for (int i = 0; i < result.size(); ++i) {
            
            //if(hDecay.at(i)) continue;
     
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, 240);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - 91.2, 2); // mass
            float rec = std::pow(recoil_fcc.mass - 125.0, 2); // recoil
            float d = 0.5*mass + 0.5*rec;
            d = mass;
            
            //cout << " idx=" << i << "  mZ = "<< result.at(i).mass << " mRec = " << recoil_fcc.mass << " oneFromHiggs=" << hDecay.at(i) << endl;
            
            //cout << " one muon from higgs = " << hDecay.at(i) << " mZ = "<< result.at(i).mass << " mRec = " << recoil_fcc.mass << endl;
            
            // MC constrained
            /*
            if(hDecay.at(i) == 0) {
                idx_min = i;
                break;
            }
            */
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }
     
        }
     
        //cout << " nReso=" << result.size() << " mZ=" << result.at(idx_min).mass << endl;
        //if(hDecay.at(idx_min) == 1) cout << " nReso=" << result.size() << " mZ=" << result.at(idx_min).mass << endl;
     
        //cout << " -> selected idx=" << idx_min << " oneFromHiggs=" << hDecay.at(idx_min) << endl;
        if(idx_min > -1) bestReso.push_back(result.at(idx_min));
        return bestReso;
    }
    else {

       // if(result.size() > 0 and hDecay.at(0)) { // return empty if one muon comes from the Higgs decay
        //    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> bestReso;
       //     return bestReso;
        //}
        
        //if(result.size() > 0 and hDecay.at(0) == 1) cout << " nReso=" << result.size() << " mZ=" << result.at(0).mass << endl;
        return result;
    }
}    

HiggsTools::resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::resonanceBuilder_mass_recoil::resonanceBuilder_mass_recoil::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
        ROOT::VecOps::RVec<int> recind,
				ROOT::VecOps::RVec<int> mcind,
				ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                ROOT::VecOps::RVec<int> parents,
                ROOT::VecOps::RVec<int> daugthers)   {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            edm4hep::ReconstructedParticleData reso;
            reso.charge = 0;
            TLorentzVector reso_lv;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
        exit(1);
    }
  
    if(result.size() > 1) {
  
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> bestReso;
        
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float d = mass + chi2_recoil_frac*rec;
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }
     
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}
HiggsTools::resonanceBuilder_mass_recoil2::resonanceBuilder_mass_recoil2(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::resonanceBuilder_mass_recoil2::resonanceBuilder_mass_recoil2::operator()(
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
    ROOT::VecOps::RVec<int> recind,
    ROOT::VecOps::RVec<int> mcind,
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
    ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
    ROOT::VecOps::RVec<int> parents,
    ROOT::VecOps::RVec<int> daugthers)
{
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the leptons
    int n = legs.size();

    if(n > 1) {
        /*
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        
        int idx_1 = -1;
        int idx_2 = -1;
        do {
            std::vector<int> pair;
            edm4hep::ReconstructedParticleData reso;
            reso.charge = 0;
            TLorentzVector reso_lv;
            int idx_1_ = -1;
            int idx_2_ = -1;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    if(idx_1_ == -1) idx_1_ = i;
                    else idx_2_ = i;
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y,
                                           mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                        leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y,
                                       legs[i].momentum.z, legs[i].mass);
                    }
                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            if(reso_lv.M() < (m_recoil_mass+1) && reso_lv.M() > (m_recoil_mass-1)) {
                idx_1 = idx_1_;
                idx_2 = idx_2_;
                break;
            }
        } while(std::next_permutation(v.begin(), v.end()));
        */

        ROOT::VecOps::RVec<bool> w(n);
        std::fill(w.end() - 2, w.end(), true); // helper variable for permutations

        do {
            std::vector<int> pair;
            edm4hep::ReconstructedParticleData reso;
            reso.charge = 0;
            TLorentzVector reso_lv; 
            for(int i = 0; i < n; ++i) {
                // if(i==idx_1 || i==idx_2) continue;
                if(w[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y,
                                           mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                        leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y,
                                       legs[i].momentum.z, legs[i].mass);
                    }
                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            // if(reso_lv.M() < (m_recoil_mass+5) && reso_lv.M() > (m_recoil_mass-5)) continue; // avoid Higgs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);
        } while(std::next_permutation(w.begin(), w.end()));
    }
    else {
        // Return dummy objects when there are insufficient legs.
        auto dummy = edm4hep::ReconstructedParticleData();
        dummy.momentum.x = 0;
        dummy.momentum.y = 0;
        dummy.momentum.z = 0;
        dummy.mass = 0;
        result.emplace_back(dummy);
        result.emplace_back(dummy);
        result.emplace_back(dummy);
        return result;
    }

    if(result.size() > 1) {
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> bestReso;
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            // Calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y,
                        result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;

            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();

            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y,
                       result.at(i).momentum.z, result.at(i).mass);

            float mass_term = std::pow(result.at(i).mass - m_resonance_mass, 2); // resonance mass deviation
            float recoil_term = std::pow(recoil_fcc.mass - m_recoil_mass, 2);    // recoil mass deviation
            float d = (1.0 - chi2_recoil_frac) * mass_term + chi2_recoil_frac * recoil_term;

            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found. RETURN DUMMY." << std::endl;
            auto dummy = edm4hep::ReconstructedParticleData();
            dummy.momentum.x = 0;
            dummy.momentum.y = 0;
            dummy.momentum.z = 0;
            dummy.mass = 0;
            result.emplace_back(dummy);
            result.emplace_back(dummy);
            result.emplace_back(dummy);
            return result;
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}


gen_sel_pdgIDInt::gen_sel_pdgIDInt(int arg_pdg, bool arg_chargeconjugate) : m_pdg(arg_pdg), m_chargeconjugate( arg_chargeconjugate )  {};

ROOT::VecOps::RVec<int>  HiggsTools::gen_sel_pdgIDInt::gen_sel_pdgIDInt::operator() (ROOT::VecOps::RVec<edm4hep::MCParticleData> in) {
  ROOT::VecOps::RVec<int> result;
  for(size_t i = 0; i < in.size(); ++i) {
    auto & p = in[i];
    if(m_chargeconjugate) {
      if(std::abs( p.PDG) == std::abs(m_pdg)) result.push_back(i);
    }
    else {
      if(p.PDG == m_pdg) result.push_back(i);
    }
  }
  return result;
}

// for a given MC index, it returns whether or not one of these muons come (indirectly) from a Higgs decay
bool HiggsTools::from_Higgsdecay(int i, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {

  bool ret = false;
  std::vector<int> res;
  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the parents, Particle#0.index

  // returns whether the particle i comes from the chain containing the Higgs
  if ( i < 0 || i >= in.size() ) return ret;
  int db = in.at(i).parents_begin;
  int de = in.at(i).parents_end;

  //std::cout << "Chain for " << in.at(i).PDG << std::endl;
  //std::cout << "Chain for " << in.at(i).PDG << std::endl;
  //std::cout << "Chain for idx=" << i << " with PDG=" << in.at(i).PDG << " having db=" << db << " and de=" << de << std::endl;

  if(db == de) return false; // top of tree

  for(int id = db; id < de; id++) { // loop over all parents
    int iparent = ind.at(id);
    //std::cout << " Analyze parent idx=" << iparent << " PDG=" << in.at(iparent).PDG << std::endl;
    if(std::abs(in.at(iparent).PDG) == 25) ret = true; // if Higgs is found
    else ret = from_Higgsdecay(iparent, in, ind); // otherwise go up in the decay tree
  }
  return ret;
}

// for a given muon collection (legs), it returns whether or not one of these muons come (indirectly) from a Higgs decay
bool HiggsTools::from_Higgsdecay(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents, ROOT::VecOps::RVec<int> daugther) {
      
  bool ret = false;
  for (size_t i = 0; i < legs.size(); ++i) {
                    
    int track_index = legs[i].tracks_begin;
    int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
    if(from_Higgsdecay(mc_index, mc, parents)) {
      ret = true;
      break;
    }
  }            
  return ret;
}

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };

ROOT::VecOps::RVec<double>  HiggsTools::coneIsolation::coneIsolation::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                                 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rp) {
    
  ROOT::VecOps::RVec<double> result;
  result.reserve(recop.size());

  std::vector<TLorentzVector> lv_reco;
  std::vector<TLorentzVector> lv_charged;
  std::vector<TLorentzVector> lv_neutral;

  for(size_t i = 0; i < rp.size(); ++i) {

    TLorentzVector tlv;
    tlv.SetPxPyPzE(rp.at(i).momentum.x, rp.at(i).momentum.y, rp.at(i).momentum.z, rp.at(i).energy);
            
    if(rp.at(i).charge == 0) lv_neutral.push_back(tlv);
    else lv_charged.push_back(tlv);
  }
  
  for(size_t i = 0; i < recop.size(); ++i) {

    TLorentzVector tlv;
    tlv.SetPxPyPzE(recop.at(i).momentum.x, recop.at(i).momentum.y, recop.at(i).momentum.z, recop.at(i).energy);
    lv_reco.push_back(tlv);
  }
  
  // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
  for (auto & lv_reco_ : lv_reco) {
        
    double sumNeutral = 0.0;
    double sumCharged = 0.0;
    
    // charged
    for (auto & lv_charged_ : lv_charged) {
      double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
      if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
    }

    // neutral
    for (auto & lv_neutral_ : lv_neutral) {

      double dr = coneIsolation::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
      if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
    }

    double sum = sumCharged + sumNeutral;
    double ratio= sum / lv_reco_.P();
    result.emplace_back(ratio);
  }
  return result;
}

coneIsolationTheta::coneIsolationTheta(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };

ROOT::VecOps::RVec<double>  HiggsTools::coneIsolationTheta::coneIsolationTheta::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> recop,
                                 ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rp) {

  ROOT::VecOps::RVec<double> result;
  result.reserve(recop.size());

  std::vector<TLorentzVector> lv_reco;
  std::vector<TLorentzVector> lv_charged;
  std::vector<TLorentzVector> lv_neutral;

  for(size_t i = 0; i < rp.size(); ++i) {

    TLorentzVector tlv;
    tlv.SetPxPyPzE(rp.at(i).momentum.x, rp.at(i).momentum.y, rp.at(i).momentum.z, rp.at(i).energy);

    if(rp.at(i).charge == 0) lv_neutral.push_back(tlv);
    else lv_charged.push_back(tlv);
  }

  for(size_t i = 0; i < recop.size(); ++i) {

    TLorentzVector tlv;
    tlv.SetPxPyPzE(recop.at(i).momentum.x, recop.at(i).momentum.y, recop.at(i).momentum.z, recop.at(i).energy);
    lv_reco.push_back(tlv);
  }

  // identical to coneIsolation, except the cone distance uses polar angle (theta) instead of eta
  for (auto & lv_reco_ : lv_reco) {

    double sumNeutral = 0.0;
    double sumCharged = 0.0;

    // charged
    for (auto & lv_charged_ : lv_charged) {
      double dr = coneIsolationTheta::deltaR(lv_reco_.Theta(), lv_reco_.Phi(), lv_charged_.Theta(), lv_charged_.Phi());
      if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
    }

    // neutral
    for (auto & lv_neutral_ : lv_neutral) {

      double dr = coneIsolationTheta::deltaR(lv_reco_.Theta(), lv_reco_.Phi(), lv_neutral_.Theta(), lv_neutral_.Phi());
      if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
    }

    double sum = sumCharged + sumNeutral;
    double ratio= sum / lv_reco_.P();
    result.emplace_back(ratio);
  }
  return result;
}

std::vector<float> HiggsTools::gen_p_from_reco(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {

  std::vector<float> result;

  for (size_t i = 0; i < legs.size(); ++i) {
           
    int track_index = legs[i].tracks_begin;
    int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                  
    TLorentzVector leg_lv;
    if(mc_index >= 0 && mc_index < mc.size() ) {
      leg_lv.SetXYZM(mc.at(mc_index ).momentum.x, mc.at(mc_index ).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index ).mass);
      result.push_back(leg_lv.P());
    }
    else {
      std::cout << "MC track not found!" << std::endl;
    }
  }
  return result;
}

ROOT::VecOps::RVec<edm4hep::MCParticleData> HiggsTools::get_photons(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc) {

  ROOT::VecOps::RVec<edm4hep::MCParticleData> result;

  for(size_t i = 0; i < mc.size(); ++i) {
                 
    auto & p = mc[i];
    if(p.PDG == 22) result.emplace_back(p);
  }
  return result;
}

ROOT::VecOps::RVec<float> HiggsTools::leptonResolution_p(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons, ROOT::VecOps::RVec<int> recind,
                                                         ROOT::VecOps::RVec<int> mcind,
                                                         ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
                                                         ROOT::VecOps::RVec<edm4hep::MCParticleData> mc){
      
  ROOT::VecOps::RVec<float> result;
  result.reserve(muons.size());
      
  for(int i = 0; i < muons.size(); ++i) {

    TLorentzVector reco_;
    reco_.SetXYZM(muons[i].momentum.x, muons[i].momentum.y, muons[i].momentum.z, muons[i].mass);
    int track_index = muons[i].tracks_begin;
    int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
    if(mc_index >= 0 && mc_index < (int)mc.size()) {
      TLorentzVector mc_;
      mc_.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
      if(mc_.P() > 20) result.push_back(reco_.P()/mc_.P());
    }
  } 
  return result;
}

std::vector<int> HiggsTools::gen_decay_list(ROOT::VecOps::RVec<int> mcin, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {

  std::vector<int> result;


  // i = index of a MC particle in the Particle block
  // in = the Particle collection
  // ind = the block with the indices for the daughters, Particle#1.index
  
  // returns a vector with the indices (in the Particle block) of the daughters of the particle i
  
  for (size_t i = 0; i < mcin.size(); ++i) {
    for (size_t j = 0; j < FCCAnalyses::MCParticle::get_list_of_particles_from_decay(mcin[i],in,ind).size(); ++j) {
      if(in[FCCAnalyses::MCParticle::get_list_of_particles_from_decay(mcin[i],in,ind)[j]].PDG != 25) {
        result.push_back(in[FCCAnalyses::MCParticle::get_list_of_particles_from_decay(mcin[i], in, ind)[j]].PDG);
       }
    }
  }
  return result;
}

float HiggsTools::Higgsstrahlungness(float mll, float mrecoil) {
  float mZ = 91.2;
  float mH = 125.;
  float chi2_recoil_frac = 0.4;
  float chiZ = std::pow(mll - mZ, 2); // mass
  float chiH = std::pow(mrecoil - mH, 2); // recoil
  float chi2 = (1.0-chi2_recoil_frac)*chiZ + chi2_recoil_frac*chiH;
  return chi2;
}

// calculate the number of foward leptons
HiggsTools::polarAngleCategorization::polarAngleCategorization(float arg_thetaMin, float arg_thetaMax) : thetaMin(arg_thetaMin), thetaMax(arg_thetaMax) {};
int polarAngleCategorization::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    
    int nFwd = 0; // number of forward leptons
    for (size_t i = 0; i < in.size(); ++i) {
        
        auto & p = in[i];
        TLorentzVector lv;
        lv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        if(lv.Theta() < thetaMin || lv.Theta() > thetaMax) nFwd += 1;
    }
    return nFwd;
}

// deltaR between two reco particles, based on eta
float HiggsTools::deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;
    
    TLorentzVector tlv1;
    tlv1.SetPxPyPzE(in.at(0).momentum.x, in.at(0).momentum.y, in.at(0).momentum.z, in.at(0).energy);

    TLorentzVector tlv2;
    tlv2.SetPxPyPzE(in.at(1).momentum.x, in.at(1).momentum.y, in.at(1).momentum.z, in.at(1).energy);

    return tlv1.DeltaR(tlv2); 
}

float HiggsTools::deltaRPrime(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;
    
    TLorentzVector tlv1;
    tlv1.SetPxPyPzE(in.at(0).momentum.x, in.at(0).momentum.y, in.at(0).momentum.z, in.at(0).energy);

    TLorentzVector tlv2;
    tlv2.SetPxPyPzE(in.at(1).momentum.x, in.at(1).momentum.y, in.at(1).momentum.z, in.at(1).energy);

    float dTheta = abs(tlv1.Theta()-tlv2.Theta());
    float dPhi = abs(tlv1.DeltaPhi(tlv2));
    return std::sqrt(dTheta*dTheta + dPhi*dPhi);
}

// deltaR between two reco particles, based on eta
float HiggsTools::deltaPhi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;
    
    TLorentzVector tlv1;
    tlv1.SetPxPyPzE(in.at(0).momentum.x, in.at(0).momentum.y, in.at(0).momentum.z, in.at(0).energy);

    TLorentzVector tlv2;
    tlv2.SetPxPyPzE(in.at(1).momentum.x, in.at(1).momentum.y, in.at(1).momentum.z, in.at(1).energy);
    
    return abs(tlv1.DeltaPhi(tlv2));
}

// deltaR between two reco particles, based on eta
float HiggsTools::deltaTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;
    
    TLorentzVector tlv1;
    tlv1.SetPxPyPzE(in.at(0).momentum.x, in.at(0).momentum.y, in.at(0).momentum.z, in.at(0).energy);

    TLorentzVector tlv2;
    tlv2.SetPxPyPzE(in.at(1).momentum.x, in.at(1).momentum.y, in.at(1).momentum.z, in.at(1).energy);
    
    return abs(tlv1.Theta()-tlv2.Theta());
}

// deltaR between two reco particles, based on eta
float HiggsTools::deltaEta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() != 2) return -1;
    
    TLorentzVector tlv1;
    tlv1.SetPxPyPzE(in.at(0).momentum.x, in.at(0).momentum.y, in.at(0).momentum.z, in.at(0).energy);

    TLorentzVector tlv2;
    tlv2.SetPxPyPzE(in.at(1).momentum.x, in.at(1).momentum.y, in.at(1).momentum.z, in.at(1).energy);

    return abs(tlv1.Eta()-tlv2.Eta()); 
}

// minimum between particles
float HiggsTools::min_deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float min = 9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = tlv1.DeltaR(tlv2);
            if(value < min) min = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaR, at least two leptons required." << std::endl;
        exit(1);
    }
    return min;
}

// minimum between particles
float HiggsTools::min_deltaRPrime(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float min = 9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float dTheta = abs(tlv1.Theta()-tlv2.Theta());
            float dPhi = abs(tlv1.DeltaPhi(tlv2));
            float value = sqrt(dTheta*dTheta + dPhi*dPhi);
            if(value < min) min = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaR, at least two leptons required." << std::endl;
        exit(1);
    }
    return min;
}

// minimum between particles
float HiggsTools::min_deltaPhi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float min = 9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = abs(tlv1.DeltaPhi(tlv2)); 
            if(value < min) min = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaPhi, at least two leptons required." << std::endl;
        exit(1);
    }
    return min;
}

// minimum between particles
float HiggsTools::min_deltaTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float min = 9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = std::abs(tlv1.Theta()-tlv2.Theta());
            if(value < min) min = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaTheta, at least two leptons required." << std::endl;
        exit(1);
    }
    return min;
}

// minimum between particles
float HiggsTools::min_deltaEta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float min = 9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = std::abs(tlv1.Eta()-tlv2.Eta());
            if(value < min) min = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaTheta, at least two leptons required." << std::endl;
        exit(1);
    }
    return min;
}


//////
// minimum between particles
float HiggsTools::max_deltaR(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float max = -9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = tlv1.DeltaR(tlv2);
            if(value > max) max > value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaR, at least two leptons required." << std::endl;
        exit(1);
    }
    return max;
}

// minimum between particles
float HiggsTools::max_deltaRPrime(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float max = -9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float dTheta = abs(tlv1.Theta()-tlv2.Theta());
            float dPhi = abs(tlv1.DeltaPhi(tlv2));
            float value = sqrt(dTheta*dTheta + dPhi*dPhi);
            if(value > max) max = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaR, at least two leptons required." << std::endl;
        exit(1);
    }
    return max;
}

// minimum between particles
float HiggsTools::max_deltaPhi(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float max = -9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = abs(tlv1.DeltaPhi(tlv2)); 
            if(value > max) max = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaPhi, at least two leptons required." << std::endl;
        exit(1);
    }
    return max;
}

// minimum between particles
float HiggsTools::max_deltaTheta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float max = -9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = std::abs(tlv1.Theta()-tlv2.Theta());
            if(value > max) max = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaTheta, at least two leptons required." << std::endl;
        exit(1);
    }
    return max;
}

// minimum between particles
float HiggsTools::max_deltaEta(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    if(in.size() < 2) return -1;
    
    float result;
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = in.size();
    float max = -9999;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            TLorentzVector tlv1;
            TLorentzVector tlv2;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                }
            }
            tlv1.SetXYZM(in[pair[0]].momentum.x, in[pair[0]].momentum.y, in[pair[0]].momentum.z, in[pair[0]].mass);
            tlv2.SetXYZM(in[pair[1]].momentum.x, in[pair[1]].momentum.y, in[pair[1]].momentum.z, in[pair[1]].mass);
            float value = std::abs(tlv1.Eta()-tlv2.Eta());
            if(value > max) max = value;
            pairs.push_back(pair);
        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: min_deltaTheta, at least two leptons required." << std::endl;
        exit(1);
    }
    return max;
}

// returns the gen particles with given PDGID (absolute) that have the e+/e- as parent, i.e. from prompt
// in Whizard, the prompt leptons from the collision have two parents, the electron and positron
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::whizard_zh_select_prompt_leptons(
                                                        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, 
                                                        ROOT::VecOps::RVec<int> recind, 
                                                        ROOT::VecOps::RVec<int> mcind, 
                                                        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, 
                                                        ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, 
                                                        ROOT::VecOps::RVec<int> parents, 
                                                        ROOT::VecOps::RVec<int> daugther) {
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    for (size_t i = 0; i < in.size(); ++i) {
        int track_index = in[i].tracks_begin;
        int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
        if(HiggsTools::whizard_zh_from_prompt(mc_index, mc, parents)) {
            result.emplace_back(in[i]);
        }
    }
    return result;
}

// for a given MC index, it returns whether or not one of these muons come (indirectly) from a Higgs decay
bool HiggsTools::whizard_zh_from_prompt(int i, ROOT::VecOps::RVec<edm4hep::MCParticleData> in, ROOT::VecOps::RVec<int> ind) {

    bool ret = false;
    // i = index of a MC particle in the Particle block
    // in = the Particle collection
    // ind = the block with the indices for the parents, Particle#0.index

    // returns whether the particle i comes from the chain containing the Higgs

    if ( i < 0 || i >= in.size() ) return ret;

    int db = in.at(i).parents_begin;
    int de = in.at(i).parents_end;
  
    //std::cout << "Chain for " << in.at(i).PDG << std::endl;
    //std::cout << "Chain for " << in.at(i).PDG << std::endl;
    //std::cout << "Chain for idx=" << i << " with PDG=" << in.at(i).PDG << " having db=" << db << " and de=" << de << std::endl;
    

    if(db == de) return true; // top of tree

   
    for(int id = db; id < de; id++) { // loop over all parents

        int iparent = ind.at(id);
        //std::cout << " Analyze parent idx=" << iparent << " PDG=" << in.at(iparent).PDG << std::endl;
        
        //if(std::abs(in.at(iparent).PDG) == 11) ret = true; // if prompt
        if(iparent == 0) return true;
        else if(std::abs(in.at(iparent).PDG) == 25) ret = false; // non prompt, from Higgs decays
        else ret = HiggsTools::whizard_zh_from_prompt(iparent, in, ind); // otherwise go up in the decay tree
    }
    
    return ret;
}

// perturb the momentum scale with a given constant
HiggsTools::lepton_momentum_scale::lepton_momentum_scale(float arg_scaleunc) : scaleunc(arg_scaleunc) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> lepton_momentum_scale::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        p.momentum.x = p.momentum.x*(1. + scaleunc);
        p.momentum.y = p.momentum.y*(1. + scaleunc);
        p.momentum.z = p.momentum.z*(1. + scaleunc);
        result.emplace_back(p);
    }
    return result;
} 

ROOT::VecOps::RVec<edm4hep::MCParticleData> HiggsTools::get_gen_pdg(ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, int pdgId, bool abs= true) {

   ROOT::VecOps::RVec<edm4hep::MCParticleData> result;
   for(size_t i = 0; i < mc.size(); ++i) {
       
        auto & p = mc[i];
        if((abs and std::abs(p.PDG) == pdgId) or (not abs and p.PDG == pdgId)) result.emplace_back(p);
   }
   return result;
}

BoostFrame_gen::BoostFrame_gen( float arg_low_energy) : low_energy (arg_low_energy) {};
ROOT::VecOps::RVec<edm4hep::MCParticleData> HiggsTools::BoostFrame_gen::operator() ( ROOT::VecOps::RVec<edm4hep::MCParticleData> in ) {
  ROOT::VecOps::RVec<edm4hep::MCParticleData> result;

  // Initialize the initial and final beams
  TLorentzVector beam1(0, 0, 120, 120);
  TLorentzVector beam2(0, 0, -120, 120);
  // To approve in the future
  TLorentzVector final_beam1(0, 0, -low_energy, low_energy);
  TLorentzVector final_beam2(0, 0, 480, 480);

  // Calculate the boost vector of the initial and final system
  TVector3 boost_initial = (beam1 + beam2).BoostVector();
  TVector3 boost_final = (final_beam1 + final_beam2).BoostVector();

  // Calculate the net boost vector
  TVector3 net_boost = boost_final - boost_initial;

  beam1.Boost(net_boost);
  beam2.Boost(net_boost);

  for ( size_t i=0; i < in.size(); ++i) {
    auto & p = in[i];
    edm4hep::MCParticleData newp = p;
    TLorentzVector tlv;
    tlv.SetXYZM(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).mass);
    tlv.Boost(net_boost);
    newp.momentum.x = tlv.Px();
    newp.momentum.y = tlv.Py();
    newp.momentum.z = tlv.Pz();
    newp.mass = tlv.M();
    result.push_back( newp );
  }
  return result;
}


BoostFrame::BoostFrame( float arg_low_energy) : low_energy (arg_low_energy) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::BoostFrame::operator() ( ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in ) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;

  // Initialize the initial and final beams
  TLorentzVector beam1(0, 0, 120, 120);
  TLorentzVector beam2(0, 0, -120, 120);
  // To approve in the future
  TLorentzVector final_beam1(0, 0, -low_energy, low_energy);
  TLorentzVector final_beam2(0, 0, 480, 480);

  // Calculate the boost vector of the initial and final system
  TVector3 boost_initial = (beam1 + beam2).BoostVector();
  TVector3 boost_final = (final_beam1 + final_beam2).BoostVector();

  // Calculate the net boost vector
  TVector3 net_boost = boost_final - boost_initial;

  beam1.Boost(net_boost);
  beam2.Boost(net_boost);

  for ( size_t i=0; i < in.size(); ++i) {
    auto & p = in[i];
    edm4hep::ReconstructedParticleData newp = p;
    TLorentzVector tlv;
    tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
    tlv.Boost(net_boost);
    newp.momentum.x = tlv.Px();
    newp.momentum.y = tlv.Py();
    newp.momentum.z = tlv.Pz();
    newp.energy = tlv.E();
    result.push_back( newp );
  }
  return result;
}


HiggsTools::resonanceBuilder_mass_recoil_boosted::resonanceBuilder_mass_recoil_boosted(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::resonanceBuilder_mass_recoil_boosted::resonanceBuilder_mass_recoil_boosted::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs,
        ROOT::VecOps::RVec<int> recind,
				ROOT::VecOps::RVec<int> mcind,
				ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco,
				ROOT::VecOps::RVec<edm4hep::MCParticleData> mc,
                ROOT::VecOps::RVec<int> parents,
                ROOT::VecOps::RVec<int> daugthers)   {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
  
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            edm4hep::ReconstructedParticleData reso;
            reso.charge = 0;
            TLorentzVector reso_lv;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = FCCAnalyses::ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass_recoil, at least two leptons required." << std::endl;
        exit(1);
    }
  
    if(result.size() > 1) {
  
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> bestReso;
        
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 450, 510);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;
      
            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();
            
            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
        
            float boost = tg.P();
            float mass = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float d = mass + chi2_recoil_frac*rec;
            
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }
     
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}

HiggsTools::recoilBuilder_boosted::recoilBuilder_boosted(float arg_sqrts) : m_sqrts(arg_sqrts) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  recoilBuilder_boosted::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in) {
  ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
  auto recoil_p4 = TLorentzVector(0, 0, 450, 510);
  for (auto & v1: in) {
    TLorentzVector tv1;
    tv1.SetXYZM(v1.momentum.x, v1.momentum.y, v1.momentum.z, v1.mass);
    recoil_p4 -= tv1;
  }
  auto recoil_fcc = edm4hep::ReconstructedParticleData();
  recoil_fcc.momentum.x = recoil_p4.Px();
  recoil_fcc.momentum.y = recoil_p4.Py();
  recoil_fcc.momentum.z = recoil_p4.Pz();
  recoil_fcc.mass = recoil_p4.M();
  result.push_back(recoil_fcc);
  return result;
};


ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::missing(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, float ecm, float p_cutoff) {

    TLorentzVector tlv_ecm(0, 0, 0, ecm);
    TLorentzVector P4sum;
    for (auto & p: in) {
      TLorentzVector tlv;
      tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      P4sum += tlv;
    }

    TLorentzVector P4miss;
    P4miss = tlv_ecm - P4sum;
    
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> ret;
    edm4hep::ReconstructedParticleData res;
    res.momentum.x = P4miss.Px();
    res.momentum.y = P4miss.Py();
    res.momentum.z = P4miss.Pz();
    res.energy = P4miss.E();
    res.mass = P4miss.M();
    ret.emplace_back(res);
    return ret;
}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> HiggsTools::visible(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, float p_cutoff) {
    TLorentzVector P4sum;
    for (auto & p: in) {
      TLorentzVector tlv;
      tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
      P4sum += tlv;
    }

    
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> ret;
    edm4hep::ReconstructedParticleData res;
    res.momentum.x = P4sum.Px();
    res.momentum.y = P4sum.Py();
    res.momentum.z = P4sum.Pz();
    res.energy = P4sum.E();
    res.mass = P4sum.M();
    ret.emplace_back(res);
    return ret;
}


float HiggsTools::ZHChi2(float mZ, float mH, float chi2_H_frac) {
  float mZ_real = 91.2;
  float mH_real = 125.;
  float chiZ = std::pow(mZ - mZ_real, 2); // mass
  float chiH = std::pow(mH - mH_real, 2); // recoil
  float chi2 = (1.0-chi2_H_frac)*chiZ + chi2_H_frac*chiH;
  return chi2;
}