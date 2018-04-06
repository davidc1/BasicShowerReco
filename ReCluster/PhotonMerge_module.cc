////////////////////////////////////////////////////////////////////////
// Class:       PhotonMerge
// Plugin Type: producer (art v2_09_06)
// File:        PhotonMerge_module.cc
//
// Generated at Wed Mar 21 07:50:29 2018 by David Caratelli using cetskelgen
// from cetlib version v3_01_03.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

// larsoft data-products
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/Utilities/FindManyInChainP.h"

#include "uboone/BasicShowerReco/TwoDimTools/Linearity.h"
#include "uboone/BasicShowerReco/TwoDimTools/Poly2D.h"

class PhotonMerge;


class PhotonMerge : public art::EDProducer {
public:
  explicit PhotonMerge(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PhotonMerge(PhotonMerge const &) = delete;
  PhotonMerge(PhotonMerge &&) = delete;
  PhotonMerge & operator = (PhotonMerge const &) = delete;
  PhotonMerge & operator = (PhotonMerge &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string fShrProducer, fVtxProducer, fPhotonProducer;

  bool _debug;

  double _wire2cm, _time2cm;
  
  double _vtxW, _vtxT;

  // shower width (opening angle, in degrees)
  double fWidth; 
  // shower length (cm)
  double fShrLen;
  // maximum fraction of shower charge to be added by re-clustering
  double fFracShrQ;
  // maximum slope difference between shower cluster and to-be merged cluster
  double fMaxSlopeAngle;

  // map connecting photon cluster index to linearity object
  std::map< size_t, twodimtools::Linearity > _photon_lin_map; 
  // map connecting photon cluster index to poly2d object
  std::map< size_t, twodimtools::Poly2D > _photon_poly_map; 

  twodimtools::Poly2D projectShower(const art::Ptr<recob::Cluster> clus);

  double slopeCompat(const twodimtools::Poly2D& shr,
		     const twodimtools::Linearity& photon);

  bool photonCrossesShower(const twodimtools::Poly2D& shr,
			   const twodimtools::Poly2D& photon);
  
};


PhotonMerge::PhotonMerge(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces<std::vector<recob::PFParticle> >();
  produces<std::vector<recob::Cluster> >();
  produces<art::Assns <recob::PFParticle, recob::Cluster>    >();
  produces<art::Assns <recob::PFParticle, recob::Hit>        >();

  fShrProducer    = p.get<std::string>("fShrProducer");
  fVtxProducer    = p.get<std::string>("fVtxProducer");
  fPhotonProducer = p.get<std::string>("fPhotonProducer");
  fWidth          = p.get<double>     ("Width");
  fShrLen         = p.get<double>     ("ShrLen");
  fFracShrQ       = p.get<double>     ("FracShrQ");
  fMaxSlopeANgle  = p.get<double>     ("MaxSlopeAngle");
  
  // get detector specific properties
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  _wire2cm = geom->WirePitch(0,1,0);
  _time2cm = detp->SamplingRate() / 1000.0 * detp->DriftVelocity( detp->Efield(), detp->Temperature() );
}

void PhotonMerge::produce(art::Event & e)
{

  // produce recob::PFParticles
  std::unique_ptr< std::vector<recob::PFParticle> > PFParticle_v(new std::vector<recob::PFParticle> );
  std::unique_ptr< std::vector<recob::Cluster>    > Cluster_v   (new std::vector<recob::Cluster>    );
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Cluster>    > PFParticle_Cluster_assn_v( new art::Assns<recob::PFParticle, recob::Cluster>   );
  std::unique_ptr< art::Assns <recob::PFParticle, recob::Hit>        > PFParticle_Hit_assn_v    ( new art::Assns<recob::PFParticle, recob::Hit>       );


  // load input showers
  auto const& shr_h = e.getValidHandle<std::vector<recob::Shower>>(fShrProducer);
  // load input vertices
  auto const& vtx_h = e.getValidHandle<std::vector<recob::Vertex>>(fVtxProducer);
  // load input photon clusters
  auto const& photon_h = e.getValidHandle<std::vector<recob::Cluster>>(fPhotonProducer);
  // grab hits associated associated with photon clusters
  art::FindManyP<recob::Hit> photon_hit_assn_v(photon_h, e, fPhotonProducer);

  // grab clusters associated with shower
  art::FindManyP<recob::Cluster> shr_clus_assn_v(shr_h, e, fShrProducer);
  // grab the hits associated to the showers
  auto shr_hit_assn_v = lar::FindManyInChainP<recob::Hit, recob::Cluster>::find(shr_h, e, fShrProducer);

  if (vtx_h->size() != 1)
    std::cout << "no vertex" << std::endl;
  if (photon_h->size() == 0)
    std::cout << "no photons" << std::endl;

  // load vertex and project on collection-plane
  auto const vtx = vtx_h->at(0);
  Double_t xyz[3] = {};
  vtx.XYZ(xyz);
  auto const* geom = ::lar::providerFrom<geo::Geometry>();
  _vtxW = geom->WireCoordinate(xyz[1],xyz[2],geo::PlaneID(0,0,2)) * _wire2cm;
  _vtxT = xyz[0];

  // create polygon objects for each photon cluster.
  for (size_t p=0; p < photon_h->size(); p++){

    if (_debug) { std::cout << "... new photon" << std::endl; }

    // get assciated hits
    std::vector< art::Ptr<recob::Hit> > photon_hit_v = photon_hit_assn_v.at(p);

    if (photon_hit_v.size() == 0) continue;

    auto plane = photon_hit_v.at(0)->WireID().Plane;

    // we are re-clustering only on the collection-plane
    if ( plane != 2) continue;
    
    // create linearity objects   
    std::vector<double> hit_w_v;
    std::vector<double> hit_t_v;

    // get polygon   
    twodimtools::Poly2D clusPoly(photon_hit_v);

    for (auto hitptr : photon_hit_v){
      hit_w_v.push_back( hitptr->WireID().Wire * _wire2cm );
      hit_t_v.push_back( hitptr->PeakTime() * _time2cm);
    }// for all hits in photon cluster                          

    twodimtools::Linearity clusLin(hit_w_v, hit_t_v);

    _photon_lin_map [ p ] = clusLin;
    _photon_poly_map[ p ] = clusPoly;

    std::cout << "cluster linearity = " << clusLin.linearity() << " with poly area : " << clusPoly.Area() << std::endl;
  }// for all clusters
  
  // loop through reconstructed showers.                                         
  for (size_t s=0; s < shr_h->size(); s++) {
    
    auto const& shr = shr_h->at(s);    
    if (_debug) { std::cout << "new shower" << std::endl; }
    
    // grab collection-plane hits and cluster associated with this shower
    std::vector< art::Ptr<recob::Cluster> > shr_clus_v = shr_clus_assn_v.at(s);
    size_t collidx = 0; // cluster index associated to collection-plane
    for (size_t c=0; c < shr_clus_v.size(); c++) {
      if (shr_clus_v.at(c)->Plane().Plane == 2) {
	collidx = c;
	break;
      }
    }
    std::vector< art::Ptr<recob::Hit> > shr_hit_v = shr_hit_assn_v.at(s);
    std::vector< art::Ptr<recob::Hit> > shr_hit_pl2_v;
    for (auto const& shr_hit : shr_hit_v)
      if (shr_hit->WireID().Plane == 2) { shr_hit_pl2_v.push_back( shr_hit ); } 
    
    auto shrPoly = projectShower(shr_clus_v.at(collidx));
    
    if (_debug) { std::cout << "hits on plane before merging : " << shr_hit_pl2_v.size() << std::endl; }
    
    if (_debug) { std::cout << "polygon has size : " << shrPoly.Size() << std::endl; }

      // if no hits associated on this plane -> skip
    if (shr_hit_pl2_v.size() == 0) continue;
    
    // loop over photons for this plane
    for (std::map<size_t,twodimtools::Poly2D>::iterator it = _photon_poly_map.begin(); it != _photon_poly_map.end(); ++it) {
      
      auto const& photonPoly = it->second;
      
      // check the number of hits in the cluster
      auto nhits = photon_hit_assn_v.at(it->first).size();
      
      if (_debug) { std::cout << "\n\n\t\t new photon with " << nhits << " hits" << std::endl; }
      
      // apply cut on fraction of photon charge (# of hits here) that can be added to existing shower
      if (nhits > (shr_clus_v.at(collidx)->NHits() * fFracShrQ) ) continue;
      
      // get linearity
      auto photonLin = _photon_lin_map[ it->first ];
      
      // do the polygons overlap?
      bool overlap = ( shrPoly.Overlap(photonPoly) || shrPoly.Contained(photonPoly) );
      
      if (overlap == false) {
	if (_debug) std::cout << "\t no overlap w/ shower cone... continue" << std::endl;
	continue;
      }
      
      // for large clusters, add extra checks
      if (nhits > 8) {
	
	// apply cut on slope agreement
	if (slopeCompat( shrPoly, photonLin) > fMaxSlopeAngle) continue;
	
	// make sure the photon doesn't extend on both sides of the shower cone
	if (photonCrossesShower(shrPoly, photonPoly) == true) continue;

      }// if cluster is large
	
      // made it this far -> merge the photon with the shower
      // get set of hits to add (removing potential duplicates)                              
    }// compatible showers     
    
  }// for all showers
  
  e.put(std::move(PFParticle_v));
  e.put(std::move(Cluster_v));
  e.put(std::move(PFParticle_Cluster_assn_v));
  e.put(std::move(PFParticle_Hit_assn_v));

}
  
twodimtools::Poly2D PhotonMerge::projectShower(const art::Ptr<recob::Cluster> clus) {
  
  auto const* detp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  double sW = clus->StartWire() * _wire2cm;
  double sT = ( clus->StartTick() - detp->TriggerOffset() ) * _time2cm;
  
  double dW = (sW-_vtxW);
  double dT = (sT-_vtxT);
  
  double dmag = sqrt((dW*dW)+(dT*dT));
  dW /= dmag;
  dT /= dmag;
  
  double oangle = 0.;
  if (fWidth > 0) oangle = fWidth * 3.14 / 180.;
  
  if (_debug) { std::cout << "direction : dW " << dW << "\t dT : " << dT << std::endl; }
  
  double eW = sW + fShrLen * dW;
  double eT = sT + fShrLen * dT;
  
  // vector of coordiantes with which to consutrct triangle polygon
  std::vector<std::pair<float,float> > triangle_coordinates;
  
  // create polygon vertex for triangle
  triangle_coordinates.push_back( std::pair<float,float>(sW, sT) );
  
  // figure out how far to go on each side.
  double shrWidth = fShrLen * fabs(tan(oangle)) * 0.5;
  
  if (_debug) { std::cout << "width : " << shrWidth << std::endl; }
  
  // unlike length, width is not stretched or compressed on projection.
  // extend end-point "left" and "right" by one width to complete triangle
  // extend in direction perpendicular to line connecting start_pl and end_pl
  double slope = -1. / ( (eT - sT) / (eW - sW) );
  if (_debug)
    {
      std::cout << "\t End pt  = [ "  << eW << ", " << eT << " ]" << std::endl;
      std::cout << "\t SLOPE = " << slope << std::endl;
      std::cout << "\t WIDTH = " << shrWidth << std::endl;
    }
  
  // find triangle base coordinate points 
  double pt1_w = eW + shrWidth * sqrt( 1. / (1 + slope*slope) );
  double pt1_t = eT + slope * shrWidth * sqrt( 1. / (1 + slope*slope) );
  triangle_coordinates.push_back( std::pair<float,float>( pt1_w, pt1_t) );
  
  double pt2_w = eW - shrWidth * sqrt( 1. / (1 + slope*slope) );
  double pt2_t = eT - slope * shrWidth * sqrt( 1. / (1 + slope*slope) );
  triangle_coordinates.push_back( std::pair<float,float>( pt2_w, pt2_t) );
  
  if (_debug){
    std::cout << "SHOWER COORDINATES : " << std::endl
	      << "\t Vertex @ [ " << sW << ", " << sT << " ]" << std::endl
	      << "\t Pt1    @ [ " << pt1_w      << ", " << pt1_t      << " ]" << std::endl
	      << "\t Pt2    @ [ " << pt2_w      << ", " << pt2_t      << " ]" << std::endl;
  }
  
  return twodimtools::Poly2D(triangle_coordinates);
}

double PhotonMerge::slopeCompat(const twodimtools::Poly2D& shr,
				const twodimtools::Linearity& photon) {
  
  
  // get slope of photon
  auto const& photon_slope_val = photon._slope;
  
  // get shower slope
  double shr_slope_val;

  double ydiff = ( (shr.Point(2).second + shr.Point(1).second) / 2. ) - shr.Point(0).second;
  double xdiff = ( (shr.Point(2).first + shr.Point(1).first) / 2. ) - shr.Point(0).first;

  shr_slope_val = ydiff/xdiff;

  double angle = fabs ( atan( ( shr_slope_val - photon_slope_val ) / ( 1 + photon_slope_val * shr_slope_val ) ) );

  angle *= (180./3.14);

  if (_debug) { std::cout << "\t photon slope angle = " << angle << std::endl; }

  return angle;

}

bool PhotonMerge::photonCrossesShower(const twodimtools::Poly2D& shr,
				      const twodimtools::Poly2D& photon) {
  
  // get the triangle points
  double Ox = shr.Point(0).first;
  double Oy = shr.Point(0).second;
  double Ax = shr.Point(1).first;
  double Ay = shr.Point(1).second;
  double Bx = shr.Point(2).first;
  double By = shr.Point(2).second;

  // check where the shower end-point lies
  double Ex = (Ax+Bx)/2.;
  double Ey = (Ay+By)/2.;

  double signEA = (Ex-Ox)*(Ay-Oy) - (Ey-Oy)*(Ax-Ox);
  double signEB = (Ex-Ox)*(By-Oy) - (Ey-Oy)*(Bx-Ox);

  // how many photon edges are on the left or right of the cone?
  int slope_left  = 0;
  int slope_right = 0;

  // being on one side or another of the segments OA and OB depends on the sign
  // of the equations used to compute signEA and signEB
  // to check if points are on either side of the cone compare their sign
  // with the signs of EA and EB knowing that these two are in the center.
  for (size_t i=0; i < photon.Size(); i++) {

    auto const& pt = photon.Point(i);

    double signA = (pt.first-Ox)*(Ay-Oy) - (pt.second-Oy)*(Ax-Ox);
    double signB = (pt.first-Ox)*(By-Oy) - (pt.second-Oy)*(Bx-Ox);

    if ( ((signA*signEA) < 0) && ((signB*signEB) > 0) ) slope_left  += 1;
    if ( ((signA*signEA) > 0) && ((signB*signEB) < 0) ) slope_right += 1;

  }// for all photon points
  
  if ( (slope_left > 0) and (slope_right > 0) ){

    if (_debug) std::cout << "\t photon crosses shower" << std::endl;
    return true;
  }

  return false;
}


void PhotonMerge::beginJob()
{
  // Implementation of optional member function here.
}

void PhotonMerge::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(PhotonMerge)
