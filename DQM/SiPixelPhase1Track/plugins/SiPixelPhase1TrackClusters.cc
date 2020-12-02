// -*- C++ -*-
//
// Package:     SiPixelPhase1TrackClusters
// Class  :     SiPixelPhase1TrackClusters
//

// Original Author: Marcel Schneider

#include "DQM/SiPixelPhase1Common/interface/SiPixelPhase1Base.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelClusterShapeCache.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "RecoTracker/Record/interface/CkfComponentsRecord.h"
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterShapeHitFilter.h"


namespace {

  class SiPixelPhase1TrackClusters final : public SiPixelPhase1Base {
    enum {
      ON_TRACK_CHARGE,
      OFF_TRACK_CHARGE,
      ON_TRACK_BIGPIXELCHARGE,
      OFF_TRACK_BIGPIXELCHARGE,
      ON_TRACK_NOTBIGPIXELCHARGE,
      OFF_TRACK_NOTBIGPIXELCHARGE,

      ON_TRACK_SIZE,
      OFF_TRACK_SIZE,
      OFF_TRACK_SIZEX,
      OFF_TRACK_SIZEY,
      ON_TRACK_SHAPE,

      ON_TRACK_NCLUSTERS,
      OFF_TRACK_NCLUSTERS,

      ON_TRACK_POSITIONB,
      OFF_TRACK_POSITION_B,
      ON_TRACK_POSITIONF,
      OFF_TRACK_POSITION_F,
      OFF_TRACK_POSITION_XZ,
      OFF_TRACK_POSITION_YZ,

      DIGIS_HITMAP_ON_TRACK,
      DIGIS_HITMAP_OFF_TRACK,
      ON_TRACK_NDIGIS,
      OFF_TRACK_NDIGIS,

      //DIGIS_OVER_CLUSTER_TOTCHARGE,
      //DIGIS_OVER_CLUSTER_TOTCHARGE_2D,

      //OFF_TRACK_READOUT_CHARGE,
      //OFF_TRACK_READOUT_NCLUSTERS,
      
      NTRACKS,
      NTRACKS_INVOLUME,

      SIZE_VS_ETA_ON_TRACK_OUTER,
      SIZE_VS_ETA_ON_TRACK_INNER,
      ON_TRACK_CHARGE_OUTER,
      ON_TRACK_CHARGE_INNER,

      ON_TRACK_SHAPE_OUTER,
      ON_TRACK_SHAPE_INNER,

      ON_TRACK_SIZE_X_OUTER,
      ON_TRACK_SIZE_X_INNER,
      ON_TRACK_SIZE_X_F,
      ON_TRACK_SIZE_Y_OUTER,
      ON_TRACK_SIZE_Y_INNER,
      ON_TRACK_SIZE_Y_F,

      ON_TRACK_SIZE_XY_OUTER,
      ON_TRACK_SIZE_XY_INNER,
      ON_TRACK_SIZE_XY_F,
      CHARGE_VS_SIZE_ON_TRACK,

      SIZE_VS_ETA_OFF_TRACK,
      CHARGE_VS_ETA_OFF_TRACK,
      CHARGE_VS_SIZE_OFF_TRACK,
      
      ENUM_SIZE
    };

  public:
    explicit SiPixelPhase1TrackClusters(const edm::ParameterSet& conf);
    //explicit SiPixelPhase1Clusters(const edm::ParameterSet& conf);
    void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
    const bool applyVertexCut_;
    edm::InputTag src_; 

    edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
    edm::EDGetTokenT<reco::VertexCollection> offlinePrimaryVerticesToken_;
    edm::EDGetTokenT<SiPixelClusterShapeCache> pixelClusterShapeCacheToken_;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> pixelSrcToken_;
    edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> tPixelDigi;
    //edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
  };

  
  SiPixelPhase1TrackClusters::SiPixelPhase1TrackClusters(const edm::ParameterSet& iConfig)
      : SiPixelPhase1Base(iConfig), applyVertexCut_(iConfig.getUntrackedParameter<bool>("VertexCut", true)) {
    tracksToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));

    offlinePrimaryVerticesToken_ =
        applyVertexCut_ ? consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))
                        : edm::EDGetTokenT<reco::VertexCollection>();

    pixelClusterShapeCacheToken_ =
        consumes<SiPixelClusterShapeCache>(iConfig.getParameter<edm::InputTag>("clusterShapeCache"));

    pixelSrcToken_ = consumes<edmNew::DetSetVector<SiPixelCluster>>(iConfig.getParameter<edm::InputTag>("clusters"));

    //src_ =  iConfighttps://github.com/cms-analysis/DPGAnalysis-SiPixelTools/blob/master/HitAnalyzer/test/PixDigisTest.cc#L1132.getParameter<edm::InputTag>( "src" );
    tPixelDigi =
      consumes<edm::DetSetVector<PixelDigi>>(iConfig.getParameter<edm::InputTag>("src"));
  }

  void SiPixelPhase1TrackClusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    if (!checktrigger(iEvent, iSetup, DCS))
      return;

    if (histo.size() != ENUM_SIZE) {
      edm::LogError("SiPixelPhase1TrackClusters")
          << "incompatible configuration " << histo.size() << "!=" << ENUM_SIZE << std::endl;
      return;
    }

    // get geometry
    edm::ESHandle<TrackerGeometry> tracker;
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    assert(tracker.isValid());
    auto const& tracker_off = *tracker;

    edm::ESHandle<TrackerTopology> tTopoHandle;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
    auto const& tkTpl = *tTopoHandle;

    edm::ESHandle<ClusterShapeHitFilter> shapeFilterH;
    iSetup.get<CkfComponentsRecord>().get("ClusterShapeHitFilter", shapeFilterH);
    auto const& shapeFilter = *shapeFilterH;

    edm::Handle<reco::VertexCollection> vertices;
    if (applyVertexCut_) {
      iEvent.getByToken(offlinePrimaryVerticesToken_, vertices);
      if (!vertices.isValid() || vertices->empty())
        return;
    }

    //get the map
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);

    if (!tracks.isValid()) {
      edm::LogWarning("SiPixelPhase1TrackClusters") << "track collection is not valid";
      return;
    }

    edm::Handle<SiPixelClusterShapeCache> pixelClusterShapeCacheH;
    iEvent.getByToken(pixelClusterShapeCacheToken_, pixelClusterShapeCacheH);
    if (!pixelClusterShapeCacheH.isValid()) {
      edm::LogWarning("SiPixelPhase1TrackClusters") << "PixelClusterShapeCache collection is not valid";
      return;
    }
    auto const& pixelClusterShapeCache = *pixelClusterShapeCacheH;

    std::unordered_set<const SiPixelCluster*> rHSiPixelClusters;

    //statistics clusters off tracks
    edm::Handle<edmNew::DetSetVector<SiPixelCluster>> inputPixel;
    iEvent.getByToken(pixelSrcToken_, inputPixel);

    if (!inputPixel.isValid())
      return;

    //bool hasClusters = false;

    edmNew::DetSetVector<SiPixelCluster>::const_iterator it;    

    for (it = inputPixel->begin(); it != inputPixel->end(); ++it) {
      auto gid = DetId(it->detId());

      const PixelGeomDetUnit* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(tracker_off.idToDet(gid));
      const PixelTopology& topol = theGeomDet->specificTopology();

      for (SiPixelCluster const& gcluster : *it) {
	//if (!(rHSiPixelClusters.find(&cluster) == rHSiPixelClusters.end())) continue;
	int checks=0;
	bool isON=false;

	for (auto const& track : *tracks) {
	  if (applyVertexCut_ &&
	      (track.pt() < 0.75 || std::abs(track.dxy((*vertices)[0].position())) > 5 * track.dxyError()))
	    continue;

	  bool isBpixtrack = false, isFpixtrack = false;//, crossesPixVol = false;

	  // find out whether track crosses pixel fiducial volume (for cosmic tracks)
	  auto d0 = track.d0(), dz = track.dz();
	  /*	  if (std::abs(d0) < 16 && std::abs(dz) < 50)
	    crossesPixVol = true;
	  */
	  auto etatk = track.eta();

	  auto const& trajParams = track.extra()->trajParams();
	  assert(trajParams.size() == track.recHitsSize());
	  auto hb = track.recHitsBegin();

	  for (unsigned int h = 0; h < track.recHitsSize(); h++) {
	    auto hit = *(hb + h);
	    if (!hit->isValid())
	      continue;
	    auto id = hit->geographicalId();

	    if (id.rawId()!=gid.rawId())
	      continue;

	    // check that we are in the pixel
	    auto subdetid = (id.subdetId());
	    if (subdetid == PixelSubdetector::PixelBarrel)
	      isBpixtrack = true;
	    if (subdetid == PixelSubdetector::PixelEndcap)
	      isFpixtrack = true;
	    if (subdetid != PixelSubdetector::PixelBarrel && subdetid != PixelSubdetector::PixelEndcap)
	      continue;
	    //	    bool iAmBarrel = subdetid == PixelSubdetector::PixelBarrel;

	    // PXB_L4 IS IN THE OTHER WAY
	    // CAN BE XORed BUT LETS KEEP THINGS SIMPLE
	    //	    bool iAmOuter = ((tkTpl.pxbLadder(id) % 2 == 1) && tkTpl.pxbLayer(id) != 4) ||
	    //  ((tkTpl.pxbLadder(id) % 2 != 1) && tkTpl.pxbLayer(id) == 4);

	    auto pixhit = dynamic_cast<const SiPixelRecHit*>(hit->hit());
	    if (!pixhit)
	      continue;

	    //auto geomdetunit = dynamic_cast<const PixelGeomDetUnit*>(pixhit->detUnit());
	    //auto const& topol = geomdetunit->specificTopology();

	    // get the cluster
	    auto clustp = pixhit->cluster();
	    if (clustp.isNull())
	      continue;
	    auto const& cluster = *clustp;

	    if (isBpixtrack or isFpixtrack){
	      if(&gcluster==&cluster){
		std::cout<<"ONTRACK!"<<std::endl;
		std::cout<<"DetID: "<< id.rawId() <<" ("<<gid.rawId()<<")"<<std::endl;
		std::cout<<"SizeX: "<< cluster.sizeX() <<" ("<<gcluster.sizeX()<<")"<<std::endl;
		std::cout<<"SizeY: "<< cluster.sizeY() <<" ("<<gcluster.sizeY()<<")"<<std::endl;
		std::cout<<"Charge: "<< cluster.charge() <<" ("<<gcluster.charge()<<")"<<std::endl;
		isON=true;
		checks++;
		break;
	      }
	      else
		checks++;
	    }
	  }
	  if (isON){
	    std::cout<<"Breaking tracks loop!"<<std::endl;
	    break;
	  }
	}
	std::cout<<"Checks Needed: "<<checks<<std::endl;

      }
    }
    
    histo[ON_TRACK_NCLUSTERS].executePerEventHarvesting(&iEvent);
    histo[ON_TRACK_NDIGIS].executePerEventHarvesting(&iEvent);
    histo[OFF_TRACK_NCLUSTERS].executePerEventHarvesting(&iEvent);
    histo[OFF_TRACK_NDIGIS].executePerEventHarvesting(&iEvent);
  }

}  // namespace

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SiPixelPhase1TrackClusters);
