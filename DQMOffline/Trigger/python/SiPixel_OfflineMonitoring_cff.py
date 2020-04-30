import FWCore.ParameterSet.Config as cms

from DQMOffline.Trigger.SiPixel_OfflineMonitoring_Cluster_cff import *
from DQMOffline.Trigger.SiPixel_OfflineMonitoring_TrackCluster_cff import *
from RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi import *

hltSiPixelClusterShapeCache = siPixelClusterShapeCache.clone(src = 'hltSiPixelClusters')

sipixelMonitorHLTsequence = cms.Sequence(
    hltSiPixelClusterShapeCache +
    hltSiPixelPhase1ClustersAnalyzer +
    hltSiPixelPhase1TrackClustersAnalyzer
)    
