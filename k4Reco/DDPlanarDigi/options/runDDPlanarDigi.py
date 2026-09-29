#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import DDPlanarDigi
from Configurables import GeoSvc
from Configurables import UniqueIDGenSvc
from Configurables import RootHistSvc
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
import os

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["/afs/cern.ch/user/k/kenakkal/k4geo/FCCee/ALFA/compact/ALFA_o1_v00/ALFA_o1_v00.xml"] #[os.environ["K4GEO"]+"/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

digi = DDPlanarDigi()
digi.SubDetectorName = "OTBar" #"Vertex"
digi.IsStrip = False
digi.ResolutionU = [0.003] # or enter same 0.03 16 times for each layer 
digi.ResolutionV = [0.003] # or enter same 0.03 16 times for each layer 
digi.SimTrackHitCollectionName = ["OTBarCollection"]
digi.SimTrkHitRelCollection = ["OTBarHitRelations"]
digi.TrackerHitCollectionName = ["OTBarHits"]

iosvc = IOSvc()
iosvc.Input = "/afs/cern.ch/user/k/kenakkal/k4Reco/k4Reco/Simulation/alfaTrackerSimulation_Single_MuonMinus_100GeV_MDI_ECAL_TRACKER.root"
iosvc.Output = "alfaTrackerDigi_Single_MuonMinus_100GeV_main_MDI_ECAL.root"

# inp.collections = [
#     "VertexBarrelCollection",
#     "EventHeader",
# ]

hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink")
root_hist_svc.FileName = "ddplanardigi_hist_Single_MuonMinus_100GeV_main_MDI_ECAL_TRACKER.root"

ApplicationMgr(TopAlg=[digi],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
               OutputLevel=INFO,
               )
