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
from Gaudi.Configuration import INFO, DEBUG
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc, DDPlanarDigi, GeoSvc, UniqueIDGenSvc,  RootHistSvc, Gaudi__Histograming__Sink__Root as RootHistoSink
import os

id_service = UniqueIDGenSvc("UniqueIDGenSvc")

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["/afs/cern.ch/user/k/kenakkal/k4geo/FCCee/ALFA/compact/ALFA_o1_v00/ALFA_o1_v00.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

digi = DDPlanarDigi()
digi.SubDetectorName = "OTBar"    #"SiWrB"
digi.IsStrip = False
#smearing resolution in mm applied along the two local directions on eachs ensor plane. 
# Here 6 values corresponds to the 6 vertex detetcors layers. will need to modify this for alfa 
digi.ResolutionU = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
digi.ResolutionV = [0.003, 0.003, 0.003, 0.003, 0.003, 0.003]
digi.SimTrackHitCollectionName = ["OTBarCollection"] #input collection name; change it according to the root file 
#output container storing the relation btw each digi hit and its corresponsding sim hit so as to trace back which MC hit produced the digi hit 
digi.SimTrkHitRelCollection = ["OTBarHitRelations"]
digi.TrackerHitCollectionName = ["OTBarHits"] # conatiner for new digi hits after smearing 

iosvc = IOSvc()
iosvc.Input = "alfaTrackerSimulation.root"
iosvc.Output = "alfaTrackerDigi.root"

# inp.collections = [
#     "VertexBarrelCollection",
#     "EventHeader",
# ]

#Histogramming services
hps = RootHistSvc("HistogramPersistencySvc")
root_hist_svc = RootHistoSink("RootHistoSink") #writes histos diretcly to a root file : ddplanardigi_hist.root
root_hist_svc.FileName = "ddplanardigi_hist.root"

#TopAlg : algo executed once per event
ApplicationMgr(TopAlg=[digi],
               EvtSel="NONE",
               EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc"), root_hist_svc],
               OutputLevel=INFO,
               )
