/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef K4RECO_PARAMETERS_H
#define K4RECO_PARAMETERS_H

#include <map>
#include <string>
#include <vector>
/* struct and class in c++ are nearly identical. The only difference being default memeber visibility. Struct : public by def, Class : private by def  
Using struct here is a style signal : this is mainly a bundle of data and not an object with tightly guarded internal state/ Every step in the python 
config becomes one of these. Parameters.h : per-step cut 
*/

struct Parameters {
  // using aliases 
  using ParMap = std::map<std::string, double>;
  using StringVec = std::vector<std::string>;

public:
/* 18 positional parameters, each coppied into a matching m_.. memeber via the initialiser list */
  Parameters(std::vector<int> const& collections, double maxCellAngle, double maxCellAngleRZ, double chi2cut,
             int minClustersOnTrack, double maxDistance, double maxSlopeZ, double highPTcut, bool highPTfit,
             bool onlyZSchi2cut, bool radialSearch, bool vertexToTracker, bool kalmanFitForward, int step, bool combine,
             bool build, bool extend, bool sortTracks)
      : m_collections(collections), m_maxCellAngle(maxCellAngle), m_maxCellAngleRZ(maxCellAngleRZ), m_chi2cut(chi2cut),
        m_minClustersOnTrack(minClustersOnTrack), m_maxDistance(maxDistance), m_maxSlopeZ(maxSlopeZ),
        m_highPTcut(highPTcut), m_highPTfit(highPTfit), m_onlyZSchi2cut(onlyZSchi2cut), m_radialSearch(radialSearch),
        m_vertexToTracker(vertexToTracker), m_kalmanFitForward(kalmanFitForward), m_step(step), m_combine(combine),
        m_build(build), m_extend(extend), m_sortTracks(sortTracks) {}

  Parameters(Parameters const&) = default; // copy constructor; creates a brand new Parameters object by copying every memeber from an existing one
  Parameters& operator=(Parameters const&) = delete; //copy assignment forbidden -> you cannot take an existing Parameter object and overwrite it with another one's values via =
  Parameters(Parameters&&) = default; // move constructor: building a new object  by stealing the contents of a temp/expiring one, chepaer than full copy; allowed
  Parameters& operator=(Parameters&&) = delete; // move assignment forbidden : same as copy assignment but from a temp/expiring object 
  ~Parameters() = default; // destructor 

  std::vector<int> m_collections; /// which collections to combine; which hits to use. vector of int indices, not hits themselves
  double m_maxCellAngle; // threshold for Cell:getAngle()
  double m_maxCellAngleRZ; // threshold for Cell:getAngleRZ()
  double m_chi2cut; // fit-quality threshold once the cells and chained and fitted
  int m_minClustersOnTrack; // min chain length (hits) to accept as a track 
  double m_maxDistance; // limits which hit pairs are even considered for a Cell in the first palce 
  double m_maxSlopeZ; /// Cut on the slope in the longitudinal plane for seeding
  double m_highPTcut; //pT/curvature cut above which a track is considered high-pt
  bool m_highPTfit;
  bool m_onlyZSchi2cut; // arc length vs z chi2 cut 
  bool m_radialSearch; // build/extend chains outwards by increasing radius (?)
  bool m_vertexToTracker; 
  bool m_kalmanFitForward;
  int m_step;
  bool m_combine;
  bool m_build;
  bool m_extend;
  bool m_sortTracks;
  double m_tightenStep = 1;

  const StringVec m_existingFunctions = {
      "CombineCollections",
      "ExtendTracks",
      "BuildNewTracks",
      "SortTracks",
  };
  const StringVec m_existingFlags = {
      "HighPTFit", "OnlyZSchi2cut", "RadialSearch", "VertexToTracker", "KalmanFitForward", "KalmanFitBackward",
  };
  const StringVec m_existingParameters = {
      "MaxCellAngle", "MaxCellAngleRZ", "Chi2Cut", "MinClustersOnTrack", "MaxDistance", "SlopeZRange", "HighPTCut",
  };

  /* tighten the algo. The 3 of the acceptanmce cuts gets multiplied by the factor, which is always < 1. 
  A smaller MaxCellAngle/MaxCellAngleRZ means fewer cell pairs are considered "compatible enough" to chain; 
  a smaller chi2cut means fewer fitted tracks pass the quality threshold. So calling this function makes the algorithm stricter */
  void tighten() {
    double factor = (10.0 - m_tightenStep) / (10.0 - (m_tightenStep - 1.0));
    m_maxCellAngle *= factor;
    m_maxCellAngleRZ *= factor;
    m_chi2cut *= factor;
    m_tightenStep += 1.0;
  }

private:
  void check(StringVec const& values, StringVec const& options, std::string const& type);
  void check(ParMap const& values, StringVec const& options, std::string const& type);
};

#endif // K4RECO_PARAMETERS_H
