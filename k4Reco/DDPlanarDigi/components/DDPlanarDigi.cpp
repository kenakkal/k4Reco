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
#include "DDPlanarDigi.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "Gaudi/Accumulators/RootHistogram.h"

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "TMath.h"

#include <cmath>
#include <fmt/format.h>

DDPlanarDigi::DDPlanarDigi(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("SimTrackHitCollectionName", {"SiWrBCollection"}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("TrackerHitCollectionName", {"SiWrBHits"}),
                        KeyValues("SimTrkHitRelCollection", {"SiWrBHitRelations"})}) {
  
  std::cout << "######## USING MY LOCAL DDPlanarDigi ########" << std::endl;
  /* service<T>(name,createIfNotExisting) : Gaudi template helper available to any Gaudi component via inheritance taht looks up a service by name 
  and interface type and returns back a SmartIF<T>*/
  m_uidSvc = service<IUniqueIDGenSvc>("UniqueIDGenSvc", true);
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
  }

  if (m_resULayer.size() != m_resVLayer.size()) {
    error() << "DDPlanarDigi - Inconsistent number of resolutions given for U and V coordinate: "
            << "ResolutionU  :" << m_resULayer.size() << " != ResolutionV : " << m_resVLayer.size();

    throw std::runtime_error("DDPlanarDigi: Inconsistent number of resolutions given for U and V coordinate");
  }

  m_histograms[hu].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{this, "hu", "smearing u", {50, -5., +5.}});
  m_histograms[hv].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{this, "hv", "smearing v", {50, -5., +5.}});
  m_histograms[hT].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{this, "hT", "smearing time", {50, -5., +5.}});

  m_histograms[diffu].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{this, "diffu", "diff u", {1000, -.1, +.1}});
  m_histograms[diffv].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{this, "diffv", "diff v", {1000, -.1, +.1}});
  m_histograms[diffT].reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "diffT", "diff time", {1000, -5., +5.}});

  m_histograms[hitE].reset(
      new Gaudi::Accumulators::StaticRootHistogram<1>{this, "hitE", "hitEnergy in keV", {1000, 0, 200}});
  m_histograms[hitsAccepted].reset(new Gaudi::Accumulators::StaticRootHistogram<1>{
      this, "hitsAccepted", "Fraction of accepted hits [%]", {201, 0, 100.5}});
}

StatusCode DDPlanarDigi::initialize() {
  info() << "######### within initialize()##########" << endmsg;
  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  if (!m_geoSvc) {
    error() << "Unable to retrieve the GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  const auto detector = m_geoSvc->getDetector();
  
  // extension<T>() : type-based plugin/extension retrieval mechanism.
  const auto surfMan = detector->extension<dd4hep::rec::SurfaceManager>();
  dd4hep::DetElement det = detector->detector(m_subDetName.value());
  surfaceMap = surfMan->map(m_subDetName.value());

  if (!surfaceMap) {
    throw std::runtime_error(fmt::format("Could not find surface map for detector: {} in SurfaceManager", det.name()));
  }

  // Get and store the name for a debug message later
  (void)this->getProperty("SimTrackHitCollectionName", m_collName);

  //builds a bitmask used later to truncate cellIDs
  /* this loop is excecuted iff the user has specified a value ffor m_cellIDBits that is not the default value of 64.
   << is the left shift operator. slides every bit in the number to the left by however many positions you specify, filling zeros from the right
   for eg> 1 in 8bit binary is 00000001, if m_cellIDBits is 4;  1<< 4 = 00010000 (16 in decimnal 1 x2^4); 
   now (1<<4) -1 = 00001111 (15 in decimal) : boorrows from the single set bit, flipping it to 0 and turning all the zeros below it to 1   */
  if (m_cellIDBits != 64) { // m_cellIDBits specifies how many low-order bits of the 64-bit cellID carry meaningful info 
    m_mask = (static_cast<std::uint64_t>(1) << m_cellIDBits) - 1; // m_mask defaults everything to 1; by default everything passes through unchanged
  }

  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
DDPlanarDigi::operator()(const edm4hep::SimTrackerHitCollection& simTrackerHits,
                         const edm4hep::EventHeaderCollection& headers) const {
  auto seed = m_uidSvc->getUniqueID(headers[0].getEventNumber(), headers[0].getRunNumber(), this->name());
  debug() << "Using seed " << seed << " for event " << headers[0].getEventNumber() << " and run "
          << headers[0].getRunNumber() << endmsg;
  auto rng_engine = TRandom2(seed);

  int nCreatedHits = 0;
  int nDismissedHits = 0;

  auto trkhitVec = edm4hep::TrackerHitPlaneCollection();
  auto thsthcol = edm4hep::TrackerHitSimTrackerHitLinkCollection();

  /* m_encodingStringVariable.value () : unwarps the property to get std::sting coz m_encodingStringVariable is a Gaudi::Property<std::string> type; default string is "GlobalTrackerReadoutID" (check the header file). 
   ->constantAsString(name) : method on geometry service interface that looks up a dd4hep xml constant by name and returns its value as string */
  std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringVariable.value());
  // dd4hep::DDSegmentation::BitFieldCoder : class that parses comma-seperated string at the construction time, building an internal table of field names. 
  // Once built, this bitFieldCoder object can then decode any raw 64-bit cellID integer into its individual named components on demand
  dd4hep::DDSegmentation::BitFieldCoder bitFieldCoder(cellIDEncodingString);

  int nSimHits = simTrackerHits.size();
  debug() << "Processing collection " << m_collName << " with " << simTrackerHits.size() << " hits ... " << endmsg;

  for (const auto& hit : simTrackerHits) {
    ++(*m_histograms[hitE])[hit.getEDep() * (dd4hep::GeV / dd4hep::keV)];

    if (hit.getEDep() < m_minEnergy) {
      debug() << "Hit with insufficient energy " << hit.getEDep() * (dd4hep::GeV / dd4hep::keV) << " keV" << endmsg;
      continue;
    }

    // hit.cellID() : returns raw cellID, a 64-bit integer that has multiple peices of the geo info packed iinto different bit ranges based on the encoding string
    
    /* this & here is bit-wise AND operator.  1 & 1 = 1, every other combo is 0; 
    AND-ing against a mask of N low 1-bits keeps the low N bits of the cellID exactly as they were, and forcibly zeroes out everything above bit N.
    This matters coz the deetectors actual cellID encoding may only occupy a certain nuumber of bits total. If the raw 64-bit integer comming out of 
    the simulation has any stray bits set above teh meaniful range of the encoding, two cellIDs that are supposed to represnet the sam ephysical sensor 
    could comapre as different. Trunctting with this mask gaurentess that only teh bits that are actually used for the cellID are kept, and any stray bits above that are zeroed out.  */   
    
    const std::uint64_t cellID = hit.getCellID() & m_mask;

    // get the measurement surface for this hit using the CellID
    // .find() doesnt return teh surfuce directly. It gives an interator to the map entry, which is a pair of (key, value) where key is the cellID and value : surfuce pointer 
    dd4hep::rec::SurfaceMap::const_iterator sI = surfaceMap->find(cellID);

    if (sI == surfaceMap->end()) {
      throw std::runtime_error(fmt::format("DDPlanarDigi::processEvent(): no surface found for cellID : {}", cellID));
    }

    const dd4hep::rec::ISurface* surf = sI->second;
    int layer = bitFieldCoder.get(cellID, "layer");

    // true, unsmeared 3D position (in global detector cordinates, mm - EDM4hep std convention) where the particle actually crossed the sensor accoding to G4 simulation 
    dd4hep::rec::Vector3D oldPos(hit.getPosition()[0], hit.getPosition()[1], hit.getPosition()[2]); 
    dd4hep::rec::Vector3D newPos;

    //  Check if Hit is inside sensitive
    /*convert the true hit position into dd4hep's internal units abd check whether it falls within the senosr's actual physical bounds. If it does not,
    the hit sits outside trhe senosr it is nominally assigned t0, then enter this block */
    if (!surf->insideBounds(dd4hep::mm * oldPos)) {
      debug() << "  hit at " << oldPos
               << " " << cellID
               << " is not on surface "
               << *surf
               << " distance: " << surf->distance(  dd4hep::mm * oldPos )
               << endmsg;
      //debug() << "hit at " << oldPos << " " << " is not on the surfuce " << *surf  << endmsg;  
      if (m_forceHitsOntoSurface) {
        debug() << "forcing the hit onto the surface " << endmsg; 
        dd4hep::rec::Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos); // mathematically project the hit position straight onto the sensor's flat surfuce; thereby removing any small perpendicular offsetr was causing it to be flagged outside bounds
        dd4hep::rec::Vector3D oldPosOnSurf = (1. / dd4hep::mm) * surf->localToGlobal(lv); // oldPosOnSurf similar to oldPos should be in plain mm and not in dd4hep::mm; hence the divide 

        debug() << " moved to " << oldPosOnSurf << " distance " << (oldPosOnSurf - oldPos).r() << endmsg;

        oldPos = oldPosOnSurf;

      } else {
        ++nDismissedHits;
        debug() << " hit dismissed coz m_foreceHitsOntoSurfuce is disabled " << endmsg;
        continue; // loop control statement; skip the rest of this loop iteration and move on to the next hit in the collection; no digi hit will be created for this sim hit
      }
    }
    debug() << " hit is on the surfuce " << endmsg;
    // Smear time of the hit and apply the time window cut if needed
    double hitT = hit.getTime(); // simulated time at which this energy deposit ocurred in the sensor (ns: EDM4hep std convention)

    if (m_resTLayer.size() && m_resTLayer[0] > 0) {
      float resT = m_resTLayer.size() > 1 ? m_resTLayer[layer] : m_resTLayer[0];

      // rng_engine.Guas(mean, sigma): draws a random sample from a Gaussain distribution (here) centered at 0 with std. dev: resT; random timing jitters in a real detector 
      double tSmear = resT > 0 ? rng_engine.Gaus(0, resT) : 0;
      /* pull/normalised residual : tSmear/resT; dividing the smear offset by the resolution that generated it. If you draw many samples from a 
      Gaussian with mean 0 and standard deviation resT, and then divide every single sample by that same resT, the result is mathematically guaranteed 
      to follow a standard normal distribution — mean 0, standard deviation exactly 1 — regardless of what resT's actual numeric value was. So if you 
      fill this histogram over thousands of events and it comes out looking like a clean, unit-width bell curve, that's a strong confirmation the 
      smearing math is implemented correctly. If it comes out systematically wider or narrower than a unit Gaussian, something's off 
      (e.g., a units bug, or the histogram binning not actually matching what's being filled).
      */
      ++(*m_histograms[hT])[resT > 0 ? tSmear / resT : 0];
      // raw smearing histo that records the un-normalised smearing offset in ns. Shows the smearing magntude in genuibne physical units
      ++(*m_histograms[diffT])[tSmear];

      hitT += tSmear;
      debug() << "smeared hit at T: " << hit.getTime() << " ns to T: " << hitT
              << " ns according to resolution: " << resT << " ns" << endmsg;
    }

    // Correcting for the propagation time
    /*This correction subtracts out the geometric travel-time delay; so the corrected hitT reflects timing info 
    closer to "when did something interesting happen," rather than "when did the signal happen to physically arrive here." */
    if (m_correctTimesForPropagation) {
      debug() << "correcting for the propagation time ..." << endmsg; 
      // travel time  = distance/speed; speed = speed of light (c); 1e6 : conversion from m/s to mm/ns 
      double dt = oldPos.r() / (TMath::C() / 1e6);
      hitT -= dt;
      debug() << "corrected hit at R: " << oldPos.r() << " mm by propagation time: " << dt << " ns to T: " << hitT
              << " ns" << endmsg;
    }

    // Skip the hit if its time is outside the acceptance time window
    if (m_useTimeWindow) {
      // TODO: Check that the length of the time window is OK
      float timeWindow_min = m_timeWindowMin.size() > 1 ? m_timeWindowMin[layer] : m_timeWindowMin[0];
      float timeWindow_max = m_timeWindowMax.size() > 1 ? m_timeWindowMax[layer] : m_timeWindowMax[0];
      if (hitT < timeWindow_min || hitT > timeWindow_max) {
        debug() << "hit at T: " << hit.getTime() << " smeared to: " << hitT
                << " is outside the time window: hit dropped" << endmsg;
        ++nDismissedHits;
        continue;
      }
    }

    // Try to smear the hit position but ensure the hit is inside the sensitive region
    /* .u(),.v() : returns the sensor's in-plane basis vectors expressed in 3D global cordinates.   
    if i take one step of length 1 along this sensor's own u direction, which way does that correspond to in the detetcor's overall global x/y/z frame */
    dd4hep::rec::Vector3D u = surf->u();
    dd4hep::rec::Vector3D v = surf->v();

    // Get local coordinates on surface
    dd4hep::rec::Vector2D lv = surf->globalToLocal(dd4hep::mm * oldPos);

    // sim hit's local position on the senosr in mm 
    /* converting to local cordinates for smearing than directly applying the smear to the global x/y/z is necessary coz spatial resolution is a ppty of a detector's 
    local geometry. The resolution characteristics are expressed in sensor's own local frame. If you wanna apply smearing to the global cordinates, you have to consatntly account
    for how each indiuvidual senosr are rotated in 3d space. Going into local cordinates sidesteps taht complexity entirely and conversion to global cordiunates later will handle 
    the rotation responsiblity internally by dd4hep than implementing it with hand.  */
    double uL = lv[0] / dd4hep::mm; // divide by the dd4hep unit to convert back to plain mm from dd4hep::mm
    double vL = lv[1] / dd4hep::mm;

    bool acceptHit = false;
    int tries = 0;

    // TODO: check lengths
    float resU = m_resULayer.size() > 1 ? m_resULayer[layer] : m_resULayer[0];
    float resV = m_resVLayer.size() > 1 ? m_resVLayer[layer] : m_resVLayer[0];

    while (tries < m_maxTries) {
      // if( tries > 0 ) debug() << "retry smearing for " <<  cellid_decoder( hit ).valueString() << " : retries " <<
      // tries << endmsg;

      double uSmear = rng_engine.Gaus(0, resU);
      double vSmear = rng_engine.Gaus(0, resV);

      dd4hep::rec::Vector3D newPosTmp;
      if (m_isStrip) {
        if (m_subDetName == "SET") {
          double xStripPos, yStripPos, zStripPos;
          // Find intersection of the strip with the z=centerOfSensor plane to set it as the center of the SET strip
          dd4hep::rec::Vector3D simHitPosSmeared =
              (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
          zStripPos = surf->origin()[2] / dd4hep::mm;
          double lineParam = (zStripPos - simHitPosSmeared[2]) / v[2];
          xStripPos = simHitPosSmeared[0] + lineParam * v[0];
          yStripPos = simHitPosSmeared[1] + lineParam * v[1];
          newPosTmp = dd4hep::rec::Vector3D(xStripPos, yStripPos, zStripPos);
        } else {
          newPosTmp = (1. / dd4hep::mm) * (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, 0.)));
        }
      } else {
        newPosTmp =
            (1. / dd4hep::mm) *
            (surf->localToGlobal(dd4hep::rec::Vector2D((uL + uSmear) * dd4hep::mm, (vL + vSmear) * dd4hep::mm))); // newPosTmp unit : plain mm. localToGloabl() : results in dd4hep internal units; divison to get in plain mm 
      }

      debug() << " hit at    : " << oldPos << " smeared to: " << newPosTmp << " uL: " << uL << " vL: " << vL
              << " uSmear: " << uSmear << " vSmear: " << vSmear << endmsg;

      if (surf->insideBounds(dd4hep::mm * newPosTmp)) {
        acceptHit = true;
        newPos = newPosTmp;

        ++(*m_histograms[hu])[uSmear / resU];
        ++(*m_histograms[hv])[vSmear / resV];

        ++(*m_histograms[diffu])[uSmear];
        ++(*m_histograms[diffv])[vSmear];

        break; // exits the while loop 
      }

      // debug() << "  hit at " << newPosTmp
      //         << " " << cellid_decoder( hit).valueString()
      //         << " is not on surface "
      //         << " distance: " << surf->distance( dd4hep::mm * newPosTmp )
      //         << endmsg;

      ++tries;
    }

    if (!acceptHit) {
      debug() << "hit could not be smeared within ladder after " << m_maxTries << "  tries: hit dropped" << endmsg;
      ++nDismissedHits;
      continue; // moves to the next hit skippinh the histo fill, output creation etc 
    }

    auto trkHit = trkhitVec.create();

    trkHit.setCellID(cellID);

    trkHit.setPosition(newPos.const_array());
    trkHit.setTime(hitT);
    trkHit.setEDep(hit.getEDep());
    // Computing the orientation angles
    /* Storing theta/phi is a compact way to record "which way did this sensor's local axes point, 
    in the detector's overall frame"*/
    float u_direction[2];
    u_direction[0] = u.theta();
    u_direction[1] = u.phi();

    float v_direction[2];
    v_direction[0] = v.theta();
    v_direction[1] = v.phi();

    debug() << " U[0] = " << u_direction[0] << " U[1] = " << u_direction[1] << " V[0] = " << v_direction[0]
            << " V[1] = " << v_direction[1] << endmsg;


    // Storing the orienattion and resolution of the output hits 
    trkHit.setU(u_direction);
    trkHit.setV(v_direction);
    trkHit.setDu(resU);

    if (m_isStrip) {
      // store the resolution from the length of the wafer - in case a fitter might want to treat this as 2d hit ....
      double stripRes = surf->length_along_v() / dd4hep::mm / std::sqrt(12);
      trkHit.setDv(stripRes);
      // TODO: Set type?
      // trkHit.setType( UTIL::set_bit( trkHit.getType() ,  UTIL::ILDTrkHitTypeBit::ONE_DIMENSIONAL ) );

    } else {
      trkHit.setDv(resV);
    }

    auto association = thsthcol.create();
    association.setTo(hit);// truth hits/ sim hits
    association.setFrom(trkHit); // digi hits 

    ++nCreatedHits;
  }

  // Filling the fraction of accepted hits in the event
  float accFraction = nSimHits > 0 ? float(nCreatedHits) / float(nSimHits) * 100.0 : 0.0;
  ++(*m_histograms[hitsAccepted])[accFraction];

  debug() << "Created " << nCreatedHits << " hits, " << nDismissedHits << " hits  dismissed" << endmsg;

  return std::make_tuple(std::move(trkhitVec), std::move(thsthcol));
}
