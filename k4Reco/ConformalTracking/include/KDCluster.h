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
#ifndef K4RECO_KDCLUSTERS_H
#define K4RECO_KDCLUSTERS_H

#include <edm4hep/TrackerHitPlane.h>

#include <cmath>
#include <memory>
#include <vector>

// ------------------------------------------------------------------------------------
// The KDCluster class is a simple hit class used in the Cellular Automaton tracking.
// They are designed to be a lightweight object containing all of the information
// needed for tracking in conformal space, to be ordered in a binary tree (KDtree,
// hence the name). The hits contain their co-ordinates in conformal space, both in
// cartesian (u-v) and polar (r-theta) notation. Cartesian positions allow fast nearest
// neighbour searches, while the polar positions allow both nearest neighbour searches
// in theta (avoiding sector definitions) and directed tracking to flow inside-out or
// outside-in. They additionally contain minimal detector information (id, layer and side)
// ------------------------------------------------------------------------------------

/*KDCluster is a tracking hit class. Every TrackerHitPlane object in ROOT file gets converted into one of this before any 
pattern recogonition happens. It is deliberately lighweight (just numbers, no heavy geometry objects) as millions of these 
get created and compared during the search*/

class KDCluster {
public:
  // Constructors, main initialisation is with tracker hit
  KDCluster() = default;
  
  /* constructor initialiser list :takes in the raw hits (x,y,z) from edm4hep:TrackerHitPlane and 
  copies them unchanged and is told whether this hit came from endcap or not */
  
  KDCluster(const edm4hep::TrackerHitPlane& hit, bool endcap, bool forward = false)
      : m_x(hit.getPosition()[0]), m_y(hit.getPosition()[1]),
        m_z(hit.getPosition()[2]), // Store the (unaltered) z position
        m_endcap(endcap), m_forward(forward) {
    // Calculate conformal position in cartesian co-ordinates
    const double radius2 = (m_x * m_x + m_y * m_y);
    const double radius2Inv = 1. / radius2;
    const double radius = sqrt(radius2);
    m_u = m_x * radius2Inv;
    m_v = m_y * radius2Inv;
    // Note the position in polar co-ordinates
    m_r = 1. / radius; /* confirmal polar radius; distance of (u,v) from origin in conformal space. m_r is not a real-space distance.   
    How far out the point has moved in the conformal space. note the inversion -> hits close to the beam (smaller radius) end up far from the origin in uv space and vice versa 
    */
    m_theta = atan2(m_v, m_u) + M_PI; // polar angle shifted by Pi to stay in the covenient positive range; angular position of the hit after the conformal mapping 
    m_radius = radius; // cartesian xy radius
    // Get the error in the conformal (uv) plane
    // This is the xy error projected. Unfortunately, the
    // dU is not always aligned with the xy plane, it might
    // be dV. Check and take the smallest
    
    /* this code block picks up the smaller of the sensor's two local measurement errors (Du, Dv)to represent 
    the transverse error and the larger one becomes teh z-direction error. This is done coz the senor's local u/v 
    is not guarenteed to align with global x-y/z . Rough way of guessing which local error coeesponds to which 
    global direction*/

    if (hit.getDv() < m_error) {
      m_error = hit.getDv(); // sensor's position uncertainity in the tranverse direction; picked to be the smaller value of the sesnsor's local measurement errors
      m_errorZ = hit.getDu();
    } else {
      m_error = hit.getDu();
      m_errorZ = hit.getDv();
    }
    // converts real-space error into conformal space error (m_errorU, m_errorV)
    const double sinTheta = sin(m_theta);
    const double cosTheta = cos(m_theta);
    if (endcap) {
      m_errorU = (m_error * std::abs(sinTheta) + m_errorZ * std::abs(cosTheta)) * (m_r * m_r);
      m_errorV = (m_error * std::abs(cosTheta) + m_errorZ * std::abs(sinTheta)) * (m_r * m_r);
      m_errorX = m_error * sinTheta;
      m_errorY = m_error * cosTheta;

      // Need to set endcap error in z!
      m_errorZ = 0.25;

    } else {
      //m_errorU = m_errorX * m_r^2 (u = x/r^2; m_r = 1/r)
      m_errorU = m_error * m_r * m_r * sinTheta;
      m_errorV = m_error * m_r * m_r * cosTheta;

      // m_error is pointing perpendicular to the hit. the two code lines below are the errors on global X,Y 
      m_errorX = m_error * sinTheta;
      m_errorY = m_error * cosTheta;
    }
  }
  // copy constructor, = delete -> forbids it entirely 
  KDCluster(const KDCluster&) = delete;
  // copy assignment forbidden 
  KDCluster& operator=(const KDCluster&) = delete;
  // move constructor and move assignment constructor; = default tells teh compiler to generate the normal version; move does not copy the data, it simply transfers the ownership from a temp/expiring object into a new one, leavaing old one empty 
  KDCluster(KDCluster&&) = default;
  KDCluster& operator=(KDCluster&&) = default;
  // destructor function 
  ~KDCluster() = default;

  double getX() const { return m_x; }
  double getY() const { return m_y; }
  double getU() const { return m_u; }
  double getV() const { return m_v; }
  double getR() const { return m_r; }
  double getRadius() const { return m_radius; } // "real" radius in xy
  double getTheta() const { return m_theta; }
  double getZ() const { return m_z; }
  float getError() const { return m_error; }
  float getErrorX() const { return m_errorX; }
  float getErrorY() const { return m_errorY; }
  float getErrorU() const { return m_errorU; }
  float getErrorV() const { return m_errorV; }
  float getErrorZ() const { return m_errorZ; }
  bool used() const { return m_used; }

  void setU(double u) { m_u = u; }
  void setV(double v) { m_v = v; }
  void setR(double r) { m_r = r; }
  void setTheta(double theta) { m_theta = theta; }
  void setZ(double z) { m_z = z; }
  void setError(float error) { m_error = error; }
  void used(bool used) { m_used = used; }

  // Subdetector information
  void setDetectorInfo(long subdet, long side, long layer, long module, long sensor) {
    m_subdet = subdet;
    m_side = side;
    m_layer = layer;
    m_module = module;
    m_sensor = sensor;
  }
  // Used long because this is what's used in DD4hep
  long getSubdetector() const { return m_subdet; }
  long getSide() const { return m_side; }
  long getLayer() const { return m_layer; }
  long getModule() const { return m_module; }
  long getSensor() const { return m_sensor; }
  bool forward() const { return m_forward; }
  bool endcap() const { return m_endcap; }

  // Check if another hit is on the same detecting layer
  // compares other hit  with this hit's own fields
  /* a real particle passing through a barrel detector normally leaves at most one hit per physical layer. 
  So if two KDClusters report the same subdetector+side+layer, they almost certainly can't both be legitimate 
  consecutive hits on the same track. This function is a cheap sanity filter used when the algorithm is deciding
  whether two hits are allowed to be linked into a Cell — reject same-layer pairings early, before doing any of 
  the more expensive u-v/R-z geometric checks.*/
  bool sameLayer(const std::shared_ptr<KDCluster> kdhit) const {
    if (kdhit->getSubdetector() == m_subdet && kdhit->getSide() == m_side && kdhit->getLayer() == m_layer)
      return true;
    return false;
  }

  // Check if another hit is on the same sensor of the same detecting layer
  /* Each == produces a bool (true/false) on its own. Chaining them with && produces one more bool — true only if 
  every comparison is true. Since the function's return type is already bool, you can just return that expression directly, 
  instead of writing an if condition */
  bool sameSensor(const std::shared_ptr<KDCluster> kdhit) const {
    return kdhit->getLayer() == m_layer && kdhit->getSubdetector() == m_subdet && kdhit->getSide() == m_side &&
           kdhit->getModule() == m_module && kdhit->getSensor() == m_sensor;
  }

private:
  // Each hit contains the conformal co-ordinates in cartesian
  // and polar notation, plus the subdetector information
  double m_x = 0.0;
  double m_y = 0.0;
  double m_u = 0.0;
  double m_v = 0.0;
  double m_r = 0.0;
  double m_radius = 0.0; // "real" radius in xy
  double m_z = 0.0;
  float m_error = 0.0;
  float m_errorX = 0.0;
  float m_errorY = 0.0;
  float m_errorU = 0.0;
  float m_errorV = 0.0;
  float m_errorZ = 0.0;
  double m_theta = 0.0;
  long m_subdet = 0;
  long m_side = 0;
  long m_layer = 0;
  long m_module = 0;
  long m_sensor = 0;
  bool m_used = false;
  bool m_endcap = false;
  bool m_forward = false;
};

/* typedef creates an alias - new name for an existing type. This line doesn't create anything new at runtime; 
it just tells the compiler "from now on, whenever you see SKDCluster in this codebase, treat it as meaning exactly 
std::shared_ptr<KDCluster>." So SKDCluster and std::shared_ptr<KDCluster> are 100% interchangeable — same type, just 
a shorter name.*/
typedef std::shared_ptr<KDCluster> SKDCluster;

/*SharedKDClusters means std::vector<std::shared_ptr<KDCluster>> — a resizable vector of shared-pointers-to-hits.*/
typedef std::vector<SKDCluster> SharedKDClusters;

#endif
