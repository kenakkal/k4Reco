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
#ifndef K4RECO_CELL_H
#define K4RECO_CELL_H

#include "KDCluster.h"

#include <cmath>
#include <memory>
#include <vector>

#ifdef CF_USE_VDT
#include <vdt/atan.h>
#endif

// ------------------------------------------------------------------------------------
// The Cell class is a simple object which connects two points in 2D space. It is
// used in Cellular Automaton tracking to create tracks, by connecting all plausable
// hit points with cells, and linking these cells together in a chain. Each cell
// needs to know what cells it is connected to, and hold a weight determined by
// its position in the chain.
// ------------------------------------------------------------------------------------

class Cell {
public:
  //Type aliases. cell to cell links uses weak_ptr
  typedef std::vector<std::weak_ptr<Cell>> WeakCells;
  typedef std::shared_ptr<Cell> SCell;
  typedef std::weak_ptr<Cell> WCell;

public:
  // Constructors, main initialisation is with two kd hits
  Cell() { m_weight = 0; }
 // copy constructs forbidden
  Cell(const Cell&) = delete;
  //copy asignment forbidden
  Cell& operator=(const Cell&) = delete;
  // move constructor defaulted
  Cell(Cell&&) = default;
  //move assignment constructor 
  Cell& operator=(Cell&&) = default;
  ~Cell() = default;

  // m_gradient = delV/delU : slope of a line ( v = mu +c) which connetcs these 2 hits in the conformal space 
  // m_gradientRZ = delR/delZ; re
  Cell(SKDCluster const& hit1, SKDCluster const& hit2)
      : m_weight(0), m_gradient((hit2->getV() - hit1->getV()) / (hit2->getU() - hit1->getU())),
        m_gradientRZ((hit2->getRadius() - hit1->getRadius()) / (hit2->getZ() - hit1->getZ())), m_start(hit1),
        m_end(hit2) {}

  // Weight of the cell (first cell in a chain has weight 0, and each subsequent link has weight +1)
  int getWeight() const { return m_weight; }
  void setWeight(int weight) { m_weight = weight; }

  // Gradient of the cell connecting two hits
  double getGradient() const { return m_gradient; }
  void setGradient(double gradient) { m_gradient = gradient; }
  double getGradientRZ() const { return m_gradientRZ; }

  // Angle between two cells. This is assumed to be less than 90 degrees. The below code block shows function overloading: getAngle() & getAngleRZ()have two versions 
  /* SCell const& cell2 : parameter is a std::shared_ptr<Cell>, taken by const reference; 
  cell2.get(): every shared_ptr has .get() method that returns a raw pointer (a plain Cell*), without affecting the ownership
  *(cell2.get()): dereferncing the raw pointer with * giving the Cell object 
  getAngle(*(cell2.get()) : now that we have an actual Cell (not a shared_ptr or raw_ptr), this calss the other overload - 
  getAngle(Cell const& cell2) the calculates the angle btw two cells using tangent substraction math 
  */
  inline double getAngle(SCell const& cell2) const { return getAngle(*(cell2.get())); }
  inline double getAngleRZ(SCell const& cell2) const { return getAngleRZ(*(cell2.get())); }

  /*tangent substraction math tan(theta2 - theta1) = tan theta1 - tan theta 2/ 1 + tan theta1 * tan theta 2*/
  inline double getAngle(Cell const& cell2) const {
#ifdef CF_USE_VDT
    return fabs(vdt::fast_atan((cell2.m_gradient - m_gradient) / (1 + m_gradient * cell2.m_gradient)));
#else
    return fabs(std::atan((cell2.m_gradient - m_gradient) / (1 + m_gradient * cell2.m_gradient)));
#endif
  }

  inline double getAngleRZ(Cell const& cell2) const {
#ifdef CF_USE_VDT
    return fabs(vdt::fast_atan((cell2.m_gradientRZ - m_gradientRZ) / (1 + m_gradientRZ * cell2.m_gradientRZ)));
#else
    return fabs(std::atan((cell2.m_gradientRZ - m_gradientRZ) / (1 + m_gradientRZ * cell2.m_gradientRZ)));
#endif
  }

  // Start and end points of the cell
  SKDCluster const& getStart() const { return m_start; }
  SKDCluster const& getEnd() const { return m_end; }

  // Increment the cell weight (usually if the chain length is extended upstream of this cell)
  void update(SCell const& cell2) {
    if ((cell2->getWeight() + 1) > m_weight)
      m_weight = cell2->getWeight() + 1;
  }

  // The cell has a memory of all cells that connect to it, and all cells that it connects to. If several cells point to
  // this cell, then the weight taken from the highest weighted of those (longest chain)
  /* setFrom() does everything update() does, plus it recors cell2 as predecessor in m_from and stores cell2's weight in m_weights */
  
  /* m_from & m_to are weaak_ptr. if m_from and m_to used shared_ptr, you'd get a reference cycle — cell A holds a 
  shared_ptr to cell B (via setTo), and cell B holds a shared_ptr back to cell A (via setFrom). Neither's reference 
  count would ever drop to zero, even after nothing else references them, so they'd never get destroyed — a classic 
  shared_ptr memory leak. weak_ptr breaks the cycle: it observes an object without owning it (doesn't increment the 
  reference count), so cells can point at each other freely without keeping each other alive artificially. The tradeoff 
  is that code has to call .lock() on a weak_ptr to actually use it (which returns a temporary shared_ptr, or a null 
  one if the object's already gone)*/

  void setFrom(SCell const& cell2) {
    m_from.push_back(WCell(cell2)); // m_from : "which cells lead into me?" predecessors, used to trace the chain backwards  once you found the highest-weight endpoint  
    m_weights.push_back(cell2->getWeight());
    if ((cell2->getWeight() + 1) > m_weight)
      m_weight = cell2->getWeight() + 1;
  }
  void setTo(SCell const& cell2) { m_to.push_back(WCell(cell2)); } // m_to = "which cells I lead into"; successors, used when building/extending chains forward
  WeakCells& getFrom() { return m_from; }
  WeakCells& getTo() { return m_to; }

  // doca : distance of closest approach 
  double doca() const {
    // solving for c in v = m.u + c using start hit's own (u,v) 
    double intercept = m_start->getV() - m_start->getU() * m_gradient;
    // point to line distance formula for a point fromn the origin to line v = mu + c; d = |c|/sqrt(m2 + 1)
    double doca = fabs(intercept) / sqrt(m_gradient * m_gradient + 1.);
    return doca;
  }

private:
  // Each cell contains a weight, a gradient, two hits which it connects
  // and a list of cells that it connects to or from
  int m_weight = 0;
  double m_gradient = 0.0;
  double m_gradientRZ = 0.0;
  SKDCluster m_start = nullptr;
  SKDCluster m_end = nullptr;
  WeakCells m_from{};
  std::vector<int> m_weights{};
  WeakCells m_to{};
};

//aliases
using SCell = Cell::SCell; // Cell::SCell  : typedef std::shared_ptr<Cell> SCell; 
/* vector of cells; . If cell A points to cell B points to cell C (via the setTo/setFrom), a cellularTrack is 
literally {A, B, C} — the ordered list of cells making up one candidate track path, before it's been fitted into a KDTrack. */
using cellularTrack = std::vector<SCell>; 
using SharedCells = std::vector<std::shared_ptr<Cell>>; // std::vector<SCell> 
/* Unlike shared_ptr (which allows many owners, reference-counted), unique_ptr allows exactly one owner —
 no reference counting overhead at all, and it literally cannot be copied (only moved), enforced by the compiler. 
 UcellularTrack is a single candidate chain, wrapped so that whoever holds it is the sole owner of that particular 
chain object.*/
using UcellularTrack = std::unique_ptr<cellularTrack>;
/* The full output of a pattern-recognition pass: a list of candidate track chains, each one uniquely owned.*/
using UniqueCellularTracks = std::vector<std::unique_ptr<cellularTrack>>;
#endif
