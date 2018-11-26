// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


/// \file  RegularSpline1D.cxx
/// \brief Implementation of RegularSpline1D class
///
/// modified by Felix Lapp 06.09.2018
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>


#include "RegularSpline1D.h"
#include <math.h>
#include <vector>
#include <iostream>

namespace ali_tpc_common {
namespace tpc_fast_transformation {


RegularSpline1D::RegularSpline1D()
  :
  mNumberOfKnots(0) //Amount of Knots needs to be 0
{
  /// Default constructor. Creates an empty uninitialised object
}

RegularSpline1D::~RegularSpline1D()
{
  mNumberOfKnots = 0;
}   


void RegularSpline1D::construct( int numberOfKnots )  
{
  /// Constructor.
  /// Initialises the spline with a grid with numberOfKnots knots in the interval [0,1]
  /// array inputKnots[] has numberOfKnots entries, ordered from 0. to 1.
  /// knots on the edges u==0. & u==1. are obligatory
  ///
  /// The number of knots created and their values may change during initialisation:
  /// - Edge knots 0.f and 1.f will be added if they are not present.
  /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
  /// - Knots which are too close to each other will be merged
  /// - At least 5 knots and at least 4 axis bins will be created for consistency reason
  ///
  /// \param numberOfKnots     Number of knots in knots[] array
  /// \param knots             Array of knots.
  /// \param numberOfAxisBins Number of axis bins to map U coordinate to
  ///                          an appropriate [knot(i),knot(i+1)] interval.
  ///                          The knot positions have a "granularity" of 1./numberOfAxisBins
  ///

  //FlatObject::startConstruction();

  mNumberOfKnots = numberOfKnots < 5 ? 5 : numberOfKnots;
  //mNumberOfKnots = numberOfKnots > 255 ? 255 : numberOfKnots;

  //FlatObject::finishConstruction( mNumberOfKnots*sizeof(int) );

}

}// namespace
}// namespace

