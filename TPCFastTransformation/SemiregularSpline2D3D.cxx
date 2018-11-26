// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


/// \file  RegularSpline2D3D.cxx
/// \brief Implementation of RegularSpline2D3D class
///
/// modified by Felix Lapp 06.09.2018
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>


#include "SemiregularSpline2D3D.h"

namespace ali_tpc_common {
  namespace tpc_fast_transformation {

    SemiregularSpline2D3D::SemiregularSpline2D3D()
      :
      FlatObject(),
      mGridV(),
      mNumberOfRows(0),
      mNumberOfKnots(0),
      mDataIndexMapOffset(0)
    {
      /// Default constructor. Creates an empty uninitialised object
    }


    void SemiregularSpline2D3D::destroy()
    {
      mNumberOfRows = 0;
      mDataIndexMapOffset = 0;
      mNumberOfKnots = 0;
      FlatObject::destroy();
    }

    void SemiregularSpline2D3D::relocateBufferPointers( const char* oldBuffer, char *actualBuffer )
    {
      /// relocate pointers from old to new buffer location
 
      /*char *bufferV = FlatObject::relocatePointer( oldBuffer, actualBuffer, mGridV.getFlatBufferPtr() );
	mGridV.setActualBufferAddress( bufferV );*/

	/*for( int i=0; i<mNumberOfRows; i++) {
	char *bufferUi = FlatObject::relocatePointer(oldBuffer, actualBuffer, mSplineArray[i].getFlatBufferPtr() );
	mSplineArray[i].setActualBufferAddress( bufferUi );
	}*/

    }


    void SemiregularSpline2D3D::cloneFromObject( const SemiregularSpline2D3D &obj, char *newFlatBufferPtr )
    {
      /// See FlatObject for description      

      const char *oldFlatBufferPtr = obj.mFlatBufferPtr;
  
      FlatObject::cloneFromObject( obj, newFlatBufferPtr );
      mNumberOfRows = obj.mNumberOfRows;
      mNumberOfKnots = obj.mNumberOfKnots;
      mGridV = obj.mGridV;
      mDataIndexMapOffset = obj.mDataIndexMapOffset;
      /* char *bufferV = FlatObject::relocatePointer( oldFlatBufferPtr, mFlatBufferPtr, obj.mGridV.getFlatBufferPtr() );
	 mGridV.cloneFromObject( obj.mGridV, bufferV );*/

      /*for( int i=0; i<mNumberOfRows; i++) {
	char *bufferUi = FlatObject::relocatePointer( oldFlatBufferPtr, mFlatBufferPtr, obj.mSplineArray[i].getFlatBufferPtr() );
	mSplineArray[i].cloneFromObject( obj.mSplineArray[i], bufferUi );
	}*/
    }
 


    void SemiregularSpline2D3D::moveBufferTo( char *newFlatBufferPtr )
    {
      /// See FlatObject for description      
      const char *oldFlatBufferPtr = mFlatBufferPtr;
      FlatObject::moveBufferTo( newFlatBufferPtr );
      relocateBufferPointers( oldFlatBufferPtr, mFlatBufferPtr );
    }
  
  
    void SemiregularSpline2D3D::setActualBufferAddress( char* actualFlatBufferPtr )
    {
      /// See FlatObject for description      
      const char *oldFlatBufferPtr = mFlatBufferPtr;
      FlatObject::setActualBufferAddress( actualFlatBufferPtr );
      relocateBufferPointers( oldFlatBufferPtr, mFlatBufferPtr );
    }

  
    void SemiregularSpline2D3D::setFutureBufferAddress( char* futureFlatBufferPtr )
    {
      /// See FlatObject for description      
      /*const char* oldFlatBufferPtr = mFlatBufferPtr;
 
	char *bufferV = relocatePointer( oldFlatBufferPtr, futureFlatBufferPtr, mGridV.getFlatBufferPtr() );
	mGridV.setFutureBufferAddress( bufferV );*/

	/*for( int i=0; i<mNumberOfRows; i++ ) {
	char *bufferUi = relocatePointer( oldFlatBufferPtr, futureFlatBufferPtr, mSplineArray[i].getFlatBufferPtr() );
	mSplineArray[i].setFutureBufferAddress( bufferUi );
	}*/


      FlatObject::setFutureBufferAddress( futureFlatBufferPtr );
    }



    void SemiregularSpline2D3D::construct(const int numberOfRows, const int numbersOfKnots[] )
    {
      /// Constructor
      ///
      /// Number of knots created and their values may differ from the input values:
      /// - Edge knots 0.f and 1.f will be added if they are not present.
      /// - Knot values are rounded to closest axis bins: k*1./numberOfAxisBins.
      /// - Knots which are too close to each other will be merged
      /// - At least 5 knots and at least 4 axis bins will be created for consistency reason
      ///
      /// \param numberOfKnotsU     U axis: Number of knots in knots[] array 
      /// \param knotsU             U axis: Array of knots.
      /// \param numberOfAxisBinsU  U axis: Number of axis bins to map U coordinate to
      ///                           an appropriate [knot(i),knot(i+1)] interval.
      ///                           The knot positions have a "granularity" of 1./numberOfAxisBins
      ///
      /// \param numberOfKnotsV     V axis: Number of knots in knots[] array 
      /// \param knotsV             V axis: Array of knots.
      /// \param numberOfAxisBinsV  V axis: Number of axis bins to map U coordinate to
      ///                           an appropriate [knot(i),knot(i+1)] interval.
      ///                           The knot positions have a "granularity" of 1./numberOfAxisBins
      /// 

      ///  float f[nKnots] = { 3.5, 2.0, 1.4, 3.8, 2.3};
      ///  spline.correctEdges( f );
      ///  spline.getSpline( f, 0. ); // == 3.5
      ///  spline.getSpline( f, 0.1 ); // == some interpolated value
      ///  spline.getSpline( f, 0.25 ); // == 2.0 
      ///  spline.getSpline( f, 0.5  ); // == 1.4 
      ///  spline.getSpline( f, 1.  ); // == 2.3 

      /*if(numberOfRows < 5) {
	return;
	}*/
  


      FlatObject::startConstruction();

      //construct regular grid for v
      mGridV.construct( numberOfRows );

      // For each x element numbersOfKnots may be a single RegularSpline1D with x knots.
      // so first create the array
      RegularSpline1D splineArray[numberOfRows];

      // And construct them
      for( int i=0; i<numberOfRows; i++) {
	RegularSpline1D spline;                          // instantiate
	spline.construct(numbersOfKnots[i]);             // construct...
	splineArray[i] = spline;                         // save
      }

      // this is the space which is taken just by the RegularSpline1D's
      mDataIndexMapOffset = numberOfRows * sizeof(RegularSpline1D);

      //The buffer size is the size of the array
      FlatObject::finishConstruction(mDataIndexMapOffset + numberOfRows*sizeof(int));

      // Array for the 1D-Splines inside the buffer
      RegularSpline1D *bufferSplines = getSplineArrayNonConst();
      
      // paste local splineArray to the buffer
      for( int i=0; i<numberOfRows; i++) {
	bufferSplines[i] = splineArray[i];
      }

      // Just calculating the total number of knots in this 2D3D spline.
      int numberOfKnots = 0;
      for( int i=0; i < numberOfRows; i++) {
	int knotsU = getGridU(i).getNumberOfKnots();
	numberOfKnots += knotsU;
      }

      //save the numberOfRows and numberOfKnots
      mNumberOfRows = numberOfRows;
      mNumberOfKnots = numberOfKnots;

    // map to save the starting data index for each v-coordinate
    int* dataIndexMap = getDataIndexMapNonConst();

    // this will count the amount of u-knots "under" a v-coordinate
    int uSum = 0;

      //count the amount of knots which are in gridU's lower than i
      for( int dv=0; dv<mNumberOfRows; dv++ ) {
	dataIndexMap[dv] = uSum;
	uSum += numbersOfKnots[dv];
      }
    }




  }// namespace
}// namespace
