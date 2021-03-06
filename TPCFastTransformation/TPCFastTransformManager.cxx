// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


/// \file  TPCFastTransformManager.cxx
/// \brief Implementation of TPCFastTransformManager class
///
/// \author  Sergey Gorbunov <sergey.gorbunov@cern.ch>


#include "TPCFastTransformManager.h"
#include "AliTPCTransform.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include "AliTPCcalibDB.h"
#include "AliHLTTPCGeometry.h"
#include "TPCFastTransform.h"
#include "TPCFastTransformSRS.h"

namespace ali_tpc_common {
namespace tpc_fast_transformation {


TPCFastTransformManager::TPCFastTransformManager()
  :
  mError(),
  mOrigTransform(nullptr)
{
}


int  TPCFastTransformManager::create( TPCFastTransform &fastTransform, AliTPCTransform *transform, Long_t TimeStamp )
{
  /// Initializes TPCFastTransform object
 
  mOrigTransform = transform;

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();
  if(!pCalib ) return storeError( -1, "TPCFastTransformManager::Init: No TPC calibration instance found");

  AliTPCParam *tpcParam = pCalib->GetParameters();
  if( !tpcParam ) return storeError( -2, "TPCFastTransformManager::Init: No TPCParam object found");

  // if there is no fastTransform: use default one given by pCalib
  if( !transform ) transform = pCalib->GetTransform();
  if( !transform ) return storeError( -3, "TPCFastTransformManager::Init: No TPC transformation found");
 
  tpcParam->Update();
  tpcParam->ReadGeoMatrices();

  const AliTPCRecoParam *rec = transform->GetCurrentRecoParam();
  if( !rec ) return storeError( -5, "TPCFastTransformManager::Init: No TPC Reco Param set in transformation");

  bool useCorrectionMap = rec->GetUseCorrectionMap();

  if ( useCorrectionMap) transform->SetCorrectionMapMode(kTRUE); //If the simulation set this to false to simulate distortions, we need to reverse it for the transformation


  // find last calibrated time bin
  
  fLastTimeBin = rec->GetLastBin();  
  
  fastTransform.startConstruction( tpcParam->GetNRowLow()+ tpcParam->GetNRowUp() );

  TPCDistortionIRS& distortion = fastTransform.getDistortionNonConst();
  
  distortion.startConstruction( tpcParam->GetNRowLow()+ tpcParam->GetNRowUp(), 1 );
  
  float tpcZlengthSideA = tpcParam->GetZLength(0);
  float tpcZlengthSideC = tpcParam->GetZLength(TPCFastTransform::getNumberOfSlices()/2);

  fastTransform.setTPCgeometry(  tpcZlengthSideA, tpcZlengthSideC );
  distortion.setTPCgeometry(  tpcZlengthSideA, tpcZlengthSideC );
 
  for( int iRow=0; iRow<fastTransform.getNumberOfRows(); iRow++){
    int sector=0, secrow=0;
    AliHLTTPCGeometry::Slice2Sector(0, iRow, sector, secrow);
    Int_t nPads = tpcParam->GetNPads( sector, secrow );
    float xRow = tpcParam->GetPadRowRadii( sector,  secrow);
    float padWidth = tpcParam->GetInnerPadPitchWidth();
    if( iRow >= tpcParam->GetNRowLow() ) padWidth = tpcParam->GetOuterPadPitchWidth();
    fastTransform.setTPCrow( iRow, xRow, nPads, padWidth );
    distortion.setTPCrow( iRow, xRow, nPads, padWidth, 0 );
  }

  fastTransform.setCalibration( -1, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f );

  IrregularSpline2D3D spline;
  {
    int nKnotsU = 15;
    int nAxisTicksU = tpcParam->GetNPads(0, 10);    
    int nKnotsV = 20;
    int nAxisTicksV = fLastTimeBin+1;
    float knotsU[nKnotsU];
    float knotsV[nKnotsV];
    for( int i=0; i<nKnotsU; i++ ) knotsU[i] = 1./(nKnotsU-1)*i;
    for( int i=0; i<nKnotsV; i++ ) knotsV[i] = 1./(nKnotsV-1)*i;
    
    // make bining similar to old HLT transformation 
    // TODO: adjust to binning in the calibration   

    double d1 = 0.6;
    double d2 = 0.9 - d1;
    double d3 = 1.-d2 - d1;
    
    /*for( int i=0; i<5; i++ ){ // 5 bins in first 6% of drift
      knotsV[i] = i / 4. * d1;
    }
    for( int i=0; i<10; i++ ){ // 10 bins for 6% <-> 90%
      knotsV[4+i] = d1 + i/9. * d2;
    }
    for( int i=0; i<5; i++ ){ // 5 bins for last 90% <-> 100%
      knotsV[13+i] = d1 + d2 + i/4. * d3;
      }*/

    spline.construct( nKnotsU, knotsU, nAxisTicksU,
		      nKnotsV, knotsV, nAxisTicksV );
  }
  distortion.setApproximationScenario( 0, spline );  
  distortion.finishConstruction();
  fastTransform.finishConstruction();

  return updateCalibration( fastTransform, TimeStamp );
}






int  TPCFastTransformManager::createSRS( TPCFastTransformSRS &fastTransform, AliTPCTransform *transform, Long_t TimeStamp )
{
  /// Initializes TPCFastTransform object
 
  mOrigTransform = transform;

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();
  if(!pCalib ) return storeError( -1, "TPCFastTransformManager::Init: No TPC calibration instance found");

  AliTPCParam *tpcParam = pCalib->GetParameters();
  if( !tpcParam ) return storeError( -2, "TPCFastTransformManager::Init: No TPCParam object found");

  // if there is no fastTransform: use default one given by pCalib
  if( !transform ) transform = pCalib->GetTransform();
  if( !transform ) return storeError( -3, "TPCFastTransformManager::Init: No TPC transformation found");
 
  tpcParam->Update();
  tpcParam->ReadGeoMatrices();

  const AliTPCRecoParam *rec = transform->GetCurrentRecoParam();
  if( !rec ) return storeError( -5, "TPCFastTransformManager::Init: No TPC Reco Param set in transformation");

  bool useCorrectionMap = rec->GetUseCorrectionMap();

  if ( useCorrectionMap) transform->SetCorrectionMapMode(kTRUE); //If the simulation set this to false to simulate distortions, we need to reverse it for the transformation


  // find last calibrated time bin
  
  fLastTimeBin = rec->GetLastBin();  
  

  fastTransform.startConstruction( tpcParam->GetNRowLow()+ tpcParam->GetNRowUp() );

  TPCDistortionSRS& distortion = fastTransform.getDistortionNonConst();

  int numberOfSlices = distortion.getNumberOfSlices();
  distortion.startConstruction( tpcParam->GetNRowLow()+ tpcParam->GetNRowUp(), fastTransform.getNumberOfRows()*numberOfSlices );
  
  float tpcZlengthSideA = tpcParam->GetZLength(0);
  float tpcZlengthSideC = tpcParam->GetZLength(TPCFastTransform::getNumberOfSlices()/2);

  AliTPCChebCorr* chebCorr = transform->GetCorrMapCacheRef();

  fastTransform.setTPCgeometry(  tpcZlengthSideA, tpcZlengthSideC );
  distortion.setTPCgeometry(  tpcZlengthSideA, tpcZlengthSideC );

  for( int iRow=0; iRow<fastTransform.getNumberOfRows(); iRow++){
    int sector=0, secrow=0;
    AliHLTTPCGeometry::Slice2Sector(0, iRow, sector, secrow);
    Int_t nPads = tpcParam->GetNPads( sector, secrow );
    float xRow = tpcParam->GetPadRowRadii( sector,  secrow);
    float padWidth = tpcParam->GetInnerPadPitchWidth();
    if( iRow >= tpcParam->GetNRowLow() ) padWidth = tpcParam->GetOuterPadPitchWidth();
    fastTransform.setTPCrow( iRow, xRow, nPads, padWidth );
    distortion.setTPCrow( iRow, xRow, nPads, padWidth );
  }   



  for (int iSlice=0; iSlice<numberOfSlices; iSlice++) {
    for( int iRow=0; iRow<fastTransform.getNumberOfRows(); iRow++){
      int sector=0, secrow=0;
      AliHLTTPCGeometry::Slice2Sector(iSlice, iRow, sector, secrow);
      Int_t nPads = tpcParam->GetNPads( sector, secrow );
      float xRow = tpcParam->GetPadRowRadii( sector,  secrow);
      float padWidth = tpcParam->GetInnerPadPitchWidth();
      if( iRow >= tpcParam->GetNRowLow() ) padWidth = tpcParam->GetOuterPadPitchWidth();

      SemiregularSpline2D3D spline;
      {      
	int numberOfSplines = 20;
	int numbersOfKnots[numberOfSplines];
	for( int i=0; i<numberOfSplines; i++) {
	  numbersOfKnots[i] = 15;
	}
	spline.construct(numberOfSplines, numbersOfKnots);
      }
      distortion.setApproximationScenario(iSlice*fastTransform.getNumberOfRows() + iRow, spline );      
    }
  }
  
  fastTransform.setCalibration( -1, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f );

  
  distortion.finishConstruction();
  fastTransform.finishConstruction();

  return updateCalibrationSRS( fastTransform, TimeStamp );
}




  
int TPCFastTransformManager::updateCalibration( TPCFastTransform &fastTransform, Long_t TimeStamp )
{
  // Update the calibration with the new time stamp

  Long_t lastTS = fastTransform.getTimeStamp();
  
  // deinitialize

  fastTransform.setTimeStamp( -1 ); 

  if( TimeStamp < 0  ) return 0; 

  // search for the calibration database

  if( !mOrigTransform ) return storeError( -1, "TPCFastTransformManager::SetCurrentTimeStamp: TPC transformation has not been set properly"); 

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();  
  if(!pCalib ) return storeError( -2, "TPCFastTransformManager::SetCurrentTimeStamp: No TPC calibration found"); 
  
  AliTPCParam *tpcParam = pCalib->GetParameters(); 
  if( !tpcParam ) return  storeError( -3, "TPCFastTransformManager::SetCurrentTimeStamp: No TPCParam object found");   

  AliTPCRecoParam *recoParam = mOrigTransform->GetCurrentRecoParamNonConst();
  if( !recoParam ) return storeError( -5, "TPCFastTransformManager::Init: No TPC Reco Param set in transformation");
  
  // calibration found, set the initialized status back

  fastTransform.setTimeStamp(lastTS);

  // less than 60 seconds from the previois time stamp, don't do anything

  if( lastTS>=0 && TMath::Abs(lastTS - TimeStamp ) <60 ) return 0; 
  
  // start the initialization
  
  bool useCorrectionMap = recoParam->GetUseCorrectionMap();
  
  if( useCorrectionMap ){
    // If the simulation set this to false to simulate distortions, we need to reverse it for the transformation
    // This is a design feature. Historically HLT code runs as a part of simulation, not reconstruction.
    mOrigTransform->SetCorrectionMapMode(kTRUE);
  }

  // set the current time stamp
  
  mOrigTransform->SetCurrentTimeStamp( static_cast<UInt_t>(TimeStamp) );
  fastTransform.setTimeStamp(TimeStamp);
  
  // find last calibrated time bin
  
  fLastTimeBin = recoParam->GetLastBin();
  
  double t0 = mOrigTransform->GetTBinOffset();
  double driftCorrPT = mOrigTransform->GetDriftCorrPT();
  double vdCorrectionTime = mOrigTransform->GetVDCorrectionTime();
  double vdCorrectionTimeGY = mOrigTransform->GetVDCorrectionTimeGY();
  double time0CorrTime = mOrigTransform->GetTime0CorrTime();

  // original formula:
  // L = (t-t0)*ZWidth*driftCorrPT*vdCorrectionTime*( 1 + yLab*vdCorrectionTimeGY )  -  time0CorrTime + 3.*tpcParam->GetZSigma();
  // Z = Z(L) - fDeltaZCorrTime
  // chebyshev distortions for xyz
  // Time-of-flight correction: ldrift += dist-to-vtx*tofCorr

  // fast transform formula:
  // L = (t-t0)*(mVdrift + mVdriftCorrY*yLab ) + mLdriftCorr 
  // Z = Z(L) +  tpcAlignmentZ
  // spline distortions for xyz
  // Time-of-flight correction: ldrift += dist-to-vtx*tofCorr

  double vDrift = tpcParam->GetZWidth()*driftCorrPT*vdCorrectionTime;
  double vdCorrY = vDrift*vdCorrectionTimeGY;
  double ldCorr = - time0CorrTime + 3*tpcParam->GetZSigma();

  double tpcAlignmentZ = - mOrigTransform->GetDeltaZCorrTime();

  double tofCorr = ( 0.01*tpcParam->GetDriftV() )/TMath::C();
  double primVtxZ = mOrigTransform->GetPrimVertex()[2];

  bool useTOFcorrection = recoParam->GetUseTOFCorrection();

  if( !useTOFcorrection ) tofCorr = 0;

  fastTransform.setCalibration( TimeStamp, t0, vDrift, vdCorrY, ldCorr, tofCorr, primVtxZ, tpcAlignmentZ);

  // now calculate distortion map: dx,du,dv = ( origTransform() -> x,u,v) - fastTransformNominal:x,u,v

  TPCDistortionIRS& distortion = fastTransform.getDistortionNonConst();
  
  // switch TOF correction off for a while

  recoParam->SetUseTOFCorrection( kFALSE );
 
  for( int slice=0; slice<distortion.getNumberOfSlices(); slice++){      
  
    for( int row=0; row<distortion.getNumberOfRows(); row++ ){

      const TPCFastTransform::RowInfo &rowInfo = fastTransform.getRowInfo( row );
    
      const IrregularSpline2D3D& spline = distortion.getSpline( slice, row );

      float *data = distortion.getSplineDataNonConst(slice,row);

      for( int knot=0; knot<spline.getNumberOfKnots(); knot++ ){

	data[3*knot+0] = 0.f;
	data[3*knot+1] = 0.f;
	data[3*knot+2] = 0.f;
	
	// x cordinate of the knot
	float x = rowInfo.x;

	// spline (su,sv) cordinates of the knot (su,sv) in (0,1)x(0,1)
	float su=0, sv=0;
	spline.getKnotUV( knot, su, sv );
	
	// x, u, v cordinates of the knot (local cartesian coord. of slice towards central electrode )
	float u=0, v=0;
	distortion.convSUVtoUV( slice, row, su, sv, u, v );

	// row, pad, time coordinates of the knot 
	float pad=0, time=0;
	fastTransform.convUVtoPadTime( slice, row, u, v, pad, time );
	
	// nominal x,y,z coordinates of the knot (without distortions and time-of-flight correction)
	float y=0, z=0;
	fastTransform.convUVtoYZ( slice, row, x, u, v, y, z );
	
	// original TPC transformation (row,pad,time) -> (x,y,z) without time-of-flight correction
	float ox=0, oy=0, oz=0;	
	{
	  int sector=0, secrow=0;
	  AliHLTTPCGeometry::Slice2Sector( slice, row, sector, secrow );
	  int is[]={sector};
	  double xx[]={ static_cast<double>(secrow), pad, time };
	  mOrigTransform->Transform(xx, is, 0, 1);
	  ox = xx[0];
	  oy = xx[1];
	  oz = xx[2];
	}
	// convert to u,v
	float ou=0, ov=0;
	fastTransform.convYZtoUV( slice, row, ox, oy, oz, ou, ov );

	// distortions in x,u,v:
	float dx = ox-x;
	float du = ou-u;
	float dv = ov-v;
	data[3*knot+0] = dx;
	data[3*knot+1] = du;
	data[3*knot+2] = dv;	
      } // knots
      
      spline.correctEdges(data);
    } // row
  } // slice
  
  
  // set back the time-of-flight correction;
  
  recoParam->SetUseTOFCorrection( useTOFcorrection );
  
  return 0;
}



int TPCFastTransformManager::updateCalibrationSRS( TPCFastTransformSRS &fastTransform, Long_t TimeStamp )
{
  // Update the calibration with the new time stamp

  Long_t lastTS = fastTransform.getTimeStamp();
  
  // deinitialize

  fastTransform.setTimeStamp( -1 ); 

  if( TimeStamp < 0  ) return 0; 

  // search for the calibration database

  if( !mOrigTransform ) return storeError( -1, "TPCFastTransformManager::SetCurrentTimeStamp: TPC transformation has not been set properly"); 

  AliTPCcalibDB* pCalib=AliTPCcalibDB::Instance();  
  if(!pCalib ) return storeError( -2, "TPCFastTransformManager::SetCurrentTimeStamp: No TPC calibration found"); 
  
  AliTPCParam *tpcParam = pCalib->GetParameters(); 
  if( !tpcParam ) return  storeError( -3, "TPCFastTransformManager::SetCurrentTimeStamp: No TPCParam object found");   

  AliTPCRecoParam *recoParam = mOrigTransform->GetCurrentRecoParamNonConst();
  if( !recoParam ) return storeError( -5, "TPCFastTransformManager::Init: No TPC Reco Param set in transformation");
  
  // calibration found, set the initialized status back

  fastTransform.setTimeStamp(lastTS);

  // less than 60 seconds from the previois time stamp, don't do anything

  if( lastTS>=0 && TMath::Abs(lastTS - TimeStamp ) <60 ) return 0; 
  
  // start the initialization
  
  bool useCorrectionMap = recoParam->GetUseCorrectionMap();
  
  if( useCorrectionMap ){
    // If the simulation set this to false to simulate distortions, we need to reverse it for the transformation
    // This is a design feature. Historically HLT code runs as a part of simulation, not reconstruction.
    mOrigTransform->SetCorrectionMapMode(kTRUE);
  }

  // set the current time stamp
  
  mOrigTransform->SetCurrentTimeStamp( static_cast<UInt_t>(TimeStamp) );
  fastTransform.setTimeStamp(TimeStamp);
  
  // find last calibrated time bin
  
  fLastTimeBin = recoParam->GetLastBin();
  
  double t0 = mOrigTransform->GetTBinOffset();
  double driftCorrPT = mOrigTransform->GetDriftCorrPT();
  double vdCorrectionTime = mOrigTransform->GetVDCorrectionTime();
  double vdCorrectionTimeGY = mOrigTransform->GetVDCorrectionTimeGY();
  double time0CorrTime = mOrigTransform->GetTime0CorrTime();

  // original formula:
  // L = (t-t0)*ZWidth*driftCorrPT*vdCorrectionTime*( 1 + yLab*vdCorrectionTimeGY )  -  time0CorrTime + 3.*tpcParam->GetZSigma();
  // Z = Z(L) - fDeltaZCorrTime
  // chebyshev distortions for xyz
  // Time-of-flight correction: ldrift += dist-to-vtx*tofCorr

  // fast transform formula:
  // L = (t-t0)*(mVdrift + mVdriftCorrY*yLab ) + mLdriftCorr 
  // Z = Z(L) +  tpcAlignmentZ
  // spline distortions for xyz
  // Time-of-flight correction: ldrift += dist-to-vtx*tofCorr

  double vDrift = tpcParam->GetZWidth()*driftCorrPT*vdCorrectionTime;
  double vdCorrY = vDrift*vdCorrectionTimeGY;

  //TODO Change here back to 3* if the semiregular splines are using 3 floats for each point again
  double ldCorr = - time0CorrTime + 3*tpcParam->GetZSigma();

  double tpcAlignmentZ = - mOrigTransform->GetDeltaZCorrTime();

  double tofCorr = ( 0.01*tpcParam->GetDriftV() )/TMath::C();
  double primVtxZ = mOrigTransform->GetPrimVertex()[2];

  bool useTOFcorrection = recoParam->GetUseTOFCorrection();

  if( !useTOFcorrection ) tofCorr = 0;

  fastTransform.setCalibration( TimeStamp, t0, vDrift, vdCorrY, ldCorr, tofCorr, primVtxZ, tpcAlignmentZ);

  // now calculate distortion map: dx,du,dv = ( origTransform() -> x,u,v) - fastTransformNominal:x,u,v

  TPCDistortionSRS& distortion = fastTransform.getDistortionNonConst();
  
  // switch TOF correction off for a while

  recoParam->SetUseTOFCorrection( kFALSE );
 
  for( int slice=0; slice<distortion.getNumberOfSlices(); slice++){      
  
    for( int row=0; row<distortion.getNumberOfRows(); row++ ){

      const TPCFastTransformSRS::RowInfo &rowInfo = fastTransform.getRowInfo( row );
    
      const SemiregularSpline2D3D& spline = distortion.getSpline( slice, row );

      float *data = distortion.getSplineDataNonConst(slice,row);

      for( int knot=0; knot<spline.getNumberOfKnots(); knot++ ){

	data[3*knot+0] = 0.f;
	data[3*knot+1] = 0.f;
	data[3*knot+2] = 0.f;
	
	// x cordinate of the knot
	float x = rowInfo.x;

	// spline (su,sv) cordinates of the knot (su,sv) in (0,1)x(0,1)
	float su=0, sv=0;
	spline.getKnotUV( knot, su, sv );
	
	// x, u, v cordinates of the knot (local cartesian coord. of slice towards central electrode )
	float u=0, v=0;
	distortion.convSUVtoUV( slice, row, su, sv, u, v );

	// row, pad, time coordinates of the knot 
	float pad=0, time=0;
	fastTransform.convUVtoPadTime( slice, row, u, v, pad, time );
	
	// nominal x,y,z coordinates of the knot (without distortions and time-of-flight correction)
	float y=0, z=0;
	fastTransform.convUVtoYZ( slice, row, x, u, v, y, z );
	
	// original TPC transformation (row,pad,time) -> (x,y,z) without time-of-flight correction
	float ox=0, oy=0, oz=0;	
	{
	  int sector=0, secrow=0;
	  AliHLTTPCGeometry::Slice2Sector( slice, row, sector, secrow );
	  int is[]={sector};
	  double xx[]={ static_cast<double>(secrow), pad, time };
	  mOrigTransform->Transform(xx, is, 0, 1);
	  ox = xx[0];
	  oy = xx[1];
	  oz = xx[2];
	}
	// convert to u,v
	float ou=0, ov=0;
	fastTransform.convYZtoUV( slice, row, ox, oy, oz, ou, ov );

	// distortions in x,u,v:
	float dx = ox-x;
	float du = ou-u;
	float dv = ov-v;
	data[3*knot+0] = dx;
	data[3*knot+1] = du;
	data[3*knot+2] = dv;	
      } // knots
      
      spline.correctEdges(data);
    } // row
  } // slice
  
  
  // set back the time-of-flight correction;
  
  recoParam->SetUseTOFCorrection( useTOFcorrection );
  
  return 0;
}



 

}} // namespaces
