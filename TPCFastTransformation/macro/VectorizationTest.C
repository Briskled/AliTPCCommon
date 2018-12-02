/*

// works only with ROOT >= 6

alienv load ROOT/latest-root6
alienv load Vc/latest

root -l  loadlibs.C
.x RegularSpline2D3DTest.C++
*/

#include "TFile.h"
#include "TRandom.h"
#include "TNtuple.h"
#include "Riostream.h"
#include "TSystem.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "SemiregularSpline2D3D.h"
#include "IrregularSpline2D3D.h"
#include "RegularSpline1D.h"
#include "IrregularSpline1D.h"

#if !defined(__CINT__) && !defined(__ROOTCINT__) && !defined(HLTCA_GPUCODE)
#include <Vc/Vc>
#endif

const int PolynomDegree = 7;

double cu[PolynomDegree+1], cv[PolynomDegree+1];

float Fx( float u, float v )
{
  u-=0.5;
  v-=0.5;
  double uu = 1.;  
  double vv = 1.;  
  double f = 0;
  for( int i=0; i<=PolynomDegree; i++ ){
    f+=cu[i]*uu;
    f+=cv[i]*vv;
    uu*=u;
    vv*=v;

  }
  return f;
}

float Fy( float u, float v ){  return v; }
float Fz( float u, float v ){  return (u-.5)*(u-.5); }


int VectorizationTest()
{
  using namespace ali_tpc_common::tpc_fast_transformation ;

  const int numberOfRows = 6;
  const int numbersOfKnots[numberOfRows] = {6, 6, 6, 6, 6, 6};

  const int numberOfKnotsU = 6;
  const int numberOfKnotsV = 6;
  const float knotsUIRS[numberOfKnotsU] = {0., 0.2, 0.4, 0.6, 0.8, 1.};
  const float knotsVIRS[numberOfKnotsV] = {0., 0.2, 0.4, 0.6, 0.8, 1.};

  const int numberOfKnots1D = 6;
  const float knots1D[numberOfKnots1D] = {0., 0.2, 0.4, 0.6, 0.8, 1.};


  //======================================================
  cout << "Constructing SRS..." << endl;
  SemiregularSpline2D3D splineSRS;
  splineSRS.construct( numberOfRows, numbersOfKnots );

  int nKnotsTotSRS = splineSRS.getNumberOfKnots();
  const RegularSpline1D &gridV = splineSRS.getGridV();

  float *dataSRS = new float[ 3*nKnotsTotSRS ];

  for( int i=0; i<gridV.getNumberOfKnots(); i++ ){
    double v = gridV.indexToU(i);
    const RegularSpline1D &gridU = splineSRS.getGridU(i);
    for( int j=0; j<gridU.getNumberOfKnots(); j++ ){  
      double u =  gridU.indexToU(j);
      int ind = splineSRS.getDataIndex(j, i);
      dataSRS[ind+0] = Fx(u,v);
      dataSRS[ind+1] = Fy(u,v);
      dataSRS[ind+2] = Fz(u,v);
      dataSRS[ind+0] = gRandom->Uniform(-1,1);//Gaus();
    }
  }

  splineSRS.correctEdges(dataSRS);



  //=======================================================
  cout << "Constructing IRS..." << endl;
  IrregularSpline2D3D splineIRS;
  splineIRS.construct( numberOfKnotsU, knotsUIRS, numberOfKnotsU-1,
		       numberOfKnotsV, knotsVIRS, numberOfKnotsV-1 );


  const IrregularSpline1D &gridU = splineIRS.getGridU();
  const IrregularSpline1D &gridVIRS = splineIRS.getGridV();      
  const int nKnotsTotIRS = splineIRS.getNumberOfKnots();

  float *dataIRS = dataSRS;
  splineIRS.correctEdges( dataIRS );


  //======================================================
  cout << "Constructing RS1D..." << endl;
  RegularSpline1D splineRS1D;
  splineRS1D.construct(numberOfKnots1D);

  float *dataRS1D = new float[ splineRS1D.getNumberOfKnots() ]; // corrected data
 
  for( int i=0; i<splineRS1D.getNumberOfKnots(); i++ ){
    dataRS1D[i] = gRandom->Uniform(-1.,1.);
  }
  
  splineRS1D.correctEdges(dataRS1D);


  //=====================================================
  cout << "Constructing IRS1D.." << endl;
  IrregularSpline1D splineIRS1D;
  splineIRS1D.construct(numberOfKnots1D, knots1D, numberOfKnots1D-1);
  
  splineIRS1D.correctEdges(dataRS1D);

  cout << "Starting to measure time..." << endl;



  float stepu = 1.e-4; // 1.e-8;
  long borderU = 1./stepu;
  float stepv = 1.e-4;
  long borderV = 1./stepu;
  float step1D = 1.e-8;
  long border1D = 1./step1D;

  cout << "Current step 2D: " << stepu << endl;
  cout << "Current step 1D: " << step1D << endl;



  cout << "Measuring time for SRS getSpline Normal..." << endl;
  float sum1=0;
  TStopwatch timerSRSnormal;
  {
    float u=0, v=0;
    for( long i=0; i<borderU; i++, u+=stepu ){
      for( long j=0; j<borderV; j++, v+=stepv ){
	float x, y, z;
	splineSRS.getSpline( dataSRS, u,v,x,y,z);
	sum1 += x+y+z;
      }
    }
    timerSRSnormal.Stop();
    cout << "Sum1: " << sum1 << endl;
  }

  cout << "Measuring time for SRS getSpline Normal unoptimized..." << endl;
  float sum2=0;
  TStopwatch timerSRSunopt;
  {
    float u=0, v=0;
    for( long i=0; i<borderU; i++, u+=stepu ){
      for( long j=0; j<borderV; j++, v+=stepv ){
	float x, y, z;
	splineSRS.getSplineUnoptimized( dataSRS, u,v,x,y,z);
	sum2 += x+y+z;
      }
    }
    timerSRSunopt.Stop();
    cout << "Sum2: " << sum2 << endl;
  }

  cout << "Measuring time for SRS getSpline Vectorized..." << endl;
  float sum3=0;
  TStopwatch timerSRSvec;
  {
    float u=0, v=0;
    for( long i=0; i<borderU; i++, u+=stepu ){
      for( long j=0; j<borderV; j++, v+=stepv ){
	float x, y, z;
	splineSRS.getSplineVec( dataSRS, u,v,x,y,z);
	sum3 += x+y+z;
      }
    }
    timerSRSvec.Stop();
    cout << "Sum3: " << sum3 << endl;
  }


  cout << "Measuring time for IRS getSpline Normal..." << endl;
  float sum4=0;
  TStopwatch timerIRSnormal;
  {
    float u=0, v=0;
    for( long i=0; i<borderU; i++, u+=stepu ){
      for( long j=0; j<borderV; j++, v+=stepv ){
	float x, y, z;
	splineIRS.getSpline( dataIRS, u,v,x,y,z);
	sum4 += x+y+z;
      }
    }
    timerIRSnormal.Stop();
    cout << "Sum4: " << sum4 << endl;
  }


  cout << "Measuring time for IRS getSpline Vectorized..." << endl;
  float sum5=0;
  TStopwatch timerIRSvec;
  {
    float u=0, v=0;
    for( long i=0; i<borderU; i++, u+=stepu ){
      for( long j=0; j<borderV; j++, v+=stepv ){
	float x, y, z;
	splineIRS.getSplineVec( dataIRS, u,v,x,y,z);
	sum5 += x+y+z;
      }
    }
  }
  timerIRSvec.Stop();
  cout << "Sum5: " << sum5 << endl;
  


  cout << "Measuring time for RS1D getSpline..." << endl;


  /*std::vector<float> tms;
    for(int i=0; i<100; i++) {
    float sum6=0;
    TStopwatch timerRS1D;
    {
    //for( float u=-0.01; u<=1.01; u+= step1D) {
    float u=0;
    for( long i=0; i<border1D; i++, u+= step1D) {
    sum6 += splineRS1D.getSpline(dataRS1D, u);
    }
    }
    timerRS1D.Stop();
    tms.push_back(timerRS1D.RealTime());
    cout << "Sum6: " << sum6 << endl;
    cout << "Time taken: " << timerRS1D.RealTime()*1.e3 << " ms" << endl;
    }

    float avr = 0;
    for(unsigned int i=0; i<tms.size(); i++) {
    avr += tms[i];
    }
    avr /= tms.size();
    cout << "AVR: " << avr << endl;*/

  float sum6=0;
  TStopwatch timerRS1D;
  {
    //for( float u=-0.01; u<=1.01; u+= step1D) {
    float u=0;
    for( long i=0; i<border1D; i++, u+= step1D) {
      sum6 += splineRS1D.getSpline(dataRS1D, u);
    }
  }
  timerRS1D.Stop();
  cout << "Sum6: " << sum6 << endl;
  


  cout << "Measuring time for IRS1D getSpline..." << endl;
  float sum7=0;
  TStopwatch timerIRS1D;
  {
    //for( float u=-0.01; u<=1.01; u+= step1D) {
    float u=0;
    for( long i=0; i<border1D; i++, u+= step1D) {
      sum7 += splineIRS1D.getSpline(dataRS1D, u);
    }
  }
  timerIRS1D.Stop();
  cout << "Sum7: " << sum7 << endl;


  
  //============================================
  

  cout << "Constructing 1D-Splines for 2D-Simulation" << endl;
  int gridUSize = 12;
  float dataSim[gridUSize*4*4];
  for( int i=0; i<gridUSize*4*4; i++) dataSim[i] = 1;

  /*  for( int i=0; i<4; i++ ){
    double v = i/3.;
    for( int j=0; j<gridUSize*4; j++ ){  
      double u = j/((double)gridUSize-1);
      int ind = i*gridUSize+j;
      dataSim[ind+0] = Fx(u,v);
      dataSim[ind+1] = Fy(u,v);
      dataSim[ind+2] = Fz(u,v);
      dataSim[ind+3] = 0;
      dataSim[ind+0] = gRandom->Uniform(-1,1);//Gaus();
    }
    }*/

  cout << "Simulating 2D SRS by using 4 RegularSplines with " << gridUSize << " knots." << endl;
  cout << "Simulating rows first..." << endl;
  float sum8=0;
  long iterRow=0;
  float stepUU = stepu;
  TStopwatch timerSimRow;
  {      
    float *dt0 = dataSim;
    float *dt1 = dt0 + gridUSize*3;
    float *dt2 = dt1 + gridUSize*3;
    float *dt3 = dt2 + gridUSize*3;
    for( int iter=0; iter<10000; iter++)
    for( float u=-0.01; u<=1.01; u+=stepUU ){
      for( int k=0; k<3; k++) {
	float x = RegularSpline1D::getSpline(1, gridUSize, dt0[0+k], dt0[3+k], dt0[6+k], dt0[9+k], u);
	float y = RegularSpline1D::getSpline(1, gridUSize, dt1[0+k], dt1[3+k], dt1[6+k], dt1[9+k], u);
	float z = RegularSpline1D::getSpline(1, gridUSize, dt2[0+k], dt2[3+k], dt2[6+k], dt2[9+k], u);
	float w = RegularSpline1D::getSpline(1, gridUSize, dt3[0+k], dt3[3+k], dt3[6+k], dt3[9+k], u);
	sum8 += x+y+z+w;
      }
    }
  }
  timerSimRow.Stop();
  cout << "Sum8: " << sum8 << endl;
  cout << "Iterations: " << iterRow << endl;






  cout << "Simulating columns first..." << endl;
  float sum9=0;
  long iterCol=0;
  TStopwatch timerSimCol;
  {
    float *dt0 = dataSim;
    float *dt1 = dt0 + gridUSize*3;
    float *dt2 = dt1 + gridUSize*3;
    float *dt3 = dt2 + gridUSize*3;
     
    for( int iter=0; iter<10000; iter++)
    for( float u=-0.01; u<=1.01; u+=stepUU ){
      for( int k=0; k<12; k++) {
	float x = RegularSpline1D::getSpline(1, 4, dt0[k], dt1[k], dt2[k], dt3[k], u);
	sum9 += x;
      }
    }
  }
  timerSimCol.Stop();
  cout << "Sum9: " << sum9 << endl;
  cout << "Iterations: " << iterCol << endl;
  

  cout << endl;
  cout << "SRS Normal:\t\t" << timerSRSnormal.RealTime()*1.e3 << " ms" << endl;
  cout << "SRS Normal Unoptimized:\t" << timerSRSunopt.RealTime()*1.e3 << " ms" << endl;
  cout << "SRS Vectorized:\t\t" << timerSRSvec.RealTime()*1.e3 << " ms" << endl;
  cout << endl;

  cout << "IRS Normal:\t\t" << timerIRSnormal.RealTime()*1.e3 << " ms" << endl;
  cout << "IRS Vectorized:\t\t" << timerIRSvec.RealTime()*1.e3 << " ms" << endl;
  cout << endl;

  cout << "RS1D:\t\t\t" << timerRS1D.RealTime()*1.e3 << " ms" << endl;
  cout << "IRS1D:\t\t\t" << timerIRS1D.RealTime()*1.e3 << " ms" << endl;
  cout << endl;

  cout << "Simulation rows first:\t" << timerSimRow.RealTime()*1.e3 << " ms" << endl;
  cout << "Simulation cols first:\t" << timerSimCol.RealTime()*1.e3 << " ms" << endl;
  cout << endl;

  cout << "SRS Optimization speedup:\t" << timerSRSunopt.RealTime()/timerSRSnormal.RealTime() << endl;
  cout << "SRS Vectorization speedup:\t" << timerSRSnormal.RealTime()/timerSRSvec.RealTime() << endl;
  cout << "IRS Vectorization speedup:\t" << timerIRSnormal.RealTime()/timerIRSvec.RealTime() << endl;
  cout << "IRS to SRS Normal speedup:\t" << timerIRSnormal.RealTime()/timerSRSnormal.RealTime() << endl;
  cout << "IRS to SRS Vectorized speedup:\t" << timerIRSvec.RealTime()/timerSRSvec.RealTime() << endl;
  cout << "IRS1D to RS1D speedup:\t\t" << timerIRS1D.RealTime()/timerRS1D.RealTime() << endl;
  cout << "1D-2D-ratio for SRS:\t\t" << timerSRSnormal.RealTime()/timerRS1D.RealTime() << endl;
  cout << "1D-2D-ratio for IRS:\t\t" << timerIRSnormal.RealTime()/timerIRS1D.RealTime() << endl;

  return 0;
}
