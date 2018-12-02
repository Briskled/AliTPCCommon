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

#include "SemiregularSpline2D3D.h"

#ifndef Spline2D3DTest_H
#define Spline2D3DTest_H

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

#endif

template<typename T>
bool equalityCheck(const string title, const T expected, const T received) {
  bool check = expected == received;
  cout << "=== "<< title << " ===" << endl;
  cout << "  expected:\t"<< expected << endl;
  cout << "  received:\t"<< received << endl;
  cout << "  check:\t"<< check << "\t\t\t" << (check ? "Success!" : "Error") << endl;
  cout << endl;
  return check;
}


int SemiregularSpline2D3DTest()
{
  using namespace ali_tpc_common::tpc_fast_transformation ;

  gRandom->SetSeed(0);
  UInt_t seed = gRandom->Integer(100000); // 605
  gRandom->SetSeed(seed);
  cout<<"Random seed: "<<seed<<" "<<gRandom->GetSeed()<<endl;

  for( int i=0; i<=PolynomDegree; i++ ){
    cu[i] = gRandom->Uniform(-1,1);
    cv[i] = gRandom->Uniform(-1,1);
  }


  const int numberOfRows = 6;
  const int numbersOfKnots[numberOfRows] = {6, 6, 6, 6, 6, 6};
  int knotAmountExp = 0; //Expected amount of knots in total
  for( int i=0; i < numberOfRows; i++) {
    knotAmountExp += numbersOfKnots[i];
  }
    
  SemiregularSpline2D3D spline;
  spline.construct( numberOfRows, numbersOfKnots );
  std::vector<bool> checker;
  
  
 
  int nKnotsTot = spline.getNumberOfKnots();
  checker.push_back(equalityCheck("Number of rows", numberOfRows, spline.getNumberOfRows()));
  checker.push_back(equalityCheck("Number of knots", knotAmountExp, nKnotsTot));  

  const RegularSpline1D &gridV = spline.getGridV();

  float *data0 = new float[ 6*nKnotsTot ];
  float *data = new float[ 6*nKnotsTot ];

  int nv = gridV.getNumberOfKnots();   //4

  checker.push_back(equalityCheck("Number of Knots in gridV", numberOfRows, nv));

  //loop through all rows (v-coordinate)
  for( int i=0; i<nv; i++ ){

    //get each v coordinate
    double v = gridV.indexToU(i);

    //get Spline for u-coordinate at that point
    const RegularSpline1D &gridU = spline.getGridU(i);

    //loop through all u-indexes
    for( int j=0; j<gridU.getNumberOfKnots(); j++ ){  
      //get the u coodrinate
      double u =  gridU.indexToU(j);
      int ind = spline.getDataIndex(j, i);
      data0[ind+0] = Fx(u,v);
      data0[ind+1] = Fy(u,v);
      data0[ind+2] = Fz(u,v);
      // just some random values
      data0[ind+0] = gRandom->Uniform(-1,1);//Gaus();
    }
  }



  for( int i=0; i<nKnotsTot*3; i++ ){
    data[i] = data0[i];
  }

  cout << "Start correcting edges..." << endl;
  spline.correctEdges( data );
  cout << "corrected edges!" << endl;

 

  TCanvas *canv = new TCanvas("cQA","2D splines  QA",1500,1500);
  canv->Draw();
  //canv->Divide(3,1);
  //canv->Update();

  TGraph2D *gknots = new TGraph2D();
  gknots->SetName("gknots");
  gknots->SetTitle("gknots");
  gknots->SetLineColor(kRed);
  gknots->SetMarkerSize(1.);
  gknots->SetMarkerStyle(8);
  gknots->SetMarkerColor(kRed);

  int gknotsN=0;
  TNtuple *knots = new TNtuple("knots","knots", "u:v:f");
  double diff=0;
  for( int i=0; i<gridV.getNumberOfKnots(); i++ ){
    const RegularSpline1D &gridU = spline.getGridU(i);
    for( int j=0; j<gridU.getNumberOfKnots(); j++ ){
      double v = gridV.indexToU(i);
      double u = gridU.indexToU(j);
      int ind = spline.getDataIndex(j, i);
      double fx0 = data0[ind+0];
      knots->Fill( u, v, fx0 );
      float x, y, z;
      spline.getSpline( data, u,v,x,y,z);    
      diff+=(fx0-x)*(fx0-x);
      gknots->SetPoint(gknotsN++,u,v,fx0);      
    }
  }
  
  checker.push_back(equalityCheck("inserted points", nKnotsTot, gknotsN));
  cout<<"mean diff at knots: "<<sqrt(diff)/gknotsN<<endl;

  knots->SetMarkerSize(1.);
  knots->SetMarkerStyle(8);
  knots->SetMarkerColor(kRed);

  TNtuple *nt = new TNtuple("nt","nt","u:v:f0:fs");


  TGraph2D *gf0 = new TGraph2D();
  gf0->SetName("gf0");
  gf0->SetTitle("gf0");
  gf0->SetLineColor(kRed);
  int gf0N=0;

  TGraph2D *gfs = new TGraph2D();
  gfs->SetName("gfs");
  gfs->SetTitle("gfs");
  gfs->SetLineColor(kBlue);
  int gfsN=0;

  TH1F *qaX = new TH1F("qaX","qaX [um]",1000,-100.,100.);
  TH1F *qaY = new TH1F("qaY","qaY [um]",1000,-10.,10.);
  TH1F *qaZ = new TH1F("qaZ","qaZ [um]",1000,-10.,10.);

  int iter=0;
  float stepu = 1.e-3;
  float stepv = 1.e-3;			       

  for( float u=-0.01; u<=1.01; u+=stepu ){
    for( float v=-0.01; v<=1.01; v+=stepv ){
      iter++;
      float x, y, z;
      spline.getSplineVec( data, u,v,x,y,z);
      if( u>=0 && v>=0 && u<=1 && v<=1 ){
	qaX->Fill( 1.e4*(x - Fx(u,v)) );
	qaY->Fill( 1.e4*(y - Fy(u,v)) );
	qaZ->Fill( 1.e4*(z - Fz(u,v)) );
	nt->Fill(u,v,Fx(u,v),x);
      }
      gf0->SetPoint(gf0N++,u,v,Fx(u,v));
      gfs->SetPoint(gfsN++,u,v,x);
    }
  }

  gStyle->SetPalette(1);
  
  //gknots->Draw("P");  
  gfs->Draw("surf");
  //gf0->Draw("surf,same");
  gknots->Draw("P,same");
  canv->Update();

  /*
    Specific drawing options can be used to paint a TGraph2D:

    "TRI"  : The Delaunay triangles are drawn using filled area.
    An hidden surface drawing technique is used. The surface is
    painted with the current fill area color. The edges of each
    triangles are painted with the current line color.
    "TRIW" : The Delaunay triangles are drawn as wire frame
    "TRI1" : The Delaunay triangles are painted with color levels. The edges
    of each triangles are painted with the current line color.
    "TRI2" : the Delaunay triangles are painted with color levels.
    "P"    : Draw a marker at each vertex
    "P0"   : Draw a circle at each vertex. Each circle background is white.
    "PCOL" : Draw a marker at each vertex. The color of each marker is
    defined according to its Z position.
    "CONT" : Draw contours
    "LINE" : Draw a 3D polyline
  */

  //Check if everything went good
  bool success = true;
  for(unsigned int i=0; i < checker.size(); i++) {
    success = success && checker[i];
    if(!success) {
      cout << "Something went wrong!" << endl;
      break;
    }
  }
  if( success ) {
    cout << "Everything worked fine" << endl;
  }


  return success;
}
