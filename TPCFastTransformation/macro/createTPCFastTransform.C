// $Id$
/**
 * @file makeTPCFastTransform.C
 * @brief A macro to create TPCFastTransform object
 *
 * <pre>
 * Usage:
 *
 * aliroot
 * //.x initTPCcalibration.C("alien://Folder=/alice/data/2015/OCDB",246984,1)
 * .x initTPCcalibration.C("~/alice/events/OCDBsim.root",246984,1)
 *  gSystem->Load("libAliTPCFastTransformation")
 * .L createTPCFastTransform.C++
 * ali_tpc_common::tpc_fast_transformation::TPCFastTransform fastTransform;
 * createTPCFastTransform(fastTransform);
 *
 * </pre>
 *
 * @author sergey gorbunov
 *
 */


#include "AliTPCcalibDB.h"
#include "Riostream.h"
#include "TStopwatch.h"

#include "TPCFastTransform.h"
#include "TPCFastTransformManager.h"
#include "TPCFastTransformQA.h"
#include <iostream>

using namespace std;
using namespace ali_tpc_common::tpc_fast_transformation;

int createTPCFastTransform(TPCFastTransform &fastTransform) {

  AliTPCcalibDB* tpcCalib=AliTPCcalibDB::Instance();
  if(!tpcCalib){
    cerr << "AliTPCcalibDB does not exist"<<endl;
    return -1;
    }
  AliTPCTransform *origTransform = tpcCalib->GetTransform();
  UInt_t timeStamp = origTransform->GetCurrentTimeStamp();

  /*TPCFastTransformManager manager;

  TStopwatch timer;
  timer.Start();  
  
  int err = manager.create( fastTransform, origTransform, timeStamp );

  timer.Stop();
  
  cout<<"\n\n Initialisation: "<<timer.CpuTime()<<" / "<<timer.RealTime()<<" sec.\n\n"<<endl;

  if( err!=0 ){
    cerr << "Cannot create fast transformation object from AliTPCcalibDB, TPCFastTransformManager returns  "<<err<<endl;
    return -1;
    } */

  // qa

  /*AliTPCChebCorr* chebCorr = origTransform->GetCorrMapCacheRef();
  chebCorr->Print("p");


  float y2x=50, z=100;
  cout << y2x << " " << z << endl;
  for(int sector=0; sector < 72; sector++) {
    cout << "sector: " << sector << endl;
    const AliCheb2DStack* stack = chebCorr->GetParam(sector, y2x, z);
    cout<<"class name: "<<stack->ClassName()<<endl;
    for( int iTPCrow=0; iTPCrow<159; iTPCrow++){
      unsigned int rows = stack->getNRows()[iTPCrow];
      unsigned int cols = stack->getNCols()[iTPCrow];
      //unsigned int rows = static_cast<unsigned int>(stack->getNRows());
      //unsigned int cols = static_cast<unsigned int>(stack->getNCols());
      cout << "rows: " << rows << "\t  cols: " << cols << endl;
    }
    }*/
  

  ali_tpc_common::tpc_fast_transformation::TPCFastTransformQA qa;
  qa.doQA( timeStamp );
  qa.doQASRS( timeStamp );
  
  return 0;  
}

