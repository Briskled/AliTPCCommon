//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTTPCCASTARTHITSFINDER_H
#define ALIHLTTPCCASTARTHITSFINDER_H

#include "AliHLTTPCCADef.h"
#include "AliHLTTPCCAHitId.h"

MEM_CLASS_PRE() class AliHLTTPCCATracker;

/**
 * @class AliHLTTPCCAStartHitsFinder
 *
 */
class AliHLTTPCCAStartHitsFinder
{
  public:
    MEM_CLASS_PRE() class AliHLTTPCCASharedMemory
    {
        friend class AliHLTTPCCAStartHitsFinder;
      public:
#if !defined(HLTCA_GPUCODE)
        AliHLTTPCCASharedMemory()
            : fIRow( 0 ), fNHits( 0 ), fNRowStartHits( 0 ) {
}

        AliHLTTPCCASharedMemory( const AliHLTTPCCASharedMemory& /*dummy*/ )
	  : fIRow( 0 ), fNHits( 0 ), fNRowStartHits( 0 ) {
	}
        AliHLTTPCCASharedMemory& operator=( const AliHLTTPCCASharedMemory& /*dummy*/ ) { return *this; }
#endif //!HLTCA_GPUCODE

      protected:
        int fIRow; // row index
        int fNHits; // n hits in the row
        int fNRowStartHits; //start hits found in the row
    };

    GPUd() static int NThreadSyncPoints() { return 2; }

    GPUd() static void Thread( int nBlocks, int nThreads, int iBlock, int iThread, int iSync,
                               MEM_LOCAL(GPUsharedref() AliHLTTPCCASharedMemory) &smem,  MEM_CONSTANT(GPUconstant() AliHLTTPCCATracker) &tracker );
};


#endif //ALIHLTTPCCASTARTHITSFINDER_H
