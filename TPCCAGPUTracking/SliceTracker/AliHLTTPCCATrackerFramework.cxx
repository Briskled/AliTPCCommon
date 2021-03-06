// @(#) $Id: AliHLTTPCCATracker.cxx 34611 2009-09-04 00:22:05Z sgorbuno $
// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de> *
//                  Ivan Kisel <kisel@kip.uni-heidelberg.de>                *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "AliHLTTPCCADef.h"
#include "AliHLTTPCCAGPUConfig.h"
#include "AliHLTTPCCATrackerFramework.h"
#include "AliHLTTPCCAGPUTracker.h"
#include "AliHLTTPCCATracker.h"
#include "AliTPCCommonMath.h"
#include "AliHLTTPCCAClusterData.h"
#include "AliGPUReconstruction.h"

#ifdef WIN32
#include <windows.h>
#include <winbase.h>
#else
#include <dlfcn.h>
#endif

#if defined(HLTCA_BUILD_O2_LIB) & defined(HLTCA_STANDALONE)
#undef HLTCA_STANDALONE //We disable the standalone application features for the O2 lib. This is a hack since the HLTCA_STANDALONE setting is ambigious... In this file it affects standalone application features, in the other files it means independence from aliroot
#endif

#ifdef HLTCA_HAVE_OPENMP
#include <omp.h>
#endif

ClassImp( AliHLTTPCCATrackerFramework )

int AliHLTTPCCATrackerFramework::InitGPU(int forceDeviceID)
{
	//Initialize GPU Tracker and determine if GPU available
	int retVal;
	if (!fGPULibAvailable)
	{
		CAGPUError("GPU Library not loaded\n");
		return(1);
	}
	if (fGPUTrackerAvailable && (retVal = ExitGPU())) return(retVal);
	retVal = fGPUTracker->InitGPU(-1, forceDeviceID);
	fUseGPUTracker = fGPUTrackerAvailable = retVal == 0;
	return(retVal);
}

int AliHLTTPCCATrackerFramework::ExitGPU()
{
	//Uninitialize GPU Tracker
	if (!fGPUTrackerAvailable) return(0);
	fUseGPUTracker = false;
	fGPUTrackerAvailable = false;
	return(fGPUTracker->ExitGPU());
}

void AliHLTTPCCATrackerFramework::SetGPUDebugLevel(int Level, std::ostream *OutFile, std::ostream *GPUOutFile)
{
	//Set Debug Level for GPU Tracker and also for CPU Tracker for comparison reasons
	fGPUTracker->SetDebugLevel(Level, GPUOutFile);
	fGPUDebugLevel = Level;
	for (int i = 0;i < fgkNSlices;i++)
	{
		fCPUTrackers[i].SetGPUDebugLevel(Level, OutFile);
	}
}

int AliHLTTPCCATrackerFramework::SetGPUTracker(bool enable)
{
	//Enable / disable GPU Tracker
	if (enable && !fGPUTrackerAvailable)
	{
		fUseGPUTracker = false;
		return(1);
	}
	fUseGPUTracker = enable;
	return(0);
}

GPUhd() void AliHLTTPCCATrackerFramework::SetOutputControl( AliHLTTPCCASliceOutput::outputControlStruct* val)
{
	//Set Output Control Pointers
	fOutputControl = val;
	fGPUTracker->SetOutputControl(val);
	for (int i = 0;i < fgkNSlices;i++)
	{
		fCPUTrackers[i].SetOutputControl(val);
	}
}

int AliHLTTPCCATrackerFramework::ProcessSlices(int firstSlice, int sliceCount, AliHLTTPCCAClusterData* pClusterData, AliHLTTPCCASliceOutput** pOutput)
{
	long long int totalNClusters = 0;
	for (int iSlice = 0;iSlice < CAMath::Min(sliceCount, fgkNSlices - firstSlice);iSlice++)
	{
		totalNClusters += pClusterData[iSlice].NumberOfClusters();
		if (totalNClusters >= ((long long int) 1) << (sizeof(int) * 8))
		{
			printf("Too many clusters for cluster indexing: %lld >= %lld\n", totalNClusters, ((long long int) 1) << (sizeof(int) * 8));
			return(1);
		}
	}
	int useGlobalTracking = fGlobalTracking;
	if (fGlobalTracking && (firstSlice || sliceCount != fgkNSlices))
	{
		CAGPUWarning("Global Tracking only available if all slices are processed!");
		useGlobalTracking = 0;
	}

	//Process sliceCount slices starting from firstslice, in is pClusterData array, out pOutput array
	if (fUseGPUTracker)
	{
		if (fGPUTracker->Reconstruct(pOutput, pClusterData, firstSlice, CAMath::Min(sliceCount, fgkNSlices - firstSlice))) return(1);
	}
	else
	{
		bool error = false;
#ifdef HLTCA_STANDALONE
		int nLocalTracks = 0, nGlobalTracks = 0, nOutputTracks = 0, nLocalHits = 0, nGlobalHits = 0;
#endif
#ifdef HLTCA_HAVE_OPENMP
		if (fOutputControl->fOutputPtr && omp_get_max_threads() > 1)
		{
			CAGPUError("fOutputPtr must not be used with OpenMP\n");
			return(1);
		}
#pragma omp parallel for
#endif
		for (int iSlice = 0;iSlice < CAMath::Min(sliceCount, fgkNSlices - firstSlice);iSlice++)
		{
			if (error) continue;
			if (fCPUTrackers[firstSlice + iSlice].ReadEvent(&pClusterData[iSlice]))
			{
				error = true;
				continue;
			}
			fCPUTrackers[firstSlice + iSlice].SetOutput(&pOutput[iSlice]);
			fCPUTrackers[firstSlice + iSlice].Reconstruct();
			fCPUTrackers[firstSlice + iSlice].CommonMemory()->fNLocalTracks = fCPUTrackers[firstSlice + iSlice].CommonMemory()->fNTracks;
			fCPUTrackers[firstSlice + iSlice].CommonMemory()->fNLocalTrackHits = fCPUTrackers[firstSlice + iSlice].CommonMemory()->fNTrackHits;
			if (!useGlobalTracking)
			{
				fCPUTrackers[firstSlice + iSlice].ReconstructOutput();
#ifdef HLTCA_STANDALONE
				nOutputTracks += (*fCPUTrackers[firstSlice + iSlice].Output())->NTracks();
				nLocalTracks += fCPUTrackers[firstSlice + iSlice].CommonMemory()->fNTracks;
#endif
				if (!fKeepData)
				{
					fCPUTrackers[firstSlice + iSlice].SetupCommonMemory();
				}
			}
		}
		if (error) return(1);

		if (useGlobalTracking)
		{
			for (int iSlice = 0;iSlice < CAMath::Min(sliceCount, fgkNSlices - firstSlice);iSlice++)
			{
				int sliceLeft = (iSlice + (fgkNSlices / 2 - 1)) % (fgkNSlices / 2);
				int sliceRight = (iSlice + 1) % (fgkNSlices / 2);
				if (iSlice >= fgkNSlices / 2)
				{
					sliceLeft += fgkNSlices / 2;
					sliceRight += fgkNSlices / 2;
				}
				fCPUTrackers[iSlice].PerformGlobalTracking(fCPUTrackers[sliceLeft], fCPUTrackers[sliceRight], fCPUTrackers[sliceLeft].NMaxTracks(), fCPUTrackers[sliceRight].NMaxTracks());
			}
			for (int iSlice = 0;iSlice < CAMath::Min(sliceCount, fgkNSlices - firstSlice);iSlice++)
			{
				fCPUTrackers[firstSlice + iSlice].ReconstructOutput();
#ifdef HLTCA_STANDALONE
				//printf("Slice %d - Tracks: Local %d Global %d - Hits: Local %d Global %d\n", iSlice, fCPUTrackers[iSlice].CommonMemory()->fNLocalTracks, fCPUTrackers[iSlice].CommonMemory()->fNTracks, fCPUTrackers[iSlice].CommonMemory()->fNLocalTrackHits, fCPUTrackers[iSlice].CommonMemory()->fNTrackHits);
				nLocalTracks += fCPUTrackers[iSlice].CommonMemory()->fNLocalTracks;
				nGlobalTracks += fCPUTrackers[iSlice].CommonMemory()->fNTracks;
				nLocalHits += fCPUTrackers[iSlice].CommonMemory()->fNLocalTrackHits;
				nGlobalHits += fCPUTrackers[iSlice].CommonMemory()->fNTrackHits;
				nOutputTracks += (*fCPUTrackers[iSlice].Output())->NTracks();
#endif
				if (!fKeepData)
				{
					fCPUTrackers[firstSlice + iSlice].SetupCommonMemory();
				}
			}
		}
		for (int iSlice = 0;iSlice < CAMath::Min(sliceCount, fgkNSlices - firstSlice);iSlice++)
		{
			if (fCPUTrackers[iSlice].GPUParameters()->fGPUError != 0)
			{
				const char* errorMsgs[] = HLTCA_GPU_ERROR_STRINGS;
				const char* errorMsg = (unsigned) fCPUTrackers[iSlice].GPUParameters()->fGPUError >= sizeof(errorMsgs) / sizeof(errorMsgs[0]) ? "UNKNOWN" : errorMsgs[fCPUTrackers[iSlice].GPUParameters()->fGPUError];
				printf("Error during tracking: %s\n", errorMsg);
				return(1);
			}
		}
#ifdef HLTCA_STANDALONE
		//printf("Slice Tracks Output %d: - Tracks: %d local, %d global -  Hits: %d local, %d global\n", nOutputTracks, nLocalTracks, nGlobalTracks, nLocalHits, nGlobalHits);
		/*for (int i = firstSlice;i < firstSlice + sliceCount;i++)
		{
			fCPUTrackers[i].DumpOutput(stdout);
		}*/
#endif
	}

	if (fGPUDebugLevel >= 6 && fUseGPUTracker)
	{
	    fUseGPUTracker = 0;
	    ProcessSlices(firstSlice, sliceCount, pClusterData, pOutput);
	    fUseGPUTracker = 1;
	}

	return(0);
}

double AliHLTTPCCATrackerFramework::GetTimer(int iSlice, int iTimer)
{
	return(fUseGPUTracker ? fGPUTracker->GetTimer(iSlice, iTimer) : fCPUTrackers[iSlice].GetTimer(iTimer));
}
void AliHLTTPCCATrackerFramework::ResetTimer(int iSlice, int iTimer)
{
	return(fUseGPUTracker ? fGPUTracker->ResetTimer(iSlice, iTimer) : fCPUTrackers[iSlice].ResetTimer(iTimer));
}

int AliHLTTPCCATrackerFramework::InitializeSliceParam(int iSlice, const AliGPUCAParam *param)
{
	//Initialize Tracker Parameters for a slice
	if (fGPUTrackerAvailable && fGPUTracker->InitializeSliceParam(iSlice, param)) return(1);
	fCPUTrackers[iSlice].Initialize(param, iSlice);
	return(0);
}

AliHLTTPCCATrackerFramework::AliHLTTPCCATrackerFramework(AliGPUReconstruction* rec) : fGPULibAvailable(false), fGPUTrackerAvailable(false), fUseGPUTracker(false), fGPUDebugLevel(0), fGPUTracker(NULL), fGPULib(NULL), fOutputControl( NULL ), fKeepData(false), fGlobalTracking(false)
{
	//Constructor
	if (rec) fGPUTracker = rec->GetTPCTracker();
	fGPULibAvailable = (fGPUTracker != NULL);

	if (fGPULibAvailable)
	{
		fUseGPUTracker = (fGPUTrackerAvailable = (fGPUTracker->InitGPU(-1, -1) == 0));
		if(fUseGPUTracker)
		{
		  CAGPUInfo("GPU Tracker Initialized and available in framework");
		}
		else
		{
		  CAGPUError("GPU Tracker NOT Initialized and NOT available in framework");
		}
	}
	else
	{
		fGPUTracker = new AliHLTTPCCAGPUTracker;
	}
}

AliHLTTPCCATrackerFramework::~AliHLTTPCCATrackerFramework()
{
}
