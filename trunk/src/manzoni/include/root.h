/**
* Manzoni
*
* Author: Cedric Hernaslteens <cedric.hernalsteens@cern.ch>
* 
* European Organization for Nuclear Research
*
* Copyright (c) 2010+ CERN. All rights reserved.
*
**/

#ifndef ROOT_H
#define ROOT_H

#define GCC_VERSION (__GNUC__ * 10000 \
		    +__GNUC_MINOR__ * 100 \
		    +__GNUC_PATCHLEVEL__)

#if GCC_VERSION >= 40200
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif    
    #include <TCanvas.h>
    #include <TROOT.h>
    #include <TStyle.h>
    #include <TH2F.h>
    #include <TProfile.h>
    #include <TGraph.h>
#if GCC_VERSION >= 40200
#pragma GCC diagnostic warning "-Wold-style-cast"
#pragma GCC diagnostic warning "-Wconversion"
#endif

#endif
