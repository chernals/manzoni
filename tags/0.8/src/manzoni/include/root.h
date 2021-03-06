#ifndef ROOT_H
#define ROOT_H

#define GCC_VERSION (__GNUC__ * 10000 \
		    +__GNUC_MINOR__ * 100 \
		    +__GNUC_PATCHLEVEL__)

#if GCC_VERSION >= 40200
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif
#pragma GCC diagnostic ignored "-Wno-long-long"
    
    #include <TCanvas.h>
    #include <TROOT.h>
    #include <TStyle.h>
    #include <TH2F.h>
    #include <TProfile.h>
    #include <TGraph.h>
    
#pragma GCC diagnostic warning "-Wno-long-long"
#if GCC_VERSION >= 40200
#pragma GCC diagnostic warning "-Wold-style-cast"
#pragma GCC diagnostic warning "-Wconversion"
#endif

#endif
