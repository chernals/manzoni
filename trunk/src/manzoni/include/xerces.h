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

#ifndef XERCES_H
#define XERCES_H

#define GCC_VERSION (__GNUC__ * 10000 \
		    +__GNUC_MINOR__ * 100 \
		    +__GNUC_PATCHLEVEL__)

#if GCC_VERSION >= 40200
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif
    
	#include <xercesc/dom/DOM.hpp>
	#include <xercesc/dom/DOMDocument.hpp>
	#include <xercesc/dom/DOMDocumentType.hpp>
	#include <xercesc/dom/DOMElement.hpp>
	#include <xercesc/dom/DOMImplementation.hpp>
	#include <xercesc/dom/DOMImplementationLS.hpp>
	#include <xercesc/dom/DOMNodeIterator.hpp>
	#include <xercesc/dom/DOMNodeList.hpp>
	#include <xercesc/dom/DOMText.hpp>
	#include <xercesc/parsers/XercesDOMParser.hpp>
	#include <xercesc/util/XMLUni.hpp>
	#include <xercesc/framework/StdOutFormatTarget.hpp>

    
#if GCC_VERSION >= 40200
#pragma GCC diagnostic warning "-Woverloaded-virtual"
#pragma GCC diagnostic warning "-Wctor-dtor-privacy"
#pragma GCC diagnostic warning "-Wold-style-cast"
#pragma GCC diagnostic warning "-Wconversion"
#pragma GCC diagnostic warning "-Wall"
#endif

#endif
