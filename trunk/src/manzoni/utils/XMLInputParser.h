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

#ifndef XML_PARSER_H
#define XML_PARSER_H

#include <sys/stat.h>

#include <string>
#include <vector>

#include <xerces.h>
#include <types/types.h>
#include <utils/TypeConversions.h>
#include <utils/Singleton.h>
#include <utils/Arguments.h>
#include <utils/MathLink.h>

class XMLParser
{
	private:
	    // For Xerces
    	xercesc::XercesDOMParser* parser;
    	xercesc::DOMDocument* xmlDoc;
    	xercesc::DOMElement* rootElement;
    	xercesc::DOMLSSerializer* serializer;
    	xercesc::DOMLSOutput* outputDesc;
    	xercesc::XMLFormatTarget* target;
        // Misc
    	std::string simulation_name;
			std::string flow_name;

	public:
		XMLParser();   
		~XMLParser();
		bool isFoundFromPath(std::string);
		unsigned short int howManyFoundFromPath(std::string);
		//
		// Getters --- Setters
		//
		// Attributes
		std::string getAttributeTextFromPath(std::string, std::string);
		double getAttribute(std::string, std::string);
		double getAttribute(const char[], const char[]);
		// Text element
		std::string getFirstTextElementFromPath(std::string);
		double getFirstTextd(std::string);
		double getFirstTextd(const char[]);
		int getFirstTexti(std::string);
		int getFirstTexti(const char[]);
		xercesc::DOMNode* getFirstNodeFromPath(std::string);
		std::vector<xercesc::DOMNode*> getNodesFromPath(std::string);
		std::string getFirstTextElementFromPath(std::string, bool& isFound);
		std::vector<std::string> getTextElementsFromPath(std::string);
		std::vector<std::string> getTextElementsFromNode(xercesc::DOMNode*, std::string);
		std::string getSimulationName();
		std::string getFlowName();
};

typedef Singleton<XMLParser> XMLInputParser;

#endif

