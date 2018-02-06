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

#include <utils/XMLInputParser.h>

XMLParser::XMLParser()
{           
    LOG(XML PARSER, Instanciating XMLInputParser);
       
    try
    {
        LOG(XML PARSER, Initializing Xerces);
        // Initialization of Xerces
        xercesc::XMLPlatformUtils::Initialize(); // Needed before any Xerces op
        parser = new xercesc::XercesDOMParser();
        parser->setValidationScheme(xercesc::XercesDOMParser::Val_Never);
        parser->setDoNamespaces(false);
        parser->setDoSchema(false);
        parser->setLoadExternalDTD(false);
        
        // Input file parameter
 	    const char* input_parameter = Arg::getInstance().get_argument('i');
	    LOG(XML PARSER, Trying to load the XML input file: );
	    LOGV(XML PARSER, input_parameter);
	    
	    // Check if the input path is OK
	    struct stat stat_file_info;
	    int ret_stat = stat(input_parameter, &stat_file_info);
	    if(ret_stat == 0)
	    {
	        LOG(XML PARSER, Input path OK.);
	    }
	    else if(ret_stat == ENOENT)
	    {
	        ERROR(XML PARSER, The input path does not exist or is empty !);
	        throw("The input path does not exist or is empty !");
	    }
	    else if(ret_stat == ENOTDIR)
	    {
	        ERROR(XML PARSER, The input path contains a component that is not a directory !);
	        throw("The input path contains a component that is not a directory !");
	    }
	    else if(ret_stat == ELOOP)
	    {
	        ERROR(XML PARSER, The input path contains too many symbolic links !);
	        LOG(XML PARSER, I never expected this exception to be encountered. Contact me <cedric.hernalsteens@cern.ch> I will pay you a beer !);
	        throw("The input path contains too many symbolic links !");
	    }
	    else if(ret_stat == EACCES)
	    {
	        ERROR(XML PARSER, You do not have the access rights to access the input path !);
	        throw("You don't have the access rights to access the input path !");
	    }
	    else if(ret_stat == ENAMETOOLONG)
	    {
	        ERROR(XML PARSER, The input path cannot be read !);
	        throw("The input path cannot be read !");
	    }
	    else
	    {
	    	ERROR(XML PARSER, Input path not found !);
	        throw("Input path not found !");
	    }
	    
	    // Load and parser the input file
	    parser->parse(input_parameter);
	    LOG(XML PARSER, Input file successfully loaded);
	    xmlDoc = parser->getDocument();
	    rootElement = xmlDoc->getDocumentElement();
	    
	    // Initialize the serializer
	    XMLCh tempStr[3] = {xercesc::chLatin_L,xercesc::chLatin_S,xercesc::chNull};
        xercesc::DOMImplementation* impl = xercesc::DOMImplementationRegistry::getDOMImplementation(tempStr);
        serializer = (reinterpret_cast<xercesc::DOMImplementationLS*>(impl))->createLSSerializer();
        outputDesc = (reinterpret_cast<xercesc::DOMImplementationLS*>(impl))->createLSOutput();
        outputDesc->setEncoding(0);
        target = new xercesc::StdOutFormatTarget();
        outputDesc->setByteStream(target);
    
	    	// Store the simulation name
	    	XMLCh* attr_sim_name = xercesc::XMLString::transcode("name");
        const XMLCh* sim_name = rootElement->getAttribute(attr_sim_name);
        this->simulation_name = xercesc::XMLString::transcode(sim_name);

				// Store the flow name
				XMLCh* attr_flow_name = xercesc::XMLString::transcode("flow");
				const XMLCh* flowname = rootElement->getAttribute(attr_flow_name);
				this->flow_name = xercesc::XMLString::transcode(flowname);
    } // Try
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Xerces initialization error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
}

XMLParser::~XMLParser()
{
    try
    {
        LOG(XML PARSER, Destructor called);
        delete target;
        delete parser;
        xercesc::XMLPlatformUtils::Terminate();
        // Make sure that no Xerces related destructor is called after this point !
    } // Trys
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Xerces teardown error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
}

std::string
XMLParser::getAttributeTextFromPath(std::string path, std::string attr)
{
    try
    {
        std::string attribute = "";
	    XMLCh* attr_name = xercesc::XMLString::transcode(attr.c_str());
	    
	    if(isFoundFromPath(path))
	    {
	        xercesc::DOMNode* node = getFirstNodeFromPath(path);
	        xercesc::DOMNamedNodeMap* attr_map = NULL;
	        attr_map = node->getAttributes();
	        if(attr_map == NULL)
	        {
	            throw("Problem while looking for attributes (NULL attr map)!");
	        } 
	        int returned_length = static_cast<int>(attr_map->getLength());
	        if(returned_length == 0)
	        {
	            throw("Problem while looking for attributes !");
	        }
	        const XMLCh* attr_ch = attr_map->getNamedItem(attr_name)->getNodeValue();
	        attribute = xercesc::XMLString::transcode(attr_ch);
	        #ifdef XML_PARSER_DEBUG
	            LOG(XML PARSER, The following XML attribute has been retrieved by XPath:);
	            LOGV(XML PARSER, attribute); 
	        #endif
	    }
	    else
	    {
	        throw("Node not found (AttributeText) !");
	    }
	    return attribute;
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Attribute from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch    
    catch(const char* excep)
	{
		ERRORV(XML PARSER CATCH, excep);
		ERROR(XML PARSER CATCH, Forwarding...);
		throw("Exception catched in XML PARSER");
	} // Catch
}

double
XMLParser::getAttribute(std::string path, std::string attr)
{
    return s2n<double>(getAttributeTextFromPath(path, attr));
}

double
XMLParser::getAttribute(const char path[], const char attr[])
{
    return s2n<double>(getAttributeTextFromPath(path, attr));
}

xercesc::DOMNode*
XMLParser::getFirstNodeFromPath(std::string path)
{
    xercesc::DOMNode* retrieved_node = NULL;
    try
    {
        XMLCh* xpathExpr = xercesc::XMLString::transcode(path.c_str());
        // Create the resolver and solve the XPath expression
        xercesc::DOMXPathNSResolver* resolver = xmlDoc->createNSResolver(rootElement);
        xercesc::DOMXPathResult* result = xmlDoc->evaluate(xpathExpr,rootElement,resolver,xercesc::DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,NULL);
    
        // Process the resulting nodes list
        XMLSize_t nLength = result->getSnapshotLength();
        if(nLength == 0)
        {
            WARN(XML PARSER, XPath expression not found:);
            WARNV(XML PARSER, path);
            throw("Node not found !");
        }
        else
        {
            result->snapshotItem(0);
            #ifdef XML_PARSER_DEBUG
                LOG(XML PARSER, The following XML snippet has been retrieved by XPath:);
                serializer->write(result->getNodeValue(),outputDesc);
                std::cout << std::endl;
            #endif
            retrieved_node = result->getNodeValue();             
        }       
        // Clean and return
        delete resolver;
        return retrieved_node;
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Node from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
    catch(xercesc::DOMXPathException& e)
    {
        ERROR(XML PARSER CATCH, Node from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch 
    catch(const char* excep)
	{
		ERRORV(XML PARSER CATCH, excep);
		ERROR(XML PARSER CATCH, Forwarding...);
		throw("Exception catched in XML PARSER");
	} // Catch
}

std::vector<xercesc::DOMNode*>
XMLParser::getNodesFromPath(std::string path)
{
    std::vector<xercesc::DOMNode*> retrieved_nodes;
    try
    {
        XMLCh* xpathExpr = xercesc::XMLString::transcode(path.c_str());
        // Create the resolver and solve the XPath expression
        xercesc::DOMXPathNSResolver* resolver = xmlDoc->createNSResolver(rootElement);
        xercesc::DOMXPathResult* result = xmlDoc->evaluate(xpathExpr,rootElement,resolver,xercesc::DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,NULL);
    
        // Process the resulting nodes list
        XMLSize_t nLength = result->getSnapshotLength();
        if(nLength == 0)
        {
            WARN(XML PARSER, XPath expression not found:);
            WARNV(XML PARSER, path);
            throw("Node not found !");
        }
        else
        {
            for(XMLSize_t i =0; i < nLength; i++)
            {
                result->snapshotItem(i);
                #ifdef XML_PARSER_DEBUG
                    LOG(XML PARSER, The following XML snippet has been retrieved by XPath:);
                    serializer->write(result->getNodeValue(),outputDesc);
                    std::cout << std::endl;
                #endif
                retrieved_nodes.push_back(result->getNodeValue());
            }
         }       
        // Clean and return
        delete resolver;
        return retrieved_nodes;
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Node from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
    catch(xercesc::DOMXPathException& e)
    {
        ERROR(XML PARSER CATCH, Node from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch 
    catch(const char* excep)
	{
		ERRORV(XML PARSER CATCH, excep);
		ERROR(XML PARSER CATCH, Forwarding...);
		throw("Exception catched in XML PARSER");
	} // Catch
}

std::string
XMLParser::getFirstTextElementFromPath(std::string path)
{
    bool isFound = false;
    return getFirstTextElementFromPath(path, isFound);
}

double
XMLParser::getFirstTextd(std::string path)
{
    return s2n<double>(getFirstTextElementFromPath(path));
}

double
XMLParser::getFirstTextd(const char path[])
{
    return s2n<double>(getFirstTextElementFromPath(path));
}

int
XMLParser::getFirstTexti(std::string path)
{
    return s2n<int>(getFirstTextElementFromPath(path));
}

int
XMLParser::getFirstTexti(const char path[])
{
    return s2n<int>(getFirstTextElementFromPath(path));
}

std::string
XMLParser::getFirstTextElementFromPath(std::string path, bool& isFound)
{
    std::string retrieved_string = "";
    xercesc::DOMNode* node = getFirstNodeFromPath(path);
    isFound = false;
    try
    {
        isFound = true;
        retrieved_string = xercesc::XMLString::transcode(node->getFirstChild()->getNodeValue());        
        return ML::getInstance().execute(retrieved_string);
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, First text element from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        Log::getInstance().error("XML PARSER CATCH", message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
    catch(xercesc::DOMXPathException& e)
    {
        ERROR(XML PARSER CATCH, First text element from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch 
}

std::vector<std::string>
XMLParser::getTextElementsFromNode(xercesc::DOMNode* node, std::string path)
{
    std::vector<std::string> retrievedStringsVector;
    try
    {
        XMLCh* xpathExpr = xercesc::XMLString::transcode(path.c_str());
        // Create the resolver and solve the XPath expression
        xercesc::DOMXPathNSResolver* resolver = xmlDoc->createNSResolver(node);
        xercesc::DOMXPathResult* result = xmlDoc->evaluate(xpathExpr,node,resolver,xercesc::DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,NULL);
    
        // Process the resulting nodes list
        XMLSize_t nLength = result->getSnapshotLength();
        if(nLength == 0)
        {
            WARN(XML PARSER, XPath expression not found);
        }
        for(XMLSize_t i =0;i<nLength;i++)
        {
            result->snapshotItem(i);
            #ifdef XML_PARSER_DEBUG
                LOG(XML PARSER, The following XML snippet has been retrieved by XPath:);
                serializer->write(result->getNodeValue(), outputDesc);
                std::cout << std::endl;
            #endif
            retrievedStringsVector.push_back(xercesc::XMLString::transcode(result->getNodeValue()->getFirstChild()->getNodeValue()));
        }
        // Clean and return
        delete resolver;
        return retrievedStringsVector;
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Text elements from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERROR(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
    catch(xercesc::DOMXPathException& e)
    {
        ERROR(XML PARSER CATCH, Text elements from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch 
}

std::vector<std::string>
XMLParser::getTextElementsFromPath(std::string path)
{
    std::vector<std::string> retrievedStringsVector;
    try
    {
        XMLCh* xpathExpr = xercesc::XMLString::transcode(path.c_str());
        // Create the resolver and solve the XPath expression
        xercesc::DOMXPathNSResolver* resolver = xmlDoc->createNSResolver(rootElement);
        xercesc::DOMXPathResult* result = xmlDoc->evaluate(xpathExpr,rootElement,resolver,xercesc::DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,NULL);
    
        // Process the resulting nodes list
        XMLSize_t nLength = result->getSnapshotLength();
        if(nLength == 0)
        {
            WARN(XML PARSER, XPath expression not found);
        }
        for(XMLSize_t i =0;i<nLength;i++)
        {
            result->snapshotItem(i);
            #ifdef XML_PARSER_DEBUG
                LOG(XML PARSER, The following XML snippet has been retrieved by XPath:);
                serializer->write(result->getNodeValue(), outputDesc);
                std::cout << std::endl;
            #endif
            retrievedStringsVector.push_back(xercesc::XMLString::transcode(result->getNodeValue()->getFirstChild()->getNodeValue()));
        }
        // Clean and return
        delete resolver;
        return retrievedStringsVector;
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, Text elements from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERROR(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
    catch(xercesc::DOMXPathException& e)
    {
        ERROR(XML PARSER CATCH, Text elements from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch 
}

unsigned short int
XMLParser::howManyFoundFromPath(std::string path)
{
    unsigned short int founds = 0;
    try
    {
        XMLCh* xpathExpr = xercesc::XMLString::transcode(path.c_str());
        // Create the resolver and solve the XPath expression
        xercesc::DOMXPathNSResolver* resolver = xmlDoc->createNSResolver(rootElement);
        xercesc::DOMXPathResult* result = xmlDoc->evaluate(xpathExpr,rootElement,resolver,xercesc::DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,NULL);
    
        // Process the resulting nodes list
        XMLSize_t nLength = result->getSnapshotLength();
        founds = static_cast<unsigned short int>(nLength);
        
        // Log it
        #ifdef XML_PARSER_DEBUG
          char* string = new char[BUFFER_SIZE];
	        sprintf(string, "Path checked by XPath: \"%s\"; found: %ld", path.c_str(),(long) founds);
	        LOGV(XML PARSER, string);
	        delete string;
	    	#endif
        
        // Clean and return
        delete resolver;
        return founds;
    }
    catch(xercesc::XMLException& e)
    {
        ERROR(XML PARSER CATCH, How many from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERROR(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch
    catch(xercesc::DOMXPathException& e)
    {
        ERROR(XML PARSER CATCH, How many from path error:);
        char* message = xercesc::XMLString::transcode(e.getMessage());
        ERRORV(XML PARSER CATCH, message);
        xercesc::XMLString::release(&message);
        ERROR(XML PARSER CATCH, Forwarding...);
        throw("Exception catched in XML PARSER");
    } // Catch 
}

bool
XMLParser::isFoundFromPath(std::string path)
{
    bool isFound = false;
    unsigned short int founds = howManyFoundFromPath(path);
    if(founds > 0)
    {
        isFound = true;
    }
    return isFound;
}

std::string
XMLParser::getSimulationName()
{
    return simulation_name;
}

std::string
XMLParser::getFlowName()
{
    return flow_name;
}
