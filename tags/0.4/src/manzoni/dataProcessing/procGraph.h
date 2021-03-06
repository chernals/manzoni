#ifndef PROC_GRAPH_H
#define PROC_GRAPH_H

#include <list>  
#include <vector>
#include <string>

#include <XMLInputParser.h>
#include <utils/Logger.h>
#include <utils/TypeConversions.h>
#include <root.h>

class ProcGraph
{
	private:            
		//
		// Data
		//
	    // Graph options
	    const std::string file;
	    const std::string title;
	    TCanvas* canvas;
	    // Graph data
	    double* coord;
	    const int number;
	    std::list<TGraph*> graphs_list;
	    // Misc
	    bool is_first;

	public:
		ProcGraph(int, std::string, std::string);
		~ProcGraph();
		void draw(double*, std::string, bool last=false);
};

#endif               
