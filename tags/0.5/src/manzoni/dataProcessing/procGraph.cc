#include <dataProcessing/procGraph.h>

ProcGraph::ProcGraph(int n, std::string f, std::string t) :
    file(f),
    title(t),
    number(n)
{
    LOG(PROC_GRAPH, Instanciating ProcGraph...);
    
    // Variables
    is_first = true;
    
    // Set ROOT options and create the canvas
    canvas = new TCanvas(title.c_str(),title.c_str(),300,300);
    canvas->SetLeftMargin(static_cast<Float_t>(0.13));
    canvas->SetRightMargin(static_cast<Float_t>(0.13));
    canvas->SetBottomMargin(static_cast<Float_t>(0.13));
    
    // Create the coordinates table     
    coord = new double[number];
    for(int i=0; i<number;i++)
    {
        coord[i] = i;
    }

    LOG(PROC_GRAPH, ProcGraph instanciated.);
}

ProcGraph::~ProcGraph()
{
    LOG(PROC_GRAPH, Destructor called.);
    // Delete all graphs in the list
    std::list<TGraph*>::iterator i;
    for(i=graphs_list.begin(); i != graphs_list.end(); i++)
    {
        delete *i;
    }
    delete canvas;
}

void
ProcGraph::draw(double* data, std::string title, bool last)
{
    LOG(PROC_GRAPH, Drawing the graph...);
    TGraph* graph = new TGraph(number, coord, data);        
    graph->SetTitle(title.c_str());
    graphs_list.push_back(graph);
    graph->Draw("APL");    
    canvas->Update();
    std::string tmp_string;
    if(is_first && !last)
    {
      tmp_string = file + ".pdf("; // The trick is in this (
        canvas->SaveAs(tmp_string.c_str());
    }
    else if(!last)
    {
        tmp_string = file + ".pdf";
        canvas->SaveAs(tmp_string.c_str());
	LOG(PROC GRAPH, Serializing the Root object);
	tmp_string = file + ".C";
	canvas->SaveAs(tmp_string.c_str());
    }
    else
    {   
      tmp_string = file + ".pdf)"; // The trick is in this )
      canvas->SaveAs(tmp_string.c_str());
    }
    LOG(PROC_GRAPH, Graph drawed.);
}
