#include <dataProcessing/procDistrHisto.h>

ProcDistrHisto::ProcDistrHisto(ParticlesSet* set, std::string path) : 
    part_set(set), 
    xpath(path)
{        
  LOG(PROC_DISTR_HISTO, Instanciating ProcDistrHisto);
  store = part_set->getStore();
  // Parameters
  readParameters();
  hist_counter = 0;
  // Set ROOT options
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(100);
  // Create the canvas
  int number_projections = projections.size();
  if(number_projections % 2 == 1)
    {
      number_projections += 1;
    }
  number_projections /= 2;
  canvas = new TCanvas(title.c_str(),title.c_str(),2*350,number_projections*350);
  canvas->SetBorderMode(0);
  canvas->SetFillColor(kWhite);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(100);
  canvas->SetBottomMargin(0.13);
  canvas->Divide(static_cast<const Int_t>(2),static_cast<Int_t>(number_projections));
  LOG(PROC_DISTR_HISTO, ProcDistrHisto instanciated.);
}

ProcDistrHisto::~ProcDistrHisto()
{
    LOG(PROC_DISTR_HISTO, Destructor called.);
    // Delete all histograms in the list
    std::list<TH2D*>::iterator i;
    for(i=hist_list.begin();i!=hist_list.end();i++)
    {
        delete *i;
    }
    delete canvas;
}

void
ProcDistrHisto::readParameters()
{
    LOG(PROC_DISTR_HISTO, Reading parameters);
    box_size = XMLInputParser::getInstance().getFirstTextd("./flows/alessandroFlow/windowSize");
    LOG(PROC_DISTR_HISTO, Box size:);
    LOGV(PROC_DISTR_HISTO, d2s(box_size));
    file = XMLInputParser::getInstance().getAttributeTextFromPath(xpath, "file");
	title = XMLInputParser::getInstance().getAttributeTextFromPath(xpath, "title");
	std::string p = xpath + "/projection";
	projections = XMLInputParser::getInstance().getTextElementsFromPath(p);             
}

void
ProcDistrHisto::save()
{
    LOG(PROC_DISTR_HISTO, Saving the canvas in the .pdf file.);
    // Save the canvas in a pdf file
    canvas->SaveAs(file.c_str(),"pdf");
}

void
ProcDistrHisto::draw()
{
    std::vector<std::string>::iterator v;
    for(v=projections.begin();v != projections.end(); v++)
    {
        if((*v) == "XXP")
        {
            LOG(PROC_DISTR_HISTO, Creating histogram for X - XP);
            createHistogram("X","XP",X_PLANE,XP_PLANE);
        }
        else if((*v) == "XY")
        {
            if(part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
            {
                LOG(PROC_DISTR_HISTO, Creating histogram for X - Y);
                createHistogram("X","Y",X_PLANE,Y_PLANE);
            }
            else
            {
                WARN(PROC_DISTR_HISTO, Invalid histogram requested for X - Y);
            }
        }
        else if((*v) == "XYP")
        {
            if(part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() ==SIX_DIMENSIONS)
            {
                LOG(PROC_DISTR_HISTO, Creating histogram for X - YP);
                createHistogram("X","YP",X_PLANE,YP_PLANE);
            }
            else
            {
                WARN(PROC_DISTR_HISTO, Invalid histogram requested for X - YP);
            }
        }
        else if((*v) == "YXP")
        {
            if(part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
            {
                LOG(PROC_DISTR_HISTO, Creating histogram for Y - XP);
                createHistogram("Y","XP",Y_PLANE,XP_PLANE);
            }
            else
            {
                WARN(PROC_DISTR_HISTO, Invalid histogram requested for Y - XP);
            }
        }
        else if((*v) == "YYP")
        {
            if(part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
            {
                LOG(PROC_DISTR_HISTO, Creating histogram for Y - YP);
                createHistogram("Y","YP",Y_PLANE,YP_PLANE);
            }
            else
            {
                WARN(PROC_DISTR_HISTO, Invalid histogram requested for Y - YP);
            }
        }
        else if((*v) == "XPYP")
        {
            if(part_set->getDims() == FOUR_DIMENSIONS || part_set->getDims() == SIX_DIMENSIONS)
            {
                LOG(PROC_DISTR_HISTO, Creating histogram for XP - YP);
                createHistogram("XP","YP",XP_PLANE,YP_PLANE);
            }
            else
            {
                WARN(PROC_DISTR_HISTO, Invalid histogram requested for XP - YP);
            }
        }
        else if((*v) == "ZZP")
        {
            if(part_set->getDims() == SIX_DIMENSIONS)
            {
                LOG(PROC_DISTR_HISTO, Creating histogram for Z - ZP);
                createHistogram("Z","ZP",4,5);
            }
            else
            {
                WARN(PROC_DISTR_HISTO, Invalid histogram requested for Z - ZP);
            }
        }
    }
    std::list<TH2D*>::iterator i;
    for(i=hist_list.begin();i!=hist_list.end();i++)
    {
        #ifdef DEBUG_FLAG
            Log::getInstance().log("PROC_DISTR_HISTO", "Adding histogram in pad", "id: ", hist_counter+1);
        #endif
	    canvas->cd(hist_counter+1)->SetGrid();
	   (*i)->Draw("contz");
        hist_counter++;
    }
}

void
ProcDistrHisto::createHistogram(std::string axis1, std::string axis2, unsigned short int a1, unsigned short int a2)
{ 
    TH2D* hist = new TH2D((axis1+axis2).c_str(),title.c_str(),300,-static_cast<double>(box_size),static_cast<double>(box_size),300,-static_cast<double>(box_size),static_cast<double>(box_size));
    hist_list.push_back(hist);
    hist->GetXaxis()->SetTitle(axis1.c_str());
    hist->GetXaxis()->CenterTitle();
    hist->GetXaxis()->SetTitleSize(static_cast<Float_t>(0.04));
    hist->GetYaxis()->SetTitle(axis2.c_str());
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleSize(static_cast<Float_t>(0.04));
    
	for(int i = 0; i < part_set->getParticlesNumber(); ++i)
	{
        hist->Fill(static_cast<double>((*store)(i,static_cast<int>(a1))),static_cast<double>((*store)(i,static_cast<int>(a2))));
	}
    
    LOG(PROC_DISTR_HISTO, Histogram created);
}

