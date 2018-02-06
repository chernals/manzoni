#include <dataProcessing/procDistrProfile.h>

ProcDistrProfile::ProcDistrProfile(ParticlesDistribution* part_set, std::string xpath) : 
    xpath_(xpath)
{        
	LOG(PROC DISTR PROFILE, Instanciating ProcDistrProfile);
	// Parameters
	_part_set = part_set; 
    store = part_set->getStore();
    readParameters();
    _profile_counter = 0;
    // Set ROOT options
    gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(500);
    // Create the canvas
    unsigned int number_profiles = static_cast<unsigned short int>(profiles_.size());
    if(number_profiles % 2 == 1)
    {
        number_profiles += 1U;
    }
    number_profiles /= 2;
    _canvas = new TCanvas(_title.c_str(),_title.c_str(),2*350,number_profiles*350);
    _canvas->SetLeftMargin(static_cast<Float_t>(0.13));
    _canvas->SetRightMargin(static_cast<Float_t>(0.13));
    _canvas->SetBottomMargin(static_cast<Float_t>(0.13));
    _canvas->Divide(2,number_profiles);
    LOG(PROC DISTR PROFILE, ProcDistrProfile instanciated.);
}

ProcDistrProfile::~ProcDistrProfile()
{
    LOG(PROC DISTR PROFILE, Destructor called.);
    // Delete all histograms in the list
    std::list<TH1D*>::iterator i;
    for(i=_profiles_list.begin();i!=_profiles_list.end();i++)
    {
        delete *i;
    }
    delete _canvas;
}

void
ProcDistrProfile::readParameters()
{
    LOG(PROC DISTR PROFILE, Reading parameters);
    box_size_ = XMLInputParser::getInstance().getFirstTextd("./flows/windowSize");
    _file = XMLInputParser::getInstance().getAttributeTextFromPath(xpath_, "file");
	_title = XMLInputParser::getInstance().getAttributeTextFromPath(xpath_, "title");
	std::string xpath = xpath_ + "/profile";
	profiles_ = XMLInputParser::getInstance().getTextElementsFromPath(xpath);             
}

void
ProcDistrProfile::save()
{
    LOG(PROC DISTR PROFILE, Saving the canvas in the .pdf file.);
    // Save the canvas in a pdf file
    _canvas->SaveAs(_file.c_str(),"pdf");
}

void
ProcDistrProfile::draw()
{
    std::vector<std::string>::iterator v;
    for(v=profiles_.begin();v != profiles_.end(); v++)
    {
        if((*v) == "X")
        {
            LOG(PROC DISTR PROFILE, Creating profile for X);
            createProfile("X",0);
        }
        else if((*v) == "Y")
        {
            if(_part_set->getDims() == 4 || _part_set->getDims() == 6)
            {
                LOG(PROC DISTR PROFILE, Creating profile for Y);
                createProfile("Y",2);
            }
            else
            {
                WARN(PROC DISTR PROFILE, Invalid profile requested for Y);
            }
        }
        else if((*v) == "Z")
        {
            if(_part_set->getDims() == 6)
            {
                LOG(PROC DISTR PROFILE, Creating profile for Z);
                createProfile("Z",4);
            }
            else
            {
                WARN(PROC DISTR PROFILE, Invalid profile requested for Z);
            }
        }
        else if((*v) == "XP")
        {
            LOG(PROC DISTR PROFILE, Creating profile for XP);
            createProfile("XP",1);
        }
        else if((*v) == "YP")
        {
            if(_part_set->getDims() == 4 || _part_set->getDims() == 6)
            {
                LOG(PROC DISTR PROFILE, Creating profile for YP);
                createProfile("YP",3);
            }
            else
            {
                WARN(PROC DISTR PROFILE, Invalid profile requested for YP);
            }
        }
        else if((*v) == "ZP")
        {
            if(_part_set->getDims() == 6)
            {
                LOG(PROC DISTR PROFILE, Creating profile for ZP);
                createProfile("ZP",5);
            }
            else
            {
                WARN(PROC DISTR PROFILE, Invalid profile requested for ZP);    
            }
        }
    }
    std::list<TH1D*>::iterator i;
    for(i=_profiles_list.begin();i!=_profiles_list.end();i++)
    {
        Log::getInstance().log("PROC DISTR PROFILE", "Adding profile in pad", "id: ", _profile_counter);
        _canvas->cd(_profile_counter+1);
        (*i)->Draw();
        _profile_counter++;
    }
}

void
ProcDistrProfile::createProfile(std::string axis, unsigned short int a)
{ 
    TH1D* pro = new TH1D((axis).c_str(),_title.c_str(),300,static_cast<Double_t>(-box_size_),static_cast<Double_t>(-box_size_));
    _profiles_list.push_back(pro);
    pro->GetXaxis()->SetTitle(axis.c_str());
    pro->GetXaxis()->CenterTitle();
    pro->GetXaxis()->SetTitleSize(static_cast<Float_t>(0.04));   
	for(int i=0;i<_part_set->getParticlesNumber();++i)
	{
	    if((*store)(i,static_cast<int>(a)) != LOST)
	    {
            pro->Fill(static_cast<double>((*store)(i,static_cast<int>(a))));
	    }
	}
    LOG(PROC DISTR PROFILE, Profile created);
}

