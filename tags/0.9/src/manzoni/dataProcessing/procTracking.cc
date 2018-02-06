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

#include <dataProcessing/procTracking.h>

ProcTracking::ProcTracking(ParticlesDistribution* p, std::string xp) : 
    spacing(13),
    precision(5)
{
    // Initialization
    LOG(PROC TRACKING, Instanciating ProcTracking);
    part_set = p;
    store = part_set->getStore();
    dims = part_set->getDims();
    xpath = xp;
    ics_count = part_set->getTrackingNumber();
    tracking_ic = part_set->getTracking();
    coherent_data = part_set->getCoherent();
    turns = XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/iterator/turns");    
    
    // Check which type of processing is needed
    if(XMLInputParser::getInstance().isFoundFromPath("./flows/henonFlow/dataProcessing/final/tracking/coherent") && part_set->isCoherent())
    {
        withCoherent = true;
        LOG(PROC TRACKING, Coherent activated.);
    }
    else
    {
        withCoherent = false;
    }
    if(XMLInputParser::getInstance().isFoundFromPath("./flows/henonFlow/dataProcessing/final/tracking/particles") && part_set->isTracking())
    {
        withParticles = true;
        LOG(PROC TRACKING, Particles tracking activated.);
    }
    else
    {
        withParticles = false;
    }

    if(withParticles || withCoherent)
    {
        // Then read the parameter 
        readParameters();
        
        // Then process the things
        process();
    }
    else
    {
        WARN(PROC TRACKING, ProcTracking called but nothing is activated !);
    }
}

ProcTracking::~ProcTracking()
{
    LOG(PROC TRACKING, Destructor called);
}

void
ProcTracking::readParameters()
{
    LOG(PROC TRACKING, Reading the input parameters...);
    // Coherent
    if(withCoherent)
    {
        file_name_coherent = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/dataProcessing/final/tracking/coherent","file");
    }
    // Particles
    if(withParticles)
    {
        file_name_particles = XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/dataProcessing/final/tracking/particles","files");
    }
    
    // General
    withFFT = (XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/dataProcessing/final/tracking","fft") == "yes");
    withPlots = (XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/dataProcessing/final/tracking","plots") == "yes");
    if(withFFT)
    {
        withTune = XMLInputParser::getInstance().isFoundFromPath("./flows/henonFlow/dataProcessing/final/tracking/fft/tune");
    }
    else
    {
        withTune = false;
    }
}

void
ProcTracking::process()
{
    LOG(PROC TRACKING, Processing...);
    write2File();
    if(withPlots && (withParticles || withCoherent))
    {
        drawPlots(X_PLANE);
        drawPlots(Y_PLANE);
    }
    if(withFFT && (withParticles || withCoherent))
    {
        fft_window = XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/dataProcessing/final/tracking/fft/window");
        fft_step = XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/dataProcessing/final/tracking/fft/step");
        fft_out_step = XMLInputParser::getInstance().getFirstTexti("./flows/henonFlow/dataProcessing/final/tracking/fft/stepOut");
        if(turns > fft_window)
        {
            computeFFT(X_PLANE); // Will also plot if needed
            computeFFT(Y_PLANE);
        }
        else
        {
            WARN(PROC TRACKING, FFT not done -- Window larger than turns);
        }
    }
}

void
ProcTracking::drawPlots(int dimension)
{
    LOG(PROC TRACKING, Drawing the graphs...);
    
    if(dimension != X_PLANE && dimension != Y_PLANE)
    {
        throw("DrawPlots called for an invalid plane");
    }
    if(dimension == Y_PLANE && part_set->getDims() == TWO_DIMENSIONS)
    {
        throw("DrawPlots called for dimensions greater than the distribution");
    }
     
    // Create the coordinates table    
    double* coord = new double[turns];
    for(int i=0; i<turns;i++)
    {
        coord[i] = i;
    }
    
    TCanvas* canvas = NULL;
    TGraph* graph1 = NULL;
    TGraph* graph2 = NULL;
    double* data1 = NULL;
    double* data2 = NULL;
    std::string file;
    
    // For the tracked particles
    if(withParticles)
    {
        // Loop for each tracked particles
        for(int p = 0; p < ics_count; p++)
        {
            // Set ROOT options and create the canvas
            canvas = new TCanvas(("Tracking"+i2s(p)).c_str(),("Tracking"+i2s(p)).c_str(),300,2*300);
            canvas->SetLeftMargin(static_cast<const Float_t>(0.13));
            canvas->SetRightMargin(static_cast<const Float_t>(0.13));
            canvas->SetBottomMargin(static_cast<const Float_t>(0.13));
            canvas->Divide(static_cast<const Int_t>(1),2);
        
            // Convert the data
            data1 = new double[turns];
            data2 = new double[turns];
            for(int t = 0; t < turns; t++)
            {
                data1[t] = tracking_ic[p][t][dimension];
                data2[t] = tracking_ic[p][t][dimension+1];
            }

            // Create the graph
            canvas->cd(1);
            graph1 = new TGraph(turns, coord, data1); 
            graph1->Draw("AP");
            canvas->Update();
            canvas->cd(2);
            graph2 = new TGraph(turns, coord, data2); 
            graph2->Draw("AP");
            canvas->Update();
        
            // Set the correct file name and save to file
            if(dimension == X_PLANE)
            {
                file = file_name_particles + "MotionHorizontal" + i2s(p) + ".pdf"; 
            }
            else if(dimension == Y_PLANE)
            {
                file = file_name_particles + "MotionVertical" + i2s(p) + ".pdf";
            }
            canvas->SaveAs(file.c_str(), "pdf");
        
            // Clean up
            delete graph1;
            delete graph2;
            delete canvas;
            delete [] data1;
            delete [] data2;
        }
    }
    
    // For the coherent motion of the beam
    if(withCoherent)
    {
        // Set ROOT options and create the canvas
        canvas = new TCanvas("Tracking","Tracking",300,2*300);
        canvas->SetLeftMargin(static_cast<Float_t>(0.13));
        canvas->SetRightMargin(static_cast<Float_t>(0.13));
        canvas->SetBottomMargin(static_cast<Float_t>(0.13));
        canvas->Divide(1,2);
    
        // Convert the data
        data1 = new double[turns];
        data2 = new double[turns];
        for(int t = 0; t < turns; t++)
        {
            data1[t] = coherent_data[t][dimension];
            data2[t] = coherent_data[t][dimension+1];
        }
        
        // Create the graph
        canvas->cd(1);
        graph1 = new TGraph(turns, coord, data1); 
        graph1->Draw("AP");
        canvas->Update();
        canvas->cd(2);
        graph2 = new TGraph(turns, coord, data2); 
        graph2->Draw("AP");
        canvas->Update();
        
        // Set the correct file name and save to file
            if(dimension == X_PLANE)
            {
                file = file_name_coherent + "MotionHorizontal" + ".pdf"; 
            }
            else if(dimension == Y_PLANE)
            {
                file = file_name_coherent + "MotionVertical" + ".pdf";
            }
            canvas->SaveAs(file.c_str(), "pdf");
        
        // Clean up
        delete graph1;
        delete graph2;
        delete canvas;
        delete [] data1;
        delete [] data2;
    }
    
    delete [] coord;
    LOG(PROC TRACKING, Graphs drawed.);
}

void
ProcTracking::computeFFT(int dimension)
{
    LOG(PROC TRACKING, Computing the FFT.);
    
    // Variables
    // Coherent
    double* data_x_coherent = NULL;
    double* data_xp_coherent = NULL;
    double* data_x_coherent_unpacked = NULL;
    double* data_xp_coherent_unpacked = NULL;
    double* data_x_coherent_tmp = NULL;
    double* data_xp_coherent_tmp = NULL;
    std::ofstream file_tune_coherent;
    // Tracking
    double** data_x_tracking = NULL;
    double** data_xp_tracking = NULL;
    double* data_x_tracking_unpacked = NULL;
    double* data_xp_tracking_unpacked = NULL;
    double* data_x_tracking_tmp = NULL;
    double* data_xp_tracking_tmp = NULL;
    std::ofstream file_tune_tracking;
    
    // Tmp string for the file names
    std::string file_name;
    
    //
    // For the coherent motion
    //
    if(withCoherent)
    {
        data_x_coherent = new double[turns];
        data_xp_coherent = new double[turns];
        data_x_coherent_tmp = new double[fft_window];
        data_xp_coherent_tmp = new double[fft_window];
        data_x_coherent_unpacked = new double[2*fft_window]; // To store both real and imaginary parts
        data_xp_coherent_unpacked = new double[2*fft_window]; // To store both real and imaginary parts
        for(int t = 0; t < turns; t++)
        {
            data_x_coherent[t] = coherent_data[t][dimension];
            data_xp_coherent[t] = coherent_data[t][dimension+1];
        }
        if(dimension == X_PLANE)
        {
            file_tune_coherent.open((XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/dataProcessing/final/tracking/fft/tune","file")+"CoherentHorizontal.dat").c_str());
        }
        else if(dimension == Y_PLANE)
        {
            file_tune_coherent.open((XMLInputParser::getInstance().getAttributeTextFromPath("./flows/henonFlow/dataProcessing/final/tracking/fft/tune","file")+"CoherentVertical.dat").c_str());
        }
    }
    
    //
    // For the tracking
    //
    data_x_tracking = new double*[ics_count];
    data_xp_tracking = new double*[ics_count];
    data_x_tracking_tmp = new double[fft_window];
    data_xp_tracking_tmp = new double[fft_window];
    data_x_tracking_unpacked = new double[2*fft_window]; // To store both real and imaginary parts
    data_xp_tracking_unpacked = new double[2*fft_window]; // To store both real and imaginary parts
    for(int tracked = 0; tracked < ics_count; tracked++)
    {
        data_x_tracking[tracked] = new double[turns];
        data_xp_tracking[tracked] = new double[turns];
        for(int t = 0; t < turns; t++)
        {
            data_x_tracking[tracked][t]  = tracking_ic[tracked][t][dimension];
            data_xp_tracking[tracked][t] = tracking_ic[tracked][t][dimension+1];
        }
    }
    
    // Allocate the tune data
    if(withCoherent && withTune)
    {
        coherent_tune_data = new double[(turns-fft_window)/fft_step];
    }
    if(withParticles && withTune)
    {
        tracking_tune_data = new double*[ics_count];
        for(int i = 0; i < ics_count; i++)
        {
            tracking_tune_data[i] = new double[(turns-fft_window)/fft_step];
        }
    }
    
    // For each window
    for(int w = 0; w < static_cast<int>((turns-fft_window)/fft_step); w++)
    {
        #ifdef DEBUG_FLAG
            char* string = new char[BUFFER_SIZE];
            sprintf(string, "FFT for turns in the range %i -- %i", w*fft_step, w*fft_step + fft_window - 1);
            LOGV(PROC TRACKING, string);
            delete string;
        #endif
        
        //
        // Coherent
        //
        if(withCoherent)
        {
            LOG(PROC TRACKING, FFT for coherent.);
            // Copy the windowed data
            for(int i = 0; i < fft_window; i++)
            {
                data_x_coherent_tmp[i] = data_x_coherent[w*fft_step+i];
                data_xp_coherent_tmp[i] = data_xp_coherent[w*fft_step+i];
            }

            // Compute the FFT
            gsl_fft_real_radix2_transform(data_x_coherent_tmp,1,fft_window);
            gsl_fft_real_radix2_transform(data_xp_coherent_tmp,1,fft_window);
        
            // Unpack the result
            gsl_fft_halfcomplex_radix2_unpack(data_x_coherent_tmp, data_x_coherent_unpacked, 1,fft_window);
            gsl_fft_halfcomplex_radix2_unpack(data_xp_coherent_tmp, data_xp_coherent_unpacked, 1,fft_window);
        
            // Set the correct file name and save to file
            if(dimension==0)
            {
                file_name = file_name_coherent + "FFTHorizontal" + i2s(w) + ".dat";
            }
            else if(dimension == 2)
            {
                file_name = file_name_coherent + "FFTVertical" + i2s(w) + ".dat";
            }
            
            if((w % fft_out_step) == 0)
            {
                std::ofstream f;
                f.open(file_name.c_str());
                LOG(PROC TRACKING, Write the FFT in .dat file);
                for(int r = 0 ; r < fft_window; r++)
                {
                    f << std::setw(spacing) << std::setprecision(precision) << std::scientific << data_x_coherent_unpacked[2*r]
                      << std::setw(spacing) << std::setprecision(precision) << std::scientific << data_xp_coherent_unpacked[2*r+1] 
                      << std::endl;
                }
            }
        
            // Obtain the real and the imaginary parts
            double* data_re = new double[fft_window];
            double* data_im = new double[fft_window];
            for(int t = 0; t < fft_window; t++)
            {
                data_re[t] = data_x_coherent_unpacked[2*t];
                data_im[t] = data_xp_coherent_unpacked[2*t+1];
            }

            // Compute the amplitude and the phase
            double* amplitude = new double[fft_window];
            double* phase = new double[fft_window];
            for(int i = 0; i < fft_window; i++)
            {
                amplitude[i] = sqrt(data_re[i]*data_re[i]+data_im[i]*data_im[i]);
                if(data_re[i] != 0)
                {
                    phase[i] = atan(data_im[i]/data_re[i]);
                }
                else
                {
                    phase[i] = 0.0;
                }
            }
            
            delete [] data_re;
            delete [] data_im;
        
            // Find the maximum of the amplitude
            double max_value = 0.0;
            int max_index = 0;
            for(int v=1; v < fft_window; v++) // Start at 1 to avoid the constant term
            {
                if(amplitude[v] > max_value)
                {
                    max_value = amplitude[v];
                    //std::cout << amplitude[v] << std::endl;
                    max_index = v;
                }
            }
           
            LOG(PROC TRACKING, Write the tune in .dat file);
            double tune_coherent = static_cast<double>(max_index)/fft_window;
            file_tune_coherent  <<  tune_coherent << std::endl;
            coherent_tune_data[w] = tune_coherent;
        
            if(withPlots && (w % fft_out_step) == 0)
            {
                // Set the correct file name and save to file
                std::string file = file_name_coherent + "FFTAmplitude" + i2s(w) + ".pdf";
                plotFFT(fft_window, amplitude, file);
                file = file_name_coherent + "FFTPhase" + i2s(w) + ".pdf";
                plotFFT(fft_window, phase, file);
            }
            delete [] amplitude;
            delete [] phase;
        } // If Coherent  
        
        //
        // Tracking
        //
        if(withParticles)
        {
            LOG(PROC TRACKING, FFT for tracking.);
            for(int t = 0; t < ics_count; t++)
            {
                // Copy the windowed data
                for(int i = 0; i < fft_window; i++)
                {
                    data_x_tracking_tmp[i] = data_x_tracking[t][w*fft_step+i];
                    data_xp_tracking_tmp[i] = data_xp_tracking[t][w*fft_step+i];
                }

                // Compute the FFT
                gsl_fft_real_radix2_transform(data_x_tracking_tmp,1,fft_window);
                gsl_fft_real_radix2_transform(data_xp_tracking_tmp,1,fft_window);
        
                // Unpack the result
                gsl_fft_halfcomplex_radix2_unpack(data_x_tracking_tmp, data_x_tracking_unpacked, 1,fft_window);
                gsl_fft_halfcomplex_radix2_unpack(data_xp_tracking_tmp, data_xp_tracking_unpacked, 1,fft_window);
        
                // Set the correct file name and save to file
                if(dimension==X_PLANE)
                {
                    file_name = file_name_particles + i2s(t) + "FFTHorizontal" + i2s(w) + ".dat";
                }
                else if(dimension == Y_PLANE)
                {
                    file_name = file_name_particles + i2s(t) + "FFTVertical" + i2s(w) + ".dat";
                }
            
                if((w % fft_out_step) == 0)
                {
                    std::ofstream f;
                    f.open(file_name.c_str());
                    LOG(PROC TRACKING, Write the FFT in .dat file);
                    for(int r = 0 ; r < fft_window; r++)
                    {
                        f << std::setw(spacing) << std::setprecision(precision) << std::scientific << data_x_tracking_unpacked[2*r]
                        << std::setw(spacing) << std::setprecision(precision) << std::scientific << data_xp_tracking_unpacked[2*r+1] 
                        << std::endl;
                    }
                }
        
                // Obtain the real and the imaginary parts
                double* data_re = new double[fft_window];
                double* data_im = new double[fft_window];
                for(int x = 0; x < fft_window; x++)
                {
                    data_re[t] = data_x_tracking_unpacked[2*x];
                    data_im[t] = data_xp_tracking_unpacked[2*x+1];
                }

                // Compute the amplitude and the phase
                double* amplitude = new double[fft_window];
                double* phase = new double[fft_window];
                for(int i = 0; i < fft_window; i++)
                {
                    amplitude[i] = sqrt(data_re[i]*data_re[i]+data_im[i]*data_im[i]);
                    if(data_re[i] != 0)
                    {
                        phase[i] = atan(data_im[i]/data_re[i]);
                    }
                    else
                    {
                        phase[i] = 0.0;
                    }
                }
            
                delete [] data_re;
                delete [] data_im;
        
                // Find the maximum of the amplitude
                double max_value = 0.0;
                int max_index = 0;
                for(int v=1; v < fft_window; v++) // Start at 1 to avoid the constant term
                {
                    if(amplitude[v] > max_value)
                    {
                        max_value = amplitude[v];
                        //std::cout << amplitude[v] << std::endl;
                        max_index = v;
                    }
                }
           
                LOG(PROC TRACKING, Write the tune in .dat file);
                double tune_tracking = 0;
                double k = static_cast<double>(max_index);
                tune_tracking = (k/fft_window) + (1/(2*M_PI))*asin(tuneFunctionA((1/fft_window)*fabs(amplitude[max_index]),(1/fft_window)*fabs(amplitude[max_index+1]),cos(2*M_PI/fft_window))*sin(2*M_PI/fft_window));
                
                file_tune_tracking  <<  tune_tracking << std::endl;
                tracking_tune_data[t][w] = tune_tracking;
        
                if(withPlots && (w % fft_out_step) == 0)
                {
                    // Set the correct file name and save to file
                    std::string file = file_name_particles + i2s(t) + "FFTAmplitude" + i2s(w) + ".pdf";
                    plotFFT(fft_window, amplitude, file);
                    file = file_name_coherent + i2s(t) + "FFTPhase" + i2s(w) + ".pdf";
                    plotFFT(fft_window, phase, file);
                }
                delete [] amplitude;
                delete [] phase;
            } // For particles
        } // If particles
        
    } // For windows
    LOG(PROC TRACKING, All windows computed.);
    
    // Plot the tunes
    if(withCoherent && withTune)
    {
        std::string file;
        if(dimension == X_PLANE)
        {
            file = file_name_coherent + "TuneHorizontal.pdf";
        }
        else if(dimension == Y_PLANE)
        {
            file = file_name_coherent + "TuneVertical.pdf";
        }
        plotComputedTune(coherent_tune_data,file);
    }
    if(withParticles && withTune)
    {   
        std::string file = "";
        for(int i = 0; i < ics_count; i++)
        {
            if(dimension == X_PLANE)
            {
                file = file_name_particles + i2s(i) + "TuneHorizontal.pdf"; 
            }
            else if(dimension == Y_PLANE)
            {
                file = file_name_particles + i2s(i) + "TuneVertical.pdf"; 
            }
            plotComputedTune(tracking_tune_data[i],file);
        }
    }
    
    // Clean
    if(data_x_coherent != NULL)
    {
        delete [] data_x_coherent;
    }
    if(data_xp_coherent != NULL)
    {
        delete [] data_xp_coherent;
    }
    if(data_x_coherent_unpacked != NULL)
    {
        delete [] data_x_coherent_unpacked;
    }
    if(data_xp_coherent_unpacked != NULL)
    {
        delete [] data_xp_coherent_unpacked;
    }
    if(data_x_coherent_tmp != NULL)
    {
        delete [] data_x_coherent_tmp;
    }
    if(data_xp_coherent_tmp != NULL)
    {
        delete [] data_xp_coherent_tmp;
    }
    file_tune_coherent.close();
    LOG(PROC TRACKING, Compute FFT done.);
    
}

void
ProcTracking::plotComputedTune(double* data, std::string file_name)
{
    LOG(PROC TRACKING, Plotting computed tune.);
    // Create the canvas
    TCanvas* canvas = NULL;
    TGraph* graph = NULL;
    canvas = new TCanvas("TUNE","TUNE",500,500);
    canvas->SetLeftMargin(static_cast<Float_t>(0.13));
    canvas->SetRightMargin(static_cast<Float_t>(0.13));
    canvas->SetBottomMargin(static_cast<Float_t>(0.13));
    
    // Create the coordinates and compute the length
    int length = (turns-fft_window)/fft_step;
    double* coord = new double[length];
    for(int i=0; i<length;i++)
    {
        coord[i] = fft_window + i*fft_step;
    }
    
    // Create the graph
    graph = new TGraph(length, coord, data); 
    graph->Draw("APL");
    
    // Save the graph
    canvas->SaveAs(file_name.c_str(), "pdf");
    
    // Clean
    delete graph;
    delete canvas;
}

void
ProcTracking::plotFFT(int length, double* data, std::string file)
{
    LOG(PROC TRACKING, Plotting the FFT coefficients);
    // Create the coordinates table    
    double* coord = new double[length];
    for(int i=0; i<length;i++)
    {
        coord[i] = i;
    }
    
    // Varibles
    TCanvas* canvas = NULL;
    TGraph* graph1 = NULL;
    
    // For the tracked particles
    if(withCoherent) // Just to be ... coherent
    {
        // Set ROOT options and create the canvas
        canvas = new TCanvas("TrackingFFT","TrackingFFT",500,500);
        canvas->SetLeftMargin(static_cast<Float_t>(0.13));
        canvas->SetRightMargin(static_cast<Float_t>(0.13));
        canvas->SetBottomMargin(static_cast<Float_t>(0.13));

        // Create the graph
        graph1 = new TGraph(length, coord, data); 
        graph1->Draw("APL");
        
        // Save the canvas 
        canvas->SaveAs(file.c_str(), "pdf");
        
        // Clean up
        delete graph1;
        delete canvas;
    }
}

void
ProcTracking::write2File()
{
    LOG(PROC TRACKING, Outputing to file);
    std::string file;
    std::ofstream f;
    
    // Tracking of particles
    if(withParticles)
    {
        for(int i = 0; i < ics_count; i++)
        {
            file = file_name_particles + i2s(i) + ".dat";
            f.open(file.c_str());
            for(int t = 0; t < turns; t+=1)
            {
                f << std::setw(spacing) << std::setprecision(precision) << std::scientific << tracking_ic[i][t][X_INDEX]
                  << std::setw(spacing) << std::setprecision(precision) << std::scientific << tracking_ic[i][t][XP_INDEX];
                if(dims == 4)
                {
                    f << std::setw(spacing) << std::setprecision(precision) << std::scientific << tracking_ic[i][t][Y_INDEX]
                      << std::setw(spacing) << std::setprecision(precision) << std::scientific << tracking_ic[i][t][YP_INDEX];
                }
                f  << std::endl;
            }
            f.close();
        }
    }
    
    // Coherent motion
    if(withCoherent)
    {   
        file = file_name_coherent + ".dat";   
        f.open(file.c_str());
        for(int t=0; t<turns; t+=1)
        {        
            f << std::setw(spacing) << std::setprecision(precision) << std::scientific << coherent_data[t][X_INDEX]
              << std::setw(spacing) << std::setprecision(precision) << std::scientific << coherent_data[t][XP_INDEX];
            if(dims == 4)
            {
                f << std::setw(spacing) << std::setprecision(precision) << std::scientific << coherent_data[t][Y_INDEX]
                  << std::setw(spacing) << std::setprecision(precision) << std::scientific << coherent_data[t][YP_INDEX];
            }
            f << std::endl;
        }
        f.close();
    }
}

double
ProcTracking::tuneFunctionA(double a, double b, double c)
{
    return (-(a+b*c)*(a-b)+b*sqrt(c*c*(a+b)*(a+b)-2*a*b*(2*c*c-c-1)))/(a*a+b*b+2*a*b*c);
}

