#include <iostream>
#include <stdio>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib>
#include <TGraph.h>
#include <TH1F.h>
#include <iostream>




void clone124(){
    
    TFile*a = new TFile("DGFitsGausb13tout_sn124.root","UPDATE");

    // hist_areas_Eg4->Clone("sn124_Eg4");

    // hist_areas_Eg5->Clone("sn124_Eg5");

    // hist_areas_Eg6->Clone("sn124_Eg6");

    // hist_areas_Eg7->Clone("sn124_Eg7");

    hist_areas_Eg8->Clone("sn124_Eg8");

    // hist_areas_Eg9->Clone("sn124_Eg9");

    // hist_areas_Eg10->Clone("sn124_Eg10")
      ;
    
    a->Write();
    a->Close();
}

void clone116(){
    
    TFile*a = new TFile("DGFitsGausb13tout_sn116.root","UPDATE");

    // hist_areas_Eg4->Clone("sn116_Eg4");

    // hist_areas_Eg5->Clone("sn116_Eg5");

    // hist_areas_Eg6->Clone("sn116_Eg6");

    // hist_areas_Eg7->Clone("sn116_Eg7");

    hist_areas_Eg8->Clone("sn116_Eg8");

    // hist_areas_Eg9->Clone("sn116_Eg9");

    // hist_areas_Eg10->Clone("sn116_Eg10");
    
    a->Write();
    a->Close();
}

void clone120(){
    
    TFile*a = new TFile("DGFitsGausb13tout_sn120.root","UPDATE");

    hist_areas_Eg4->Clone("sn120_Eg4");

    hist_areas_Eg5->Clone("sn120_Eg5");

    hist_areas_Eg6->Clone("sn120_Eg6");

    hist_areas_Eg7->Clone("sn120_Eg7");

    hist_areas_Eg8->Clone("sn120_Eg8");

    hist_areas_Eg9->Clone("sn120_Eg9");

    hist_areas_Eg10->Clone("sn120_Eg10");
    
    a->Write();
    a->Close();
}

