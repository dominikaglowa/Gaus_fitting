#include <iostream>
#include <stdio>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib>
#include <TGraph.h>
#include <TH1F.h>
#include <iostream>




void divide(){
    
    TFile*a = new TFile("ratio_test22.04.root","UPDATE");
    
    //sn116_Eg8->Clone("ratio_Eg8");
    //sn124_Eg4->Clone("h4b")
    //h4a->Sumw2();
    //h4b->Sumw2();
    TH1F*ratio=new TH1F("ratioEg5","ratio",300,0,6);
    ratio->Sumw2();
    ratio->SetAxisRange(0,5,"Y");
    ratio->Divide(sn116_Eg5,sn124_Eg5,0.71,1.,"B");
    
    //sn116_Eg9->Clone("ratio_Eg9");
    //ratio_Eg9->Divide(sn116_Eg9,sn124_Eg9,0.71,1.,"B");

    // sn116_Eg5->Clone(“h5a”);
    // sn124_Eg5->Clone(“h5b”);
    // h5a->Sumw2();
    // h5b->Sumw2();
    // TH1F*ratioEg5=h5b;
    // ratioEg5->Divide(h5b,h5a,1.,1.,”B”);
    
    // sn116_Eg6->Clone(“h6a”);
    // sn124_Eg6->Clone(“h6b”);
    // h6a->Sumw2();
    // h6b->Sumw2();
    // TH1F*ratioEg6=h6b;
    // ratioEg6->Divide(h6b,h6a,1.,1.,”B”);
    
    // sn116_Eg7->Clone(“h7a”);
    // sn124_Eg7->Clone(“h7b”);
    // h7a->Sumw2();
    // h7b->Sumw2();
    // TH1F*ratioEg7=h7b;
    // ratioEg7->Divide(h7b,h7a,1.,1.,”B”);
    
    // sn116_Eg8->Clone(“h8a”);
    // sn124_Eg8->Clone(“h8b”);
    // h8a->Sumw2();
    // h8b->Sumw2();
    // TH1F*ratioEg8=h8b;
    // ratioEg8->Divide(h8b,h8a,1.,1.,”B”);
    
    // sn116_Eg9->Clone(“h9a”);
    // sn124_Eg9->Clone(“h9b”);
    // h9a->Sumw2();
    // h9b->Sumw2();
    // TH1F*ratioEg9=h9b;
    // ratioEg9->Divide(h9b,h9a,1.,1.,”B”);
    
    // sn116_Eg10->Clone(“h10a”);
    // sn124_Eg10->Clone(“h10b”);
    // h10a->Sumw2();
    // h10b->Sumw2();
    // TH1F*ratioEg10=h10b;
    // ratioEg10->Divide(h10b,h10a,1.,1.,”B”);
  
    a->Write();
    a->Close();
}


