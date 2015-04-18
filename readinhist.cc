#include <iostream>
#include <stdio>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib>
#include <TGraph.h>
#include <TH1F.h>
#include <iostream>





using namespace std;



void readinhist(){

   ifstream infile("sn124p_250.out");

   Double_t t1[1800],t2[1800],t3[1800],t4[1800];
   char dummy;

   for(int j = 0; j< 1800; j++){
      t1[j] = .0;
      t2[j] = .0;
      t3[j] = .0;
      t4[j] = .0;
   }


   if(!infile){
      cout << "Cannot open file!!!" << endl;
   }
   else {

      for(int i=0; i< 1800; i++){
      //while(!infile.eof()){
	infile  >> dummy >> t1[i] >> t2[i] >> t3[i] >> t4[i]; 
      } 

   }


   for(i=0; i<1800; i++) cout << "t1[" << i << "]: " << t1[i] << ", t2[" << i << "]: " << t2[i] <<endl;
 
  
   TH1F *gr = new TH1F("kamalov","kamalov",1800,0.05,180.05);
   for(Int_t j=0; j<1800; j++){
     //gr->SetPoint(j+1,t1[j],t2[j]);
     cout<<t1[j]<<" = "<<gr->GetBinCenter(j+1)<<endl;
      gr->SetBinContent(j+1,t2[j]);
     }

  TFile *fout = NULL;
  if( fout )delete fout;
  fout = new TFile("sn124p_250.root","RECREATE");
  gr->Write();

  TCanvas *c1 = NULL;
  if( c1 )delete c1;
  c1 = new TCanvas("c1");
  gr->Draw("");
}


void scalehist(){

  TFile*a = new TFile("sn116n_200220.root","UPDATE");
  kamalov->Clone("sn116n");
  //sn120p->Scale(0.2324); //sn120: 0.2324, sn118: 0.1728, sn116: 0.1037
  sn116n->Write();
  a->Close();
}


void npcon(){

  TFile*a = new TFile("sn116-124_new.root","UPDATE");
  sn116n->Clone("sn116ncs");
  sn116ncs->Scale(0.689); //sn124: 0.5968, sn120: 0.5833, sn118: 0.5763, sn116: 0.5689, li7:0.5714
  sn116ncs->Write();
  sn116p->Clone("sn116pcs");
  sn116pcs->Scale(0.4310); //sn124: 0.4032, sn120: 0.4156, sn118: 4237, sn116: 0.4310, li7: 0.4286
  sn116pcs->Write();
  a->Close();
}


void npadd(){

  TFile*a = new TFile("sn116_220250.root","UPDATE");
  sn116ncs->Clone("sn116_220250");
  sn116_220250->Add(sn116pcs,1);
  sn116_220250->Write();
  a->Close();
}

void two(){

  TFile*a = new TFile("sn_220250.root","UPDATE");
  sn116->Clone("two");
  two->Add(sn118,1);
  two->Write();
  a->Close();
}

void total(){

  TFile*a = new TFile("sn_200220.root","UPDATE");
  two->Clone("total");
  total->Add(sn120,1);
  total->Write();
  a->Close();
}

void addhist(){
    
    TFile*a = new TFile("sn116-124_new.root","UPDATE");
    Sn124_Eg4->Clone("sn124p2");
    sn124p2->Add(Sn124_Eg5,1);
    sn124p2->Write();
    sn124p2->Clone("sn124p3");
    sn124p3->Add(Sn124_Eg8,1);
    sn124p3->Write();
    sn124p3->Clone("sn124tot");
    sn124tot->Add(Sn124_Eg9,1);
    sn124tot->Write();
    a->Close();
}
