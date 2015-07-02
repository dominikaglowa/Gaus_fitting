#include <iostream>
#include <stdio>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib>
#include <TGraph.h>
#include <TH1F.h>
#include <iostream>
#include <TRandom2.h>
#include <TMath.h>




using namespace std;



void readinhist(){

   ifstream infile("sn124p_250.out");
   
   
ifstream infile("gampi_sn120_n_fsu30_187.out");

Double_t t1[1800],t2[1800],t3[1800],t4[1800],t5[1800],t6[1800];
char dummy;

for(int j = 0; j< 1800; j++){
  t1[j] = .0;
  t2[j] = .0;
  t3[j] = .0;
  t4[j] = .0;
  t6[j] = 0.0;
 }


if(!infile){
  cout << "Cannot open file!!!" << endl;
 }
 else {

   for(int i=0; i< 1800; i++){
     //while(!infile.eof()){
     infile  >> dummy >> t1[i] >> t2[i] >> t3[i] >> t4[i]; 
     t5[i] = t1[i] - 0.05;
   } 

 }


for(i=0; i<1800; i++) cout << "t1[" << i << "]: " << t1[i] << ", t2[" << i << "]: " << t2[i] <<endl;
 
  
TGraph *gr1 = new TGraph(1800,t5,t2);
  
TRandom2 *fRandom = new TRandom2(0);
  
Double_t sigma_res = 2.0;
Double_t mean_res = 90.05; // just to have everything inside bin center is bin 901
TH1F *testH1 = new TH1F("testH1","testH1",1800,0.0,180.0);
Int_t n_trials = 100000;
Double_t val_random;  
//Creating random distribution with the right sigma 
for (int i=0 ; i<n_trials; i++) {
  val_random = fRandom->Gaus(mean_res,sigma_res);
  testH1->Fill(val_random);
 }
  
Int_t n_counts = int(3*sigma_res/testH1->GetBinWidth(1)); // I will count from -3sigma to 3sigma
  
for (int i=0; i<1800; i++) {
  for (int j=-n_counts; j<n_counts ; j++) {
    if ((i+j)>=0 && (i+j)<1800) {
      t6[i] = t6[i] + t2[i+j] * double(testH1->GetBinContent(901+j)) / double(n_trials);
    }
  }
 }
  
TGraph *gr2 = new TGraph(1800,t5,t6);
  

  
  
  TFile *fout = NULL;
  if( fout )delete fout;
  fout = new TFile("sn124p_250.root","RECREATE");
  gr2->Write();

  TCanvas *c1 = NULL;
  if( c1 )delete c1;
  c1 = new TCanvas("c1");
  gr2->Draw("");
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
