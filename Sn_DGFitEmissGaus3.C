Double_t Gmean=1;
Double_t Gwidth=1;
Double_t QQ=0;
Int_t PlotCount =1;


Double_t Gaus2Gaus2(Double_t *x, Double_t *par) { // 2 gaussian fit 
 //standard resolution function
  // if(par[0]<0||par[3]<0) {return -1E9; }
  //  if(par[2]<0||par[5]<0) {return -1E9;}
  // if(par[5]<1.2||par[5]>4)  {return -1E9;}
  //standard fit
   return TMath::Gaus(x[0],Gmean,Gwidth)*par[0]/(2.5066*Gwidth)+TMath::Gaus(x[0],par[1],par[2])*(par[3])/(2.5066*par[2]);
}

Double_t Gaus3(Double_t *x, Double_t *par) { // single gaussian fit for use near coherent maxima
  return TMath::Gaus(x[0],par[1],par[2])*par[0]/(2.5066*par[2]);
}


TH1F* hnewA;
TH1F* hnewB;
TH1F* hclaireA;

void RunFitEmiss(){
 
  TFile *_file0 = TFile::Open("Physics_CB_888_renamed.root"); // th input file
 
  TCanvas *c1 = new TCanvas("c1","c1",800,800);  // setup canvases for plotting
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  TCanvas *c4 = new TCanvas("c4","c4",800,800);
  // TCanvas *c5 = new TCanvas("c5","c5",800,800);
  // TCanvas *c6 = new TCanvas("c6","c6",800,800);
  // TCanvas *c7 = new TCanvas("c7","c7",800,800);
  // TCanvas *c8 = new TCanvas("c8","c8",800,800);
  // TCanvas *c9 = new TCanvas("c9","c9",800,800);

  // TCanvas *c10 = new TCanvas("c10","c10",800,800);
  // TCanvas *c11 = new TCanvas("c11","c11",800,800);
  // TCanvas *c12 = new TCanvas("c12","c12",800,800);
  // TCanvas *c13 = new TCanvas("c13","c13",800,800);
  // TCanvas *c14 = new TCanvas("c14","c14",800,800);
  // TCanvas *c15 = new TCanvas("c15","c15",800,800);
  // TCanvas *c16 = new TCanvas("c16","c16",800,800);

 c1->cd();c1->Divide(5,5);
 c2->cd();c2->Divide(5,4);
 c3->cd();c3->Divide(5,4);
 c4->cd();c4->Divide(5,4);
 // c5->cd();c5->Divide(5,4);
 // c6->cd();c6->Divide(5,4);
 // c7->cd();c7->Divide(5,4);
 // c8->cd();c8->Divide(5,4);
 // c9->cd();c9->Divide(5,4);
 // c10->cd();c10->Divide(5,4);
 // c11->cd();c11->Divide(5,4);
 // c12->cd();c12->Divide(5,4);
 // c13->cd();c13->Divide(5,4);
 // c14->cd();c14->Divide(5,4);
 // c15->cd();c15->Divide(5,4);


 c1->cd(1);

 TString biname="hist_areas_Eg8"; // base name for the histograms to fit
                                  /// chnange this for each bin you want to fit

 // setup histograms for dumping the fit results into 

 TH1F *hist_areas;         // the binning here must be as for the goat output
 hclaireA = new TH1F("hist_areas_Eg8","hist_areas_Eg8",300,0,6.0);
 hnewA=(TH1F*)hclaireA->Clone();
 hnewA->Reset();
 hnewA->SetNameTitle("newA","newA");
 hnewB=(TH1F*)hclaireA->Clone();
 hnewB->Reset();
 hnewB->SetNameTitle("newB","newB");
 hnewC=(TH1F*)hclaireA->Clone();
 hnewC->Reset();
 hnewC->SetNameTitle("newC","newC");
 hnewD=(TH1F*)hclaireA->Clone();
 hnewD->Reset();
 hnewD->SetNameTitle("newD","newD");
 hnewE=(TH1F*)hclaireA->Clone();
 hnewE->Reset();
 hnewE->SetNameTitle("newE","newE");
 hnewF=(TH1F*)hclaireA->Clone();
 hnewF->Reset();
 hnewF->SetNameTitle("newF","newF");
  hnewG=(TH1F*)hclaireA->Clone();
 hnewG->Reset();
 hnewG->SetNameTitle("newG","newG");
 hnewch=(TH1F*)hclaireA->Clone();
 hnewch->Reset();
 hnewch->SetNameTitle("newch","newch");
 TH1F* hEmiss=NULL;
 char name1[100];

 float high=10;
 float low=0;
 TH2F* h2dfile_pr,h2dfile_r ;

 h2dFile_pr=(TH2F*)_file0->Get("pimissen_pr_DE_E_8");
 h2dFile_r=(TH2F*)_file0->Get("pimissen_r_DE_E_8");

 double random_factor = 0.055; // CHECK WITH DAN
 h2dFile_pr->Sumw2(1);
 h2dFile_pr->Add(h2dFile_r,-random_factor);

 // We have 300 bins in q, but we want a range of 20MeV in q, so we need to have a range of 4 before we do the projection (75 bins between 0 and 1.5, in order to have the same separation as 300 bins from 0 to 6 in q) . If you want to use the full binning, you will need to modify the name of the file, the GetSignal function in order to add all the smaller bins in order to cover the coherent peak, and the factor of 4 at the bottom here.
 int bin_q_claire = hclaireA->GetNbinsX()/4 -1;

//do fits with free polynomial parameters
  for(Int_t i=bin_q_claire;i>-1;i--){ //loop over all q bins
   hEmiss=NULL;
   low=(float)i/50;    // calculate low and high limits for q rages in histogram
   high=(float)(i+1)/50;  // this is bin size dependent!
   sprintf(name1,"pimissen_q_8_%4.3f_%4.3f",low,high);
   //  cout << name1 << endl;
 
   h2dFile_pr->GetYaxis()->SetRange(4*i,4*i+3);
   hEmiss=(TH1F*)h2dFile_pr->ProjectionX();
   hEmiss->Rebin(2);
   hEmiss->SetName(name1);
  }


 GetSignal();  //  fix Gaussian mean and Gaussian width by fitting in coherent peak
               // all future fits will take this for the coherent gaussian properties 
// his->Draw("");

 c1->cd(2);

 TF1 *fg2p2 = new TF1("fg2p2",Gaus2Gaus2,-60,60,4); //fit function
 
 //setup some fit ranges and limits for the fitting 
 // fg2p2->SetParLimits(0,0.,1E6);
// fg2p2->SetParLimits(3,0.,1E6);
 Int_t Lim1Min=-70; 
 Int_t Lim1Max=70; 
 Int_t Lim2Min=10; 
 Int_t Lim2Max=120;
 fg2p2->SetParLimits(0,0,1E6);
 fg2p2->SetParLimits(1,Lim1Min,Lim1Max);
 fg2p2->SetParLimits(2,Lim2Min,Lim2Max);
 fg2p2->SetParLimits(3,0,100000);
 // fg2p2->SetParLimits(4,Gmean-2,Gmean+2);
// fg2p2->FixParameter(3,0);

 fg2p2->SetParameters(100,10,0,Lim2Min*1.5);
 

//do fits with free polynomial parameters
  for(Int_t i=bin_q_claire;i>-1;i--){ //loop over all q bins
   hEmiss=NULL;
   low=(float)i/50;    // calculate low and high limits for q rages in histogram
   high=(float)(i+1)/50;  // this is bin size dependent!
   sprintf(name1,"pimissen_q_8_%4.3f_%4.3f",low,high);
   cout << name1 << endl;
   QQ=(low+high)/2;  //centre of the q bin e.g bin 0.22_0.23 is 0.225
   // cout<<name1<<endl;
   //cout << "low, high QQ " << low << " " << high << " " << QQ << endl;
   h2dFile_pr->GetYaxis()->SetRange(4*i,4*i+3);
   hEmiss=(TH1F*)h2dFile_pr->ProjectionX();
   hEmiss->SetName(name1);
   hEmiss->SetTitle(name1);
   hEmiss->Rebin(2);
   //   hEmiss=(TH1F*)_file0->Get(name1); // get the emiss spectra from the file               OLD hemis definition (was getting it from file)
   if(hEmiss){                  
     cout << "Fitting histogram .. " << name1;
     if(hEmiss->Integral()<1)continue; // bomb out if < 1 counts in spectra
     hEmiss->GetXaxis()->UnZoom();  
     c1->cd(3);
     //hEmiss->Rebin(4); // rebin the missing energy spectra
     Int_t result= hEmiss->Fit("fg2p2","RM0"); // fit the emiss spectra with the 2 gaus function
     hEmiss->GetListOfFunctions()->Add(fg2p2->Clone());  
     Double_t p0=  fg2p2->GetParameter(0);  
     fg2p2->SetParameter(0,0); // zero parameter after stored in p0
     hEmiss->GetListOfFunctions()->Add(fg2p2->Clone()); // add fit function to histogram
     hnewA->SetBinContent(i+1,p0); // set bin content to param 0 of fit function (coherent gaussian area)
     hnewA->SetBinError(i+1,fg2p2->GetParError(0)); 
     // hnewA->SetBinError(i+1,fg2p2->Integral(Gmean-2*Gwidth,Gmean+2*Gwidth));

     // Do some clumsy code to put plots in different pads
     cout << "PlotCount" << PlotCount;
     if(PlotCount<20){c1->cd();c1->cd(PlotCount);}  // plot results of fits for all q bins 
     if(PlotCount>20 && PlotCount<40){c2->cd();c2->cd(PlotCount-20);} // use PLotcount which increments
     if(PlotCount>40 && PlotCount<60){c3->cd();c3->cd(PlotCount-40);}  // with each bin
     if(PlotCount>60 && PlotCount<80){c4->cd();c4->cd(PlotCount-60);}
     if(PlotCount>80 && PlotCount<100){c5->cd();c5->cd(PlotCount-80);}
     if(PlotCount>100 && PlotCount<120){c6->cd();c6->cd(PlotCount-100);}
     if(PlotCount>120 && PlotCount<140){c7->cd();c7->cd(PlotCount-120);}
     if(PlotCount>140 && PlotCount<160){c8->cd();c8->cd(PlotCount-140);}
     if(PlotCount>160 && PlotCount<180){c9->cd();c9->cd(PlotCount-160);}
     if(PlotCount>180 && PlotCount<200){c9->cd();c9->cd(PlotCount-180);}
     if(PlotCount>200 && PlotCount<220){c10->cd();c10->cd(PlotCount-200);}
     if(PlotCount>220 && PlotCount<240){c11->cd();c11->cd(PlotCount-220);}
     if(PlotCount>240 && PlotCount<260){c12->cd();c12->cd(PlotCount-240);}
     if(PlotCount>260 && PlotCount<280){c13->cd();c13->cd(PlotCount-260);}
     if(PlotCount>280)continue; // add more later ?
     PlotCount++;
     hEmiss->Draw("PE"); // draw the missing energy spectra with the fits
     //hnewA->Draw("same"); 


     if(fg2p2->GetParameter(2)>Lim2Max-1||fg2p2->GetParameter(2)<Lim2Min+1)continue; // if parameters of backg gaussian outside of
     if(fg2p2->GetParameter(1)>Lim1Max-1||fg2p2->GetParameter(1)<Lim1Min+1)continue; // limits dont write to histogram
     
     hnewB->SetBinContent(i+1,fg2p2->GetParameter(1)); //hnew B contains param 2
     hnewB->SetBinError(i+1,fg2p2->GetParError(1));    // from fit
       
     hnewC->SetBinContent(i+1,fg2p2->GetParameter(2));//hnew c contains param 3
     hnewC->SetBinError(i+1,fg2p2->GetParError(2));   // from fit
       // hnewF->SetBinContent(i+1,fg2p2->GetParameter(4));
       //hnewF->SetBinError(i+1,fg2p2->GetParError(4));
     hnewch->SetBinContent(i+1,fg2p2->GetChisquare()/fg2p2->GetNDF()); // hnewch stores the chi square
 
   }
 }
 //   return;
 // Write unconstrained results

 //to save uncomment out this
   TFile* fitfile=new TFile("DGFitsGausa13tout.root","update");
   hnewA->SetNameTitle(biname,biname);
   hnewA->Write("",TObject::kOverwrite);
   fitfile->Close();

 //Fit the q dependence of the backgroun gaussian parameters (mean and width)
 TF1* myp=new TF1("myp","pol2(0)",0.4,3);
 myp->SetParameters(-24.3009,113.588,-158.613,112.821,-33.7803); // set initial parameters for the fits?
 TF1* myp4=new TF1("myp4","pol4(0)",0.3,3);
 myp4->SetParameters(-24.3009,113.588,-158.613,112.821,-33.7803);
 c1->cd(4);
 hnewB->Fit("myp","R");  // 2nd order Polynomial fit of the extracted parameters from the fit par 1 (back gaus mean)
 TF1* fp1=myp->Clone("fp1");
 c1->cd(5);
 hnewC->Fit("myp","R");  // 4th order Polynomial fit of the extracted parameters from the fit par 2 (back gauss width)
 TF1* fp2=myp->Clone("fp2");
 // hnewF->Fit("myp","R");
 // TF1* fp4=myp->Clone("fp4");
 //return;

 //redo fits with constrained polynomial parameters for background gaussian
 //i.e. use smooth function to eliminate erroneous bin-to-bin fluctuations

 for(Int_t i=0;i<bin_q_claire+1;i++){
   hEmiss=NULL;
   low=(float)i/50;
   high=(float)(i+1)/50;
   sprintf(name1,"pimissen_q_8_%4.3f_%4.3f",low,high);
   QQ=(low+high)/2;  // q bin centre
   cout<<name1<<endl;
   h2dFile_pr->GetYaxis()->SetRange(4*i,4*i+3);
   hEmiss=(TH1F*)h2dFile_pr->ProjectionX();
   hEmiss->SetName(name1);
   hEmiss->Rebin(2);
   //  hEmiss=(TH1F*)_file0->Get(name1); // old
   if(hEmiss){
     hEmiss->GetXaxis()->UnZoom();
     fg2p2->FixParameter(1,fp1->Eval(QQ));   //fix background gaus parameters to
                                             //value @ bin centre of fitted function
     //fg2p2->FixParameter(2,fp2->Eval(QQ));
     fg2p2->SetParameter(2,fp2->Eval(QQ));

     //fg2p2->FixParameter(4,fp4->Eval(QQ));
 
     c1->cd(6);
     hEmiss->Fit("fg2p2","QR");  // redo fit 
     hEmiss->GetListOfFunctions()->Add(fg2p2->Clone()); // add fit to the histogram
     Double_t p0=  fg2p2->GetParameter(0);
     fg2p2->SetParameter(0,0); // zero the coherent gaussian after area stored in p0
     hEmiss->GetListOfFunctions()->Add(fg2p2->Clone()); // draw this as well (i.e. only the background gaussian)
     hnewE->SetBinContent(i+1,p0); // hnewE is size of coherent gaussian
     hnewE->SetBinError(i+1,fg2p2->GetParError(0));
     //hnewD->SetBinContent(i+1,fg2p2->Integral(Gmean-2*Gwidth,Gmean+2*Gwidth));
     //hnewD->SetBinError(i+1,fg2p2->GetParError(3));
     if(fg2p2->GetParameter(3)>1){
       hnewD->SetBinContent(i+1,fg2p2->GetParameter(3)); // hnewD is the width? of Background gaussian
       hnewD->SetBinError(i+1,fg2p2->GetParError(3));
     }
 
   }
 }
 //to save uncomment out this
 TFile* fitfile=new TFile("DGFitsGausb13tout.root","update"); //save with different filename
  hnewE->SetNameTitle(biname,biname);
  hnewE->Write("",TObject::kOverwrite);    // write the histogra of size of coh gaussian.
  fitfile->Close()  ;

  hnewD->Fit("myp4","QR"); // fit q distribution of fit parameters with polynomial
 TF1* fp3=myp4->Clone("fp3");
 ///return;

 
 for(Int_t i=0;i<bin_q_claire+1;i++){
    hEmiss=NULL;
   low=(float)i/50;
   high=(float)(i+1)/50;
   sprintf(name1,"pimissen_q_8_%4.3f_%4.3f",low,high); // change this if change egamma bin
   //  cout << name1 << endl;
   QQ=(low+high)/2;
   cout<<name1<<endl;
   h2dFile_pr->GetYaxis()->SetRange(i,i);
   hEmiss=(TH1F*)h2dFile_pr->ProjectionX();
   hEmiss->SetName(name1);
   hEmiss->Rebin(2);
   //  hEmiss=(TH1F*)_file0->Get(name1); // old
   if(hEmiss){
     //hEmiss->Rebin(2);
     //hEmiss->Scale(0.5);
     hEmiss->GetXaxis()->UnZoom();
     fg2p2->FixParameter(1,fp1->Eval(QQ));   //Take background from polynomial rather than individual fits
     //  fg2p2->FixParameter(2,fp2->Eval(QQ));
     fg2p2->SetParLimits(2,fp2->Eval(QQ),fp2->Eval(QQ)*1.5); // [2] from function value to 1.5*value??
     fg2p2->SetParameter(2,fp2->Eval(QQ)); // [2] fixed to polynomial fit (line above redundant?
     // fg2p2->SetParLimits(2,5,20);
     if(fp3->Eval(QQ)>0)fg2p2->FixParameter(3,fp3->Eval(QQ)); //[3] fixed
     else fg2p2->FixParameter(3,0);
     //fg2p2->FixParameter(4,fp4->Eval(QQ));
     // cout << "PlotCount" << PlotCount;
     // if(PlotCount<20){c1->cd();c1->cd(PlotCount);}  // plot results of fits for all q bins 
     // if(PlotCount>20 && PlotCount<40){c2->cd();c2->cd(PlotCount-20);} // use PLotcount which increments
     // if(PlotCount>40 && PlotCount<60){c3->cd();c3->cd(PlotCount-40);}  // with each bin
     // if(PlotCount>60 && PlotCount<80){c4->cd();c4->cd(PlotCount-60);}
     // if(PlotCount>80 && PlotCount<100){c5->cd();c5->cd(PlotCount-80);}
     // if(PlotCount>100 && PlotCount<120){c6->cd();c6->cd(PlotCount-100);}
     // if(PlotCount>120 && PlotCount<140){c7->cd();c7->cd(PlotCount-120);}
     // if(PlotCount>140 && PlotCount<160){c8->cd();c8->cd(PlotCount-140);}
     // if(PlotCount>160 && PlotCount<180){c9->cd();c9->cd(PlotCount-160);}
     // if(PlotCount>160)continue; // add more later ?
     hEmiss->Fit("fg2p2","QR0");    // fit emiss with background params fixed to polynomial
     hEmiss->GetListOfFunctions()->Add(fg2p2->Clone());
     Double_t p0=  fg2p2->GetParameter(0); // get area of coherent gaussian (when bgrnd fixed polynomial)
     fg2p2->SetParameter(0,0); // reset param 0
     hEmiss->GetListOfFunctions()->Add(fg2p2->Clone()); // ?
     hnewG->SetBinContent(i+1,p0);  // hnew G is fit to coherent with background from polynomial fit of gauss
     hnewG->SetBinError(i+1,fg2p2->GetParError(0));
     hnewF->SetBinContent(i+1,fg2p2->GetParameter(2));
     hnewF->SetBinError(i+1,fg2p2->GetParError(2));
     
     //c2->cd();
     //c2->cd(PlotCount);
     //  PlotCount++;
     // hEmiss->Draw("");
     //hnewG->Draw("same"); 
   }
  }
  TFile* fitfile=new TFile("DGFitsGausc13tout.root","update");
  hnewG->SetNameTitle(biname,biname);
  hnewG->Write("",TObject::kOverwrite);   // write cohpi with background from polynomial
  fitfile->Close()  ;

}



void GetSignal(){
  //Get signal by fitting coherent peak
  // TF1 *fg3 = new TF1("fg3",Gaus3,-7,4,3);//low  e
  TF1 *fg3 = new TF1("fg3",Gaus3,-25,20,12);
  //Fix parameters to 0.5 value;
  fg3->SetParameters(1000,-2,5);

 //TH1D* his=((TH1D*)gDirectory->Get("DeltaE_Missmom50_BeamE8"));
  //TH1D* his=((TH1D*)gDirectory->Get("pimissen_q_1_0.450_0.455"));
  TH1D* his=((TH1D*)gDirectory->Get("pimissen_q_8_0.460_0.480"));
  his->SetName("hisb");
  his->SetTitle("hisb");
 //  his->SetNameTitle("hisb","hisb");
   his->Add(((TH1F*)gDirectory->Get("pimissen_q_8_0.480_0.500")));
   his->Add(((TH1F*)gDirectory->Get("pimissen_q_8_0.500_0.520")));
   his->Add(((TH1F*)gDirectory->Get("pimissen_q_8_0.440_0.460")));
   his->Add(((TH1F*)gDirectory->Get("pimissen_q_8_0.520_0.540")));
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.465_0.470"))); 
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.470_0.475")));
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.475_0.480")));
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.480_0.485")));
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.485_0.490")));
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.490_0.495")));
//   his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.500_0.505")));
//  // //
//  his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.510_0.515")));
//  his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.520_0.525")));
//  his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.530_0.535")));
//  his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.540_0.545")));
// his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.550_0.555")));
//  his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.560_0.565")));

  // TH1F* his=((TH1F*)gDirectory->Get("pimissen_q_1_0.590_0.600")); //change this if change egamma bin
  // his->SetNameTitle("his","his");
  // //his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.390_0.400")));
  // //changed to 40Ca peak
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.600_0.610")));
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.610_0.620")));
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.620_0.630")));
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.630_0.640")));
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.640_0.650")));
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.650_0.660")));
  // his->Add(((TH1F*)gDirectory->Get("pimissen_q_1_0.660_0.670")));

  his->Draw("PE");
  his->Fit("fg3","P");
  Gmean=fg3->GetParameter(1);
  Gwidth=fg3->GetParameter(2);

  cout << "fitted in coherent peak and gaussian mean, width = " << Gmean << "  "  << Gwidth << endl;

}
