#include <TMath.h>
#include <TF1.h>

//0.399 is the 1/sqrt(2/pi) factor to normalize the gaussian function!

Double_t gaus_linear(Double_t *x, Double_t *par){
  Double_t functionvalue;
  functionvalue = 0.399*(par[0]/par[2])*exp(-pow(x[0]-par[1],2.0)/(2.*pow(par[2],2.0)))+par[3]+par[4]*x[0];
  return functionvalue; 
}

double double_gaus_linear(double *x, double *par){

  double functionvalue;
  functionvalue = 0.399*(par[0]/par[2])*exp(-pow(x[0]-par[1],2.0)/(2.*pow(par[2],2.0)))+par[3]+par[4]*x[0] + 0.399*(par[5]/par[7])*exp(-pow(x[0]-par[6],2.0)/(2.*pow(par[7],2.0)));
  return functionvalue;

}

double double_gaus_linear_unresolvable(double *x, double *par){

  double functionvalue;
  functionvalue = 0.399*(par[0]/par[2])*exp(-pow(x[0]-par[1],2.0)/(2.*pow(par[2],2.0)))+par[3]+par[4]*x[0] + 0.399*(par[5]*par[0]/par[7])*exp(-pow(x[0]-par[6],2.0)/(2.*pow(par[7],2.0)));
  return functionvalue;

}

double triple_gaus_linear(double *x, double *par){

double functionvalue;
functionvalue = 0.399*(par[0]/par[2])*exp(-pow(x[0]-par[1],2.0)/(2.*pow(par[2],2.0)))+par[3]+par[4]*x[0] + 0.399*(par[5]/par[7])*exp(-pow(x[0]-par[6],2.0)/(2.*pow(par[7],2.0))) +
    0.399*(par[8]/par[10])*exp(-pow(x[0]-par[9],2.0)/(2.*pow(par[10],2.0)));
return functionvalue;
}

double double_gaus(double *x, double *par){

  double functionvalue;
  functionvalue = 0.399*(par[0]/par[2])*exp(-pow(x[0]-par[1],2.0)/(2.*pow(par[2],2.0))) + 0.399*(par[3]/par[5])*exp(-pow(x[0]-par[4],2.0)/(2.*pow(par[7],2.0)));
  return functionvalue;

}

void calibrate_Ge_4438(const string &inData, TString detName, TString calibSource){
  
  vector<double> data;
  ifstream fin;
  fin.open(inData);
  double x,y,liveTime, realTime;
  while(fin>> x >>y){
    if(x==-2) liveTime = y;
    else if(x==-1) realTime=y;
    else data.push_back(y);
  }
  fin.close();
  
  
  TH1F *h1 = new TH1F("h1","Calibration data",8192,-0.5,8191.5);

  
  for(int i=0;i<data.size();i++){
    h1->SetBinContent(i+1,data[i]);
  }
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1112);
  
  
  Int_t width=1500;
  Int_t hight=900;
  TCanvas *c1 = new TCanvas("c1","c1",width, hight);
  c1->SetWindowSize(width + (width - c1->GetWw()), hight + (hight - c1->GetWh()));
  //h1->Draw();
  TString pdfName = "pdf_results/" + detName + "_Calibration_with_"+ calibSource + "_" + inData.substr(inData.find_last_of('/') + 1) + ".pdf";
  TString pdfFileOpen = pdfName + "[" ;
  c1->Print(pdfFileOpen.Data());
  ofstream output;
  TString peaksSummaryFileName = calibSource + detName +  ".peaksSummary" ;
  output.open(peaksSummaryFileName);
  
  TString rootName = "pdf_results/" + detName + "_Calibration_with_"+ calibSource + "_" + inData.substr(inData.find_last_of('/') + 1) + ".root";
  
  TFile *file = new TFile(rootName, "RECREATE");
  
  std::vector<double> channel, channel_e, sigma, sigma_e, energy;
  
  //preliminary relations of channel vs energy from previous knowledge of the detector
  double m_c = 1.107;
  double b_c = 1.137;
  double m_s = 0.0003626; //slope of the energy vs sigma relation (ch/E) unit
  double b_s = 1.2; // intercept of the nergy vs sigma relation (ch)
  
  double x_lo = -0.5/m_c - b_c/m_c;
  double x_up = 8191.5/m_c - b_c/m_c;
  
  TH1F *h2 = new TH1F("h2","Calibration data",8192,x_lo,x_up);
  double peak_c, peak_s,peak_c2,peak_s2,peak_c3,peak_s3,xlow,xhigh;
  
  for(int i=0;i<data.size();i++){
    h2->SetBinContent(i+1,data[i]);
  }
  
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);
    
  if(calibSource == "Th228"){
    
    //this peak can hardly be seen above the background 
    //Pb212 115 kev Peak
    peak_c = 115.183*m_c + b_c;  //channel 128.6
    peak_s = 115.183*m_s + b_s;
    xlow = peak_c - 5.;
    xhigh = peak_c + 5.;
    TF1* fit115 = new TF1("fit115",gaus_linear,xlow,xhigh,5);
    fit115->SetParNames("A","#mu","#sigma","b","m");
    fit115->SetParameters(30,peak_c,peak_s,247,0.099);
    h1->Fit("fit115","E","",xlow,xhigh);
    h1->SetTitle("Pb212 115.183 keV peak; ADC Channel; Counts");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    h1->Write();
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit115->GetParameter(1));
    channel_e.push_back(fit115->GetParError(1));
    sigma.push_back(fit115->GetParameter(2));
    sigma_e.push_back(fit115->GetParError(2));
    energy.push_back(115.183);
    output << 115.183 << "\t" << fit115->GetParameter(1) << "\t" << fit115->GetParError(1) << "\t" << fit115->GetParameter(2) << "\t" << fit115->GetParError(2) << endl;

    //Th228 131 keV peak not much observable
    peak_c = 131.613*m_c + b_c;
    peak_s = 131.613*m_s + b_s;
    xlow = peak_c - 15.;
    xhigh = peak_c + 15.;
    TF1* fit131 = new TF1("fit131",gaus_linear,xlow,xhigh,5);
    fit131->SetParNames("A","#mu","#sigma","b","m");
    fit131->SetParameters(9e2,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit131","E","",xlow,xhigh);
    h1->SetTitle("Th228 131.613 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    h1->Write();
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit131->GetParameter(1));
    channel_e.push_back(fit131->GetParError(1));
    sigma.push_back(fit131->GetParameter(2));
    sigma_e.push_back(fit131->GetParError(2));
    energy.push_back(131.613);
    
    /*
    //complicated fitting bad chi square
    //238 kev peak from Pb212 and 240 KeV peak from Ra224
    peak_c = 238.632*m_c + b_c;
    peak_s = 238.632*m_s + b_s;
    peak_c2 = 240.986*m_c + b_c;
    peak_s2 = 240.986*m_s + b_s;
    xlow = (peak_c + peak_c2)/2. - 40.;
    xhigh = (peak_c + peak_c2)/2. + 30.;
    TF1* fit238_240 = new TF1("fit238_240",double_gaus_linear_unresolvable,xlow,xhigh,8);
    fit238_240->SetParNames("A_{Pb212_238keV}","#mu_{Pb212_238keV}","#sigma_{Pb212_238keV}","b","m","I_{Relative Intensity Ratio}","#mu_{Ra224_240keV}","#sigma_{Ra224_240keV}","BR_{Pb212_238keV}/A_{Ra224_240keV}");
    fit238_240->SetParameters(9e5,683.5,3.33,3478,-4.668,0.1,689.9,3.2);
    //fit238_240->FixParameter(1,peak_c);
    //fit238_240->FixParameter(2,peak_s);
    //fit238_240->FixParameter(6,peak_c2);
    //fit238_240->FixParameter(7,peak_s2);
    //fit238_240->SetParLimits(1,680,684);
    //fit238_240->SetParLimits(6,684,700);
    //fit238_240->SetParLimits(7,2.5,4.5);
    fit238_240->FixParameter(5,0.0410/0.436);
    h1->Fit("fit238_240","E","",xlow,xhigh);
    h1->SetTitle("Pb212 238.632keV , Ra224 240.986 keV peak; Entries; ADC Channel");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit238_240->GetParameter(1));
    channel_e.push_back(fit238_240->GetParError(1));
    channel.push_back(fit238_240->GetParameter(6));
    channel_e.push_back(fit238_240->GetParError(6));
    sigma.push_back(fit238_240->GetParameter(2));
    sigma_e.push_back(fit238_240->GetParError(2));
    sigma.push_back(fit238_240->GetParameter(7));
    sigma_e.push_back(fit238_240->GetParError(7));
    energy.push_back(238.632);
    energy.push_back(240.986);
    output << 238.632 << "\t" << fit238_240->GetParameter(1) << "\t" << fit238_240->GetParError(1) << "\t" << fit238_240->GetParameter(2) << "\t" << fit238_240->GetParError(2) << endl;
    output << 240.986 << "\t" << fit238_240->GetParameter(6) << "\t" << fit238_240->GetParError(6) << "\t" << fit238_240->GetParameter(7) << "\t" << fit238_240->GetParError(7) << endl;
    */

    /*
    //252 kev gamma from Tl208
    peak_c = 252.61*m_c + b_c;
    peak_s = 252.61*m_s + b_s;
    xlow = peak_c - 5.;
    xhigh = peak_c + 5.;
    TF1* fit252 = new TF1("fit252",gaus_linear,xlow,xhigh,5);
    fit252->SetParNames("A","#mu","#sigma","b","m");
    fit252->SetParameters(3e3,peak_c,peak_s,266,0.999);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit252","E","",xlow,xhigh);
    h1->SetTitle("Tl208 252.61 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit252->GetParameter(1));
    channel_e.push_back(fit252->GetParError(1));
    sigma.push_back(fit252->GetParameter(2));
    sigma_e.push_back(fit252->GetParError(2));
    energy.push_back(252.61);
    output << 252.61 << "\t" << fit252->GetParameter(1) << "\t" << fit252->GetParError(1) << "\t" << fit252->GetParameter(2) << "\t" << fit252->GetParError(2) << endl;
    */
    
    //277 kev gamma from Tl208
    peak_c = 277.371*m_c + b_c;
    peak_s = 277.371*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 11.;
    TF1* fit277 = new TF1("fit277",double_gaus_linear,xlow,xhigh,8);
    fit277->SetParNames("A","#mu","#sigma","b","m","A2","#mu2","#sigma2");
    fit277->SetParameters(4790,peak_c,peak_s,1172.,-0.7423,150.,167.,0.69);
    fit277->FixParameter(7,0.69);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit277","E","",xlow,xhigh);
    h1->SetTitle("Tl208 277.371 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit277->GetParameter(1));
    channel_e.push_back(fit277->GetParError(1));
    sigma.push_back(fit277->GetParameter(2));
    sigma_e.push_back(fit277->GetParError(2));
    energy.push_back(277.371);
    output << 277.371 << "\t" << fit277->GetParameter(1) << "\t" << fit277->GetParError(1) << "\t" << fit277->GetParameter(2) << "\t" << fit277->GetParError(2) << endl;
    

    
    //300 kev peak from Pb212
    peak_c = 300.087*m_c + b_c;
    peak_s = 300.087*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit300 = new TF1("fit300",gaus_linear,xlow,xhigh,5);
    fit300->SetParNames("A","#mu","#sigma","b","m");
    fit300->SetParameters(8953.,peak_c,peak_s,200.,-3.6);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit300","E","",xlow,xhigh);
    h1->SetTitle("Pb212 300.087 peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit300->GetParameter(1));
    channel_e.push_back(fit300->GetParError(1));
    sigma.push_back(fit300->GetParameter(2));
    sigma_e.push_back(fit300->GetParError(2));
    energy.push_back(300.087);
    output << 300.087 << "\t" << fit300->GetParameter(1) << "\t" << fit300->GetParError(1) << "\t" << fit300->GetParameter(2) << "\t" << fit300->GetParError(2) << endl;
    

      
     // it will have a background from 511 keV electron capture
    //Tl208 510 keV peak
    peak_c = 510.77*m_c + b_c;
    peak_s = 510.77*m_s + b_s;
    xlow = peak_c - 30.;
    xhigh = peak_c + 20.;
    TF1* fit510 = new TF1("fit510",gaus_linear,xlow,xhigh,5);
    fit510->SetParNames("A","#mu","#sigma","b","m");
    fit510->SetParameters(6e3,peak_c,peak_s,55.,-0.494);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit510","E","",xlow,xhigh);
    h1->SetTitle("Tl208 510.770 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit510->GetParameter(1));
    channel_e.push_back(fit510->GetParError(1));
    sigma.push_back(fit510->GetParameter(2));
    sigma_e.push_back(fit510->GetParError(2));
    energy.push_back(510.770);
    

    //Tl208 583 kev peak
    peak_c = 583.191*m_c + b_c;
    peak_s = 583.191*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit583 = new TF1("fit583",gaus_linear,xlow,xhigh,5);
    fit583->SetParNames("A","#mu","#sigma","b","m");
    fit583->SetParameters(5e4,peak_c,peak_s,30,-0.5);//1669
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit583","E","",xlow,xhigh);
    h1->SetTitle("Tl208 583.191 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit583->GetParameter(1));
    channel_e.push_back(fit583->GetParError(1));
    sigma.push_back(fit583->GetParameter(2));
    sigma_e.push_back(fit583->GetParError(2));
    energy.push_back(583.191);
    output << 583.191 << "\t" << fit583->GetParameter(1) << "\t" << fit583->GetParError(1) << "\t" << fit583->GetParameter(2) << "\t" << fit583->GetParError(2) << endl;
    
    
    //727 kev gamma from Bi212
    peak_c = 727.33*m_c + b_c;
    peak_s = 727.33*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit727 = new TF1("fit727",gaus_linear,xlow,xhigh,5);
    fit727->SetParNames("A","#mu","#sigma","b","m");
    fit727->SetParameters(1e4,peak_c,peak_s,1322,-2.673);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit727","E","",xlow,xhigh);
    h1->SetTitle("Bi212 727.33 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit727->GetParameter(1));
    channel_e.push_back(fit727->GetParError(1));
    sigma.push_back(fit727->GetParameter(2));
    sigma_e.push_back(fit727->GetParError(2));
    energy.push_back(727.33);
    output << 727.33 << "\t" << fit727->GetParameter(1) << "\t" << fit727->GetParError(1) << "\t" << fit727->GetParameter(2) << "\t" << fit727->GetParError(2) << endl;
   
  
    /*
    // too close and smaller than 785 keV peak
    //763 kev gamma from Tl208
    peak_c = 763.13*m_c + b_c;
    peak_s = 763.13*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit763 = new TF1("fit763",gaus_linear,xlow,xhigh,5);
    fit763->SetParNames("A","#mu","#sigma","b","m");
    fit763->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit763","E","",xlow,xhigh);
    h1->SetTitle("Tl208 763.13 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit763->GetParameter(1));
    channel_e.push_back(fit763->GetParError(1));
    sigma.push_back(fit763->GetParameter(2));
    sigma_e.push_back(fit763->GetParError(2));
    energy.push_back(763.13);
    output << 763.13 << "\t" << fit763->GetParameter(1) << "\t" << fit763->GetParError(1) << "\t" << fit763->GetParameter(2) << "\t" << fit763->GetParError(2) << endl;
    */

    /*
    //785 kev gamma from Bi212
    peak_c = 785.37*m_c + b_c;
    peak_s = 785.37*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit785 = new TF1("fit785",gaus_linear,xlow,xhigh,5);
    fit785->SetParNames("A","#mu","#sigma","b","m");
    fit785->SetParameters(1622,peak_c,peak_s,26.5,-0.2344); //(500,2247.5,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit785","E","",xlow,xhigh);
    h1->SetTitle("Bi212 785.37 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit785->GetParameter(1));
    channel_e.push_back(fit785->GetParError(1));
    sigma.push_back(fit785->GetParameter(2));
    sigma_e.push_back(fit785->GetParError(2));
    energy.push_back(785.37);
    output << 785.37 << "\t" << fit785->GetParameter(1) << "\t" << fit785->GetParError(1) << "\t" << fit785->GetParameter(2) << "\t" << fit785->GetParError(2) << endl;
    */


   //Tl208 860 kev peak
    peak_c = 860.564*m_c + b_c;
    peak_s = 860.564*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit860 = new TF1("fit860",gaus_linear,xlow,xhigh,5);
    fit860->SetParNames("A","#mu","#sigma","b","m");
    fit860->SetParameters(6913,peak_c,peak_s,574.,-0.878);//(1e3,2464,peak_s,7780.,-3.1)
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit860","E","",xlow,xhigh);
    h1->SetTitle("Tl208 860.564 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit860->GetParameter(1));
    channel_e.push_back(fit860->GetParError(1));
    sigma.push_back(fit860->GetParameter(2));
    sigma_e.push_back(fit860->GetParError(2));
    energy.push_back(860.564);
    output << 860.564 << "\t" << fit860->GetParameter(1) << "\t" << fit860->GetParError(1) << "\t" << fit860->GetParameter(2) << "\t" << fit860->GetParError(2) << endl;
    

    /*
    //Bi212 893 kev peak
    peak_c = 893.408*m_c + b_c;
    peak_s = 893.408*m_s + b_s;
    xlow = peak_c - 8.;
    xhigh = peak_c + 8.;
    TF1* fit893 = new TF1("fit893",gaus_linear,xlow,xhigh,5);
    fit893->SetParNames("A","#mu","#sigma","b","m");
    fit893->SetParameters(507.,peak_c,peak_s,3.,0.995);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit893","E","",xlow,xhigh);
    h1->SetTitle("Bi212 893.408 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit893->GetParameter(1));
    channel_e.push_back(fit893->GetParError(1));
    sigma.push_back(fit893->GetParameter(2));
    sigma_e.push_back(fit893->GetParError(2));
    energy.push_back(893.408);
    output << 893.408 << "\t" << fit893->GetParameter(1) << "\t" << fit893->GetParError(1) << "\t" << fit893->GetParameter(2) << "\t" << fit893->GetParError(2) << endl;
    */

    /*
    //Bi212 1078 kev peak
    peak_c = 1078.62*m_c + b_c;
    peak_s = 1078.62*m_s + b_s;
    xlow = peak_c - 40.;
    xhigh = peak_c + 30.;
    TF1* fit1078_v1 = new TF1("fit1078",gaus_linear,xlow,xhigh,5);
    fit1078_v1->SetParNames("A","#mu","#sigma","b","m");
    fit1078_v1->SetParameters(5.4e4,peak_c,peak_s,7780.,-3.1);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1078","E","",xlow,xhigh);
    h1->SetTitle("Bi212 1078.62 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1078_v1->GetParameter(1));
    channel_e.push_back(fit1078_v1->GetParError(1));
    sigma.push_back(fit1078_v1->GetParameter(2));
    sigma_e.push_back(fit1078_v1->GetParError(2));
    energy.push_back(1078.62);
    output << 1078.62 << "\t" << fit1078_v1->GetParameter(1) << "\t" << fit1078_v1->GetParError(1) << "\t" << fit1078_v1->GetParameter(2) << "\t" << fit1078_v1->GetParError(2) << endl;
    */
    
    /*
    //Tl208 1093 kev peak
    peak_c = 1093.9*m_c + b_c;
    peak_s = 1093.9*m_s + b_s;
    xlow = peak_c - 9.;
    xhigh = peak_c + 9.;
    TF1* fit1093 = new TF1("fit1093",gaus_linear,xlow,xhigh,5);
    fit1093->SetParNames("A","#mu","#sigma","b","m");
    fit1093->SetParameters(1046,peak_c,peak_s,21.7,0.1076);
    fit1093->SetParLimits(2,1.5,5.);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1093","E","",xlow,xhigh);
    h1->SetTitle("Tl208 1093.9 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1093->GetParameter(1));
    channel_e.push_back(fit1093->GetParError(1));
    sigma.push_back(fit1093->GetParameter(2));
    sigma_e.push_back(fit1093->GetParError(2));
    energy.push_back(1093.9);
    output << 1093.9 << "\t" << fit1093->GetParameter(1) << "\t" << fit1093->GetParError(1) << "\t" << fit1093->GetParameter(2) << "\t" << fit1093->GetParError(2) << endl;
    */


    /*
   //Bi212 1078 kev and Tl208 1093.9 kev peak together
    peak_c = 1078.62*m_c + b_c;
    peak_s = 1078.62*m_s + b_s;
    peak_c2 = 1093.9*m_c + b_c;
    peak_s2 = 1093.9*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c2 + 10.;
    TF1* fit1078 = new TF1("fit1078",double_gaus_linear,xlow,xhigh,8);
    fit1078->SetParNames("A","#mu","#sigma","b","m","A2","#mu2","#sigma2");
    fit1078->SetParameters(650.,peak_c,peak_s,500.,-0.87,1046,peak_c2,peak_s2);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1078","E","",xlow,xhigh);
    h1->SetTitle("Bi212 1078.62 kev and Tl208 1093.9 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
   //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1078->GetParameter(1));
    channel_e.push_back(fit1078->GetParError(1));
    sigma.push_back(fit1078->GetParameter(2));
    sigma_e.push_back(fit1078->GetParError(2));
    energy.push_back(1078.62);
    output << 1078.62 << "\t" << fit1078->GetParameter(1) << "\t" << fit1078->GetParError(1) << "\t" << fit1078->GetParameter(2) << "\t" << fit1078->GetParError(2) << endl;
    channel.push_back(fit1078->GetParameter(6));
    channel_e.push_back(fit1078->GetParError(6));
    sigma.push_back(fit1078->GetParameter(7));
    sigma_e.push_back(fit1078->GetParError(7));
    energy.push_back(1093.9);
    output << 1093.9 << "\t" << fit1078->GetParameter(6) << "\t" << fit1078->GetParError(6) << "\t" << fit1078->GetParameter(7) << "\t" << fit1078->GetParError(7) << endl;
    */
    
    
    //Bi212 1512 kev peak
    peak_c = 1512.7*m_c + b_c;
    peak_s = 1512.7*m_s + b_s;
    xlow = peak_c - 40.;
    xhigh = peak_c + 20.;
    TF1* fit1512 = new TF1("fit1512",gaus_linear,xlow,xhigh,5);
    fit1512->SetParNames("A","#mu","#sigma","b","m");
    fit1512->SetParameters(600.,peak_c,peak_s,100.,-0.1);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1512","E","",xlow,xhigh);
    h1->SetTitle("Bi212 1512.7 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1512->GetParameter(1));
    channel_e.push_back(fit1512->GetParError(1));
    sigma.push_back(fit1512->GetParameter(2));
    sigma_e.push_back(fit1512->GetParError(2));
    energy.push_back(1512.7);
    output << 1512.7 << "\t" << fit1512->GetParameter(1) << "\t" << fit1512->GetParError(1) << "\t" << fit1512->GetParameter(2) << "\t" << fit1512->GetParError(2) << endl;
    

    /*
    //Tl208 double escape 1592 kev peak : the pair production both of them escapes
    peak_c = 1592.533*m_c + b_c; //0.511 or 0.533?
    peak_s = 1592.533*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit1592 = new TF1("fit1592",gaus_linear,xlow,xhigh,5);
    fit1592->SetParNames("A","#mu","#sigma","b","m");
    fit1592->SetParameters(823.,peak_c,peak_s,-139.,0.2356);//
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1592","E","",xlow,xhigh);
    h1->SetTitle("Tl208 1592.511 kev DE Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1592->GetParameter(1));
    channel_e.push_back(fit1592->GetParError(1));
    sigma.push_back(fit1592->GetParameter(2));
    sigma_e.push_back(fit1592->GetParError(2));
    energy.push_back(1592.533);
    output << 1592.533 << "\t" << fit1592->GetParameter(1) << "\t" << fit1592->GetParError(1) << "\t" << fit1592->GetParameter(2) << "\t" << fit1592->GetParError(2) << endl;
    */
    /*
    //Bi212 1620 kev peak
    peak_c = 1620.5*m_c + b_c;
    peak_s = 1620.5*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit1620 = new TF1("fit1620",gaus_linear,xlow,xhigh,5);
    fit1620->SetParNames("A","#mu","#sigma","b","m");
    fit1620->SetParameters(1794,peak_c,peak_s,-217.8,0.3186);//(1000,4634,peak_s,25.79,-0.002964)
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1620","E","",xlow,xhigh);
    h1->SetTitle("Bi212 1620.5 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1620->GetParameter(1));
    channel_e.push_back(fit1620->GetParError(1));
    sigma.push_back(fit1620->GetParameter(2));
    sigma_e.push_back(fit1620->GetParError(2));
    energy.push_back(1620.5);
    output << 1620.5 << "\t" << fit1620->GetParameter(1) << "\t" << fit1620->GetParError(1) << "\t" << fit1620->GetParameter(2) << "\t" << fit1620->GetParError(2) << endl;
    */

    /*
    //Tl208 single escape 2103 kev peak : the pair production one of them escapes
    peak_c = 2103.533*m_c + b_c;
    peak_s = 2103.533*m_s + b_s;
    xlow = peak_c - 10.;
    xhigh = peak_c + 10.;
    TF1* fit2103 = new TF1("fit2103",gaus_linear,xlow,xhigh,5);
    fit2103->SetParNames("A","#mu","#sigma","b","m");
    fit2103->SetParameters(1000,peak_c,peak_s,500.,-1.445);//(500,6010,peak_s,-3.,0.5)
    //fit2103->SetParLimits(1,6005,6015);
    //fit2103->SetParLimits(2,2.5,5.5);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit2103","E","",xlow,xhigh);
    h1->SetTitle("Tl208 2103.511 kev SE Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit2103->GetParameter(1));
    channel_e.push_back(fit2103->GetParError(1));
    sigma.push_back(fit2103->GetParameter(2));
    sigma_e.push_back(fit2103->GetParError(2));
    energy.push_back(2103.533);
    output << 2103.533 << "\t" << fit2103->GetParameter(1) << "\t" << fit2103->GetParError(1) << "\t" << fit2103->GetParameter(2) << "\t" << fit2103->GetParError(2) << endl;
    */
    
    //Tl208 2614 keV peak
    peak_c = 2614.533*m_c + b_c;
    peak_s = 2614.533*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 30.;

    TF1* fit2614 = new TF1("fit2614",gaus_linear,xlow,xhigh,5);
    fit2614->SetParNames("A","#mu","#sigma","b","m");
    fit2614->SetParameters(4000.,2907,2,0.005,0.005);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit2614","E","",xlow,xhigh);
    h1->SetTitle("Tl208 2614.533 kev Peak");
    h1->GetXaxis()->SetRangeUser(xlow,xhigh);
    //h1->Fit("gaus");
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit2614->GetParameter(1));
    channel_e.push_back(fit2614->GetParError(1));
    sigma.push_back(fit2614->GetParameter(2));
    sigma_e.push_back(fit2614->GetParError(2));
    energy.push_back(2614.533);
    
    
  }
  else if (calibSource == "Ba133"){
     //Ba133 53 keV peak
    peak_c = 53.1622*m_c + b_c;
    peak_s = 53.1622*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit53 = new TF1("fit53",gaus_linear,xlow,xhigh,5);
    fit53->SetParNames("A","#mu","#sigma","b","m");
    fit53->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit53","E","",xlow,xhigh);
    h1->SetTitle("Ba133 53.1622 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit53->GetParameter(1));
    channel_e.push_back(fit53->GetParError(1));
    sigma.push_back(fit53->GetParameter(2));
    sigma_e.push_back(fit53->GetParError(2));
    energy.push_back(53.1622);
    output << 53.1622 << "\t" << fit53->GetParameter(1) << "\t" << fit53->GetParError(1) << "\t" << fit53->GetParameter(2) << "\t" << fit53->GetParError(2) << endl;
    
    //Ba133 223 keV peak
    peak_c = 223.2368*m_c + b_c;
    peak_s = 223.2368*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit223 = new TF1("fit223",gaus_linear,xlow,xhigh,5);
    fit223->SetParNames("A","#mu","#sigma","b","m");
    fit223->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit223","E","",xlow,xhigh);
    h1->SetTitle("Ba133 223.2368 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit223->GetParameter(1));
    channel_e.push_back(fit223->GetParError(1));
    sigma.push_back(fit223->GetParameter(2));
    sigma_e.push_back(fit223->GetParError(2));
    energy.push_back(223.2368);
    output << 223.2368 << "\t" << fit223->GetParameter(1) << "\t" << fit223->GetParError(1) << "\t" << fit223->GetParameter(2) << "\t" << fit223->GetParError(2) << endl;
    
    //Ba133 276 keV peak
    peak_c = 276.3989*m_c + b_c;
    peak_s = 276.3989*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit276 = new TF1("fit276",gaus_linear,xlow,xhigh,5);
    fit276->SetParNames("A","#mu","#sigma","b","m");
    fit276->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit276","E","",xlow,xhigh);
    h1->SetTitle("Ba133 276.3989 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit276->GetParameter(1));
    channel_e.push_back(fit276->GetParError(1));
    sigma.push_back(fit276->GetParameter(2));
    sigma_e.push_back(fit276->GetParError(2));
    energy.push_back(276.3989);
    output << 276.3989 << "\t" << fit276->GetParameter(1) << "\t" << fit276->GetParError(1) << "\t" << fit276->GetParameter(2) << "\t" << fit276->GetParError(2) << endl;
    
     //Ba133 302 keV peak
    peak_c = 302.8508*m_c + b_c;
    peak_s = 302.8508*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit302 = new TF1("fit302",gaus_linear,xlow,xhigh,5);
    fit302->SetParNames("A","#mu","#sigma","b","m");
    fit302->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit302","E","",xlow,xhigh);
    h1->SetTitle("Ba133 302.8508 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit302->GetParameter(1));
    channel_e.push_back(fit302->GetParError(1));
    sigma.push_back(fit302->GetParameter(2));
    sigma_e.push_back(fit302->GetParError(2));
    energy.push_back(302.8508);
    output << 302.8508 << "\t" << fit302->GetParameter(1) << "\t" << fit302->GetParError(1) << "\t" << fit302->GetParameter(2) << "\t" << fit302->GetParError(2) << endl;
    
    //Ba133 356 keV peak
    peak_c = 356.0129*m_c + b_c;
    peak_s = 356.0129*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit356 = new TF1("fit356",gaus_linear,xlow,xhigh,5);
    fit356->SetParNames("A","#mu","#sigma","b","m");
    fit356->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit356","E","",xlow,xhigh);
    h1->SetTitle("Ba133 356.0129 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit356->GetParameter(1));
    channel_e.push_back(fit356->GetParError(1));
    sigma.push_back(fit356->GetParameter(2));
    sigma_e.push_back(fit356->GetParError(2));
    energy.push_back(356.0129);
    output << 356.0129 << "\t" << fit356->GetParameter(1) << "\t" << fit356->GetParError(1) << "\t" << fit356->GetParameter(2) << "\t" << fit356->GetParError(2) << endl;
    
    //Ba133 383 keV peak
    peak_c = 383.8485*m_c + b_c;
    peak_s = 383.8485*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit383 = new TF1("fit383",gaus_linear,xlow,xhigh,5);
    fit383->SetParNames("A","#mu","#sigma","b","m");
    fit383->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit383","E","",xlow,xhigh);
    h1->SetTitle("Ba133 383.8485 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    h1->Write();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit383->GetParameter(1));
    channel_e.push_back(fit383->GetParError(1));
    sigma.push_back(fit383->GetParameter(2));
    sigma_e.push_back(fit383->GetParError(2));
    energy.push_back(383.8485);
    output << 383.8485 << "\t" << fit383->GetParameter(1) << "\t" << fit383->GetParError(1) << "\t" << fit383->GetParameter(2) << "\t" << fit383->GetParError(2) << endl;
    
  }

  else if(calibSource == "Ra226"){
    /*
     //46 kev peak from Pb2210 and 53 KeV peak from Pb214
    peak_c = 46.539*m_c + b_c;
    peak_s = 46.539*m_s + b_s;
    peak_c2 = 53.2284*m_c + b_c;
    peak_s2 = 53.2284*m_s + b_s;
    xlow = (peak_c + peak_c2)/2. - 30.;
    xhigh = (peak_c + peak_c2)/2. + 30.;
    TF1* fit46_53 = new TF1("fit46_53",double_gaus_linear,xlow,xhigh,8);
    fit46_53->SetParNames("A","#mu","#sigma","b","m","A_2","#mu_2","#sigma_2");
    fit46_53->SetParameters(9e2,peak_c,peak_s,10.,0.1,2e2,peak_c2,peak_s2);
    fit46_53->FixParameter(1,peak_c);
    //fit238_240->FixParameter(2,peak_s);
    fit46_53->FixParameter(6,peak_c2);
    //fit238_240->FixParameter(7,peak_s2);
    h1->GetXaxis()->SetRange(0,0);
    //h1->Fit("fit46_53","E","",xlow,xhigh);
    h1->Fit("fit46_53","E","",110.,162.);
    h1->SetTitle("Pb210 46.539keV , Pb214 53.2284 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit46_53->GetParameter(1));
    channel_e.push_back(fit46_53->GetParError(1));
    channel.push_back(fit46_53->GetParameter(6));
    channel_e.push_back(fit46_53->GetParError(6));
    sigma.push_back(fit46_53->GetParameter(2));
    sigma_e.push_back(fit46_53->GetParError(2));
    sigma.push_back(fit46_53->GetParameter(7));
    sigma_e.push_back(fit46_53->GetParError(7));
    energy.push_back(46.539);
    energy.push_back(53.2284);
    output << 46.539 << "\t" << fit46_53->GetParameter(1) << "\t" << fit46_53->GetParError(1) << "\t" << fit46_53->GetParameter(2) << "\t" << fit46_53->GetParError(2) << endl;
    output << 53.2284 << "\t" << fit46_53->GetParameter(6) << "\t" << fit46_53->GetParError(6) << "\t" << fit46_53->GetParameter(7) << "\t" << fit46_53->GetParError(7) << endl;
    */
    //Rn222 186 keV peak
    peak_c = 186.211*m_c + b_c;
    peak_s = 186.211*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit186 = new TF1("fit186",gaus_linear,xlow,xhigh,5);
    fit186->SetParNames("A","#mu","#sigma","b","m");
    fit186->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit186","E","",xlow,xhigh);
    h1->SetTitle("Rn222 186.211 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit186->GetParameter(1));
    channel_e.push_back(fit186->GetParError(1));
    sigma.push_back(fit186->GetParameter(2));
    sigma_e.push_back(fit186->GetParError(2));
    energy.push_back(186.211);
    output << 186.211 << "\t" << fit186->GetParameter(1) << "\t" << fit186->GetParError(1) << "\t" << fit186->GetParameter(2) << "\t" << fit186->GetParError(2) << endl;
    
    //Pb214 242 keV peak
    peak_c = 241.9950*m_c + b_c;
    peak_s = 241.9950*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit242 = new TF1("fit242",gaus_linear,xlow,xhigh,5);
    fit242->SetParNames("A","#mu","#sigma","b","m");
    fit242->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit242","E","",xlow,xhigh);
    h1->SetTitle("Pb214 241.9950 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit242->GetParameter(1));
    channel_e.push_back(fit242->GetParError(1));
    sigma.push_back(fit242->GetParameter(2));
    sigma_e.push_back(fit242->GetParError(2));
    energy.push_back(241.9950);
    output << 241.995 << "\t" << fit242->GetParameter(1) << "\t" << fit242->GetParError(1) << "\t" << fit242->GetParameter(2) << "\t" << fit242->GetParError(2) << endl;
    
    //Pb214 258 keV peak
    peak_c = 258.86*m_c + b_c;
    peak_s = 258.86*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit258 = new TF1("fit258",gaus_linear,xlow,xhigh,5);
    fit258->SetParNames("A","#mu","#sigma","b","m");
    fit258->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit258","E","",xlow,xhigh);
    h1->SetTitle("Pb214 258.86 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit258->GetParameter(1));
    channel_e.push_back(fit258->GetParError(1));
    sigma.push_back(fit258->GetParameter(2));
    sigma_e.push_back(fit258->GetParError(2));
    energy.push_back(258.86);
    output << 258.86 << "\t" << fit258->GetParameter(1) << "\t" << fit258->GetParError(1) << "\t" << fit258->GetParameter(2) << "\t" << fit258->GetParError(2) << endl;
    
    //Pb214 352 keV peak Bi214 352 keV peak might also be buried under
    peak_c = 351.9321*m_c + b_c;
    peak_s = 351.9321*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit352 = new TF1("fit352",gaus_linear,xlow,xhigh,5);
    fit352->SetParNames("A","#mu","#sigma","b","m");
    fit352->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit352","E","",xlow,xhigh);
    h1->SetTitle("Pb214 351.9321 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit352->GetParameter(1));
    channel_e.push_back(fit352->GetParError(1));
    sigma.push_back(fit352->GetParameter(2));
    sigma_e.push_back(fit352->GetParError(2));
    energy.push_back(351.9321);
    output << 351.9321 << "\t" << fit352->GetParameter(1) << "\t" << fit352->GetParError(1) << "\t" << fit352->GetParameter(2) << "\t" << fit352->GetParError(2) << endl;
    
    //Rn222 510 keV peak
    peak_c = 510*m_c + b_c;
    peak_s = 510*m_s + b_s;
    xlow = peak_c - 30.;
    xhigh = peak_c + 30.;
    TF1* fit510 = new TF1("fit510",gaus_linear,xlow,xhigh,5);
    fit510->SetParNames("A","#mu","#sigma","b","m");
    fit510->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit510","E","",xlow,xhigh);
    h1->SetTitle("Rn222 510 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit510->GetParameter(1));
    channel_e.push_back(fit510->GetParError(1));
    sigma.push_back(fit510->GetParameter(2));
    sigma_e.push_back(fit510->GetParError(2));
    energy.push_back(510.0);
    output << 510.0 << "\t" << fit510->GetParameter(1) << "\t" << fit510->GetParError(1) << "\t" << fit510->GetParameter(2) << "\t" << fit510->GetParError(2) << endl;
    
    //Bi214 609 keV peak
    peak_c = 609.32*m_c + b_c;
    peak_s = 609.32*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit609 = new TF1("fit609",gaus_linear,xlow,xhigh,5);
    fit609->SetParNames("A","#mu","#sigma","b","m");
    fit609->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit609","E","",xlow,xhigh);
    h1->SetTitle("Bi214 609.32 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit609->GetParameter(1));
    channel_e.push_back(fit609->GetParError(1));
    sigma.push_back(fit609->GetParameter(2));
    sigma_e.push_back(fit609->GetParError(2));
    energy.push_back(609.32);
    output << 609.32 << "\t" << fit609->GetParameter(1) << "\t" << fit609->GetParError(1) << "\t" << fit609->GetParameter(2) << "\t" << fit609->GetParError(2) << endl;
    
    //Bi214 665 keV peak
    peak_c = 665.447*m_c + b_c;
    peak_s = 665.447*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit665 = new TF1("fit665",gaus_linear,xlow,xhigh,5);
    fit665->SetParNames("A","#mu","#sigma","b","m");
    fit665->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit665","E","",xlow,xhigh);
    h1->SetTitle("Bi214 665.447 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit665->GetParameter(1));
    channel_e.push_back(fit665->GetParError(1));
    sigma.push_back(fit665->GetParameter(2));
    sigma_e.push_back(fit665->GetParError(2));
    energy.push_back(665.447);
    output << 665.447 << "\t" << fit665->GetParameter(1) << "\t" << fit665->GetParError(1) << "\t" << fit665->GetParameter(2) << "\t" << fit665->GetParError(2) << endl;
    
     //Po210 803 keV peak
    peak_c = 803.06*m_c + b_c;
    peak_s = 803.06*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 40.;
    TF1* fit803 = new TF1("fit803",gaus_linear,xlow,xhigh,5);
    fit803->SetParNames("A","#mu","#sigma","b","m");
    fit803->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit803","E","",xlow,xhigh);
    h1->SetTitle("Bi214 803.06 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit803->GetParameter(1));
    channel_e.push_back(fit803->GetParError(1));
    sigma.push_back(fit803->GetParameter(2));
    sigma_e.push_back(fit803->GetParError(2));
    energy.push_back(803.06);
    output << 803.06 << "\t" << fit803->GetParameter(1) << "\t" << fit803->GetParError(1) << "\t" << fit803->GetParameter(2) << "\t" << fit803->GetParError(2) << endl;
    
    //Bi214 1238 keV peak
    peak_c = 1238.122*m_c + b_c;
    peak_s = 1238.122*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit1238 = new TF1("fit1238",gaus_linear,xlow,xhigh,5);
    fit1238->SetParNames("A","#mu","#sigma","b","m");
    fit1238->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1238","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1238.122 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit1238->GetParameter(1));
    channel_e.push_back(fit1238->GetParError(1));
    sigma.push_back(fit1238->GetParameter(2));
    sigma_e.push_back(fit1238->GetParError(2));
    energy.push_back(1238.122);
    output << 1238.122 << "\t" << fit1238->GetParameter(1) << "\t" << fit1238->GetParError(1) << "\t" << fit1238->GetParameter(2) << "\t" << fit1238->GetParError(2) << endl;
    
    //1377 and 1385 kev peak from Bi214
    peak_c = 1377.669*m_c + b_c;
    peak_s = 1377.669*m_s + b_s;
    peak_c2 = 1385.31*m_c + b_c;
    peak_s2 = 1385.31*m_s + b_s;
    xlow = (peak_c + peak_c2)/2. - 30.;
    xhigh = (peak_c + peak_c2)/2. + 30.;
    TF1* fit1377_1385 = new TF1("fit1377_1385",double_gaus_linear,xlow,xhigh,8);
    fit1377_1385->SetParNames("A","#mu","#sigma","b","m","A_2","#mu_2","#sigma_2");
    fit1377_1385->SetParameters(9e5,peak_c,peak_s,10.,0.1,7e4,peak_c2,peak_s2);
    //fit1377_1385->FixParameter(1,peak_c);
    //fit238_240->FixParameter(2,peak_s);
    //fit1377_1385->FixParameter(6,peak_c2);
    //fit238_240->FixParameter(7,peak_s2);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1377_1385","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1377.669keV , 1385.31 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1377_1385->GetParameter(1));
    channel_e.push_back(fit1377_1385->GetParError(1));
    channel.push_back(fit1377_1385->GetParameter(6));
    channel_e.push_back(fit1377_1385->GetParError(6));
    sigma.push_back(fit1377_1385->GetParameter(2));
    sigma_e.push_back(fit1377_1385->GetParError(2));
    sigma.push_back(fit1377_1385->GetParameter(7));
    sigma_e.push_back(fit1377_1385->GetParError(7));
    energy.push_back(1377.669);
    energy.push_back(1385.31);
    output << 1377.669 << "\t" << fit1377_1385->GetParameter(1) << "\t" << fit1377_1385->GetParError(1) << "\t" << fit1377_1385->GetParameter(2) << "\t" << fit1377_1385->GetParError(2) << endl;
    output << 1385.31 << "\t" << fit1377_1385->GetParameter(6) << "\t" << fit1377_1385->GetParError(6) << "\t" << fit1377_1385->GetParameter(7) << "\t" << fit1377_1385->GetParError(7) << endl;
    
    //1401 and 1407 kev peak from Bi214
    peak_c = 1401.515*m_c + b_c;
    peak_s = 1401.515*m_s + b_s;
    peak_c2 = 1407.988*m_c + b_c;
    peak_s2 = 1407.988*m_s + b_s;
    xlow = (peak_c + peak_c2)/2. - 30.;
    xhigh = (peak_c + peak_c2)/2. + 30.;
    TF1* fit1401_1407 = new TF1("fit1401_1407",double_gaus_linear,xlow,xhigh,8);
    fit1401_1407->SetParNames("A","#mu","#sigma","b","m","A_2","#mu_2","#sigma_2");
    fit1401_1407->SetParameters(3e3,peak_c,peak_s,10.,0.1,6e3,peak_c2,peak_s2);
    //fit1401_1407->FixParameter(1,peak_c);
    //fit238_240->FixParameter(2,peak_s);
    //fit1401_1407->FixParameter(6,peak_c2);
    //fit238_240->FixParameter(7,peak_s2);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1401_1407","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1401.551 keV , 1407.988 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1401_1407->GetParameter(1));
    channel_e.push_back(fit1401_1407->GetParError(1));
    channel.push_back(fit1401_1407->GetParameter(6));
    channel_e.push_back(fit1401_1407->GetParError(6));
    sigma.push_back(fit1401_1407->GetParameter(2));
    sigma_e.push_back(fit1401_1407->GetParError(2));
    sigma.push_back(fit1401_1407->GetParameter(7));
    sigma_e.push_back(fit1401_1407->GetParError(7));
    energy.push_back(1401.515);
    energy.push_back(1407.988);
    output << 1401.515 << "\t" << fit1401_1407->GetParameter(1) << "\t" << fit1401_1407->GetParError(1) << "\t" << fit1401_1407->GetParameter(2) << "\t" << fit1401_1407->GetParError(2) << endl;
    output << 1407.988 << "\t" << fit1401_1407->GetParameter(6) << "\t" << fit1401_1407->GetParError(6) << "\t" << fit1401_1407->GetParameter(7) << "\t" << fit1401_1407->GetParError(7) << endl;
    
    //Bi214 1509 keV peak
    peak_c = 1509.210*m_c + b_c;
    peak_s = 1509.210*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit1509 = new TF1("fit1509",gaus_linear,xlow,xhigh,5);
    fit1509->SetParNames("A","#mu","#sigma","b","m");
    fit1509->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1509","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1509.210 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit1509->GetParameter(1));
    channel_e.push_back(fit1509->GetParError(1));
    sigma.push_back(fit1509->GetParameter(2));
    sigma_e.push_back(fit1509->GetParError(2));
    energy.push_back(1509.210);
    output << 1509.210 << "\t" << fit1509->GetParameter(1) << "\t" << fit1509->GetParError(1) << "\t" << fit1509->GetParameter(2) << "\t" << fit1509->GetParError(2) << endl;
    
    //Bi214 1729 keV peak
    peak_c = 1729.595*m_c + b_c;
    peak_s = 1729.595*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit1729 = new TF1("fit1729",gaus_linear,xlow,xhigh,5);
    fit1729->SetParNames("A","#mu","#sigma","b","m");
    fit1729->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1729","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1729.595 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit1729->GetParameter(1));
    channel_e.push_back(fit1729->GetParError(1));
    sigma.push_back(fit1729->GetParameter(2));
    sigma_e.push_back(fit1729->GetParError(2));
    energy.push_back(1729.595);
    output << 1729.595 << "\t" << fit1729->GetParameter(1) << "\t" << fit1729->GetParError(1) << "\t" << fit1729->GetParameter(2) << "\t" << fit1729->GetParError(2) << endl;
    
    //Bi214 1764 keV peak
    peak_c = 1764.491*m_c + b_c;
    peak_s = 1764.491*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit1764 = new TF1("fit1764",gaus_linear,xlow,xhigh,5);
    fit1764->SetParNames("A","#mu","#sigma","b","m");
    fit1764->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1764","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1764.491 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit1764->GetParameter(1));
    channel_e.push_back(fit1764->GetParError(1));
    sigma.push_back(fit1764->GetParameter(2));
    sigma_e.push_back(fit1764->GetParError(2));
    energy.push_back(1764.491);
    output << 1764.491 << "\t" << fit1764->GetParameter(1) << "\t" << fit1764->GetParError(1) << "\t" << fit1764->GetParameter(2) << "\t" << fit1764->GetParError(2) << endl;
    /*
    //1838 and 1847 kev peak from Bi214
    peak_c = 1838.36*m_c + b_c;
    peak_s = 1838.36*m_s + b_s;
    peak_c2 = 1847.429*m_c + b_c;
    peak_s2 = 1847.429*m_s + b_s;
    xlow = (peak_c + peak_c2)/2. - 50.;
    xhigh = (peak_c + peak_c2)/2. + 50.;
    TF1* fit1838_1847 = new TF1("fit1838_1847",double_gaus_linear,xlow,xhigh,8);
    fit1838_1847->SetParNames("A","#mu","#sigma","b","m","A_2","#mu_2","#sigma_2");
    fit1838_1847->SetParameters(1e2,peak_c,peak_s,10.,0.1,1e3,peak_c2,peak_s2);
    fit1838_1847->FixParameter(1,peak_c);
    //fit238_240->FixParameter(2,peak_s);
    fit1838_1847->FixParameter(6,peak_c2);
    //fit238_240->FixParameter(7,peak_s2);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit1838_1847","E","",xlow,xhigh);
    h1->SetTitle("Bi214 1838.36 keV , 1847.429 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit295->GetParameter(1) << endl;
    channel.push_back(fit1838_1847->GetParameter(1));
    channel_e.push_back(fit1838_1847->GetParError(1));
    channel.push_back(fit1838_1847->GetParameter(6));
    channel_e.push_back(fit1838_1847->GetParError(6));
    sigma.push_back(fit1838_1847->GetParameter(2));
    sigma_e.push_back(fit1838_1847->GetParError(2));
    sigma.push_back(fit1838_1847->GetParameter(7));
    sigma_e.push_back(fit1838_1847->GetParError(7));
    energy.push_back(1838.36);
    energy.push_back(1847.429);
    output << 1838.36 << "\t" << fit1838_1847->GetParameter(1) << "\t" << fit1838_1847->GetParError(1) << "\t" << fit1838_1847->GetParameter(2) << "\t" << fit1838_1847->GetParError(2) << endl;
    output << 1847.429 << "\t" << fit1838_1847->GetParameter(6) << "\t" << fit1838_1847->GetParError(6) << "\t" << fit1838_1847->GetParameter(7) << "\t" << fit1838_1847->GetParError(7) << endl;
    */
    //Bi214 2204.059 keV peak
    peak_c = 2204.059*m_c + b_c;
    peak_s = 2204.059*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit2204 = new TF1("fit2204",gaus_linear,xlow,xhigh,5);
    fit2204->SetParNames("A","#mu","#sigma","b","m");
    fit2204->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit2204","E","",xlow,xhigh);
    h1->SetTitle("Bi214 2204.059 keV peak");
    h1->GetXaxis()->SetRange(xlow,xhigh);
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit2204->GetParameter(1));
    channel_e.push_back(fit2204->GetParError(1));
    sigma.push_back(fit2204->GetParameter(2));
    sigma_e.push_back(fit2204->GetParError(2));
    energy.push_back(2204.059);
    output << 2204.059 << "\t" << fit2204->GetParameter(1) << "\t" << fit2204->GetParError(1) << "\t" << fit2204->GetParameter(2) << "\t" << fit2204->GetParError(2) << endl;
    
  }
  else if(calibSource == "AmBe"){
    // attempt to fit the 4438 keV peak from the strong AmBe source 
    peak_c = 4438.0*m_c + b_c;
    peak_s = 4438.0*m_s + b_s;
    xlow = peak_c - 20.;
    xhigh = peak_c + 20.;
    TF1* fit4438 = new TF1("fit4438",gaus_linear,xlow,xhigh,5);
    fit4438->SetParNames("A","#mu","#sigma","b","m");
    fit4438->SetParameters(9e4,peak_c,peak_s,247,0.099);
    h1->GetXaxis()->SetRange(0,0);
    h1->Fit("fit4438","E","",xlow,xhigh);
    h1->SetTitle("AmBe 4438 keV peak");
    h1->GetXaxis()->SetRange(2500,5500);
    //c1->SetLogy();
    c1->Update();
    c1->Print(pdfName.Data());
    //cout << fit242->GetParameter(1) << endl;
    channel.push_back(fit4438->GetParameter(1));
    channel_e.push_back(fit4438->GetParError(1));
    sigma.push_back(fit4438->GetParameter(2));
    sigma_e.push_back(fit4438->GetParError(2));
    energy.push_back(4438.0);
    output << 4438 << "\t" << fit4438->GetParameter(1) << "\t" << fit4438->GetParError(1) << "\t" << fit4438->GetParameter(2) << "\t" << fit4438->GetParError(2) << endl;
  }
  else if(calibSource == "print"){
    // print the full spectrum, no calibration, no fit 
    // the last two pages will be blank, but the first 2 pages will have the full spectrum
    // ideal? no but  ¯\_(ツ)_/¯
    h1->GetXaxis()->SetRange(0,0);
    h1->SetTitle(("Full spectrum - " + inData).c_str());
    h1->GetXaxis()->SetRange(0,5500);
    h1->GetXaxis()->SetTitle("Channel");
    h1->Draw();
    h1->Write();
    c1->Update();
    c1->Print(pdfName.Data());
    // print with log scale too 
    c1->SetLogy();
    c1->Update();
    c1->Print(pdfName.Data());
  }


  int n = channel.size();
  //Double_t energy[6] = {727.33,238.632,240.986,300.087,510.770,583.191,860.564};//,2614.533};
  TGraphErrors *relation = new TGraphErrors(n,&energy[0],&channel[0],0,&channel_e[0]);
  relation->SetTitle("ADC channel vs Energy");
  relation->GetYaxis()->SetTitle("ADC [channel]");
  relation->GetXaxis()->SetTitle("Energy [keV]");
  relation->Draw("A*");
  relation->Fit("pol1");
  relation->SetTitle("ADC vs Energy relation");
  c1->Update();
  relation->Write();
  c1->Print(pdfName.Data());
  
  TGraphErrors *relation2 = new TGraphErrors(n,&energy[0],&sigma[0],0,&sigma_e[0]);
  relation2->SetTitle("Sigma vs Energy");
  relation2->GetYaxis()->SetTitle("Sigma [Channel]");
  relation2->GetXaxis()->SetTitle("Energy [keV]");
  relation2->Draw("A*");
  relation2->Fit("pol1");
  relation2->SetTitle("Resolution: Sigma(channel) vs Energy(KeV)");
  c1->Update();
  c1->Print(pdfName.Data());
  TString pdfFileClose = pdfName + "]";
  c1->Print(pdfFileClose.Data());
  h1->Draw();
  relation2->Write();
  auto c2 = new TCanvas("c2","c2",width,hight);
  h2->Draw();
  h2->Write();

  file->Close();

}
