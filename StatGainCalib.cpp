#include "Waveform.h"

void Init(void);
void Event(Int_t nPhe, Waveform* wf);
// Double_t YCut(TGraph* gr);
// void YCut(TGraph* gr,Double_t *par);
void LoadSimConfig();
Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err);
Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err);
// void Differentiate(Double_t *wf,Int_t ndiff);
Double_t fitfunc(Double_t *x, Double_t *par);


TGraph* grQNQWF=new TGraph();
Double_t grQNQange[2] = {-100, 1000};

Int_t Nrep=5;

Double_t WFT[gNbin]    = {};

void Init(void) {
   LoadSimConfig();
   gFSctime->SetParameter(0,ScintDecay);
   gFSingle->SetParameter(0,SPwidth);
   gFSingle->SetParameter(1,0);
   gFSingle->FixParameter(2,Gain);
   gAPtime ->FixParameter(0,APtimeconstant);
   // for (int iBin = 0; iBin < gNbin; iBin++) WFT[iBin] = (grQNQange[1] - grQNQange[0]) /gNbin * iBin + grQNQange[0];
}

void StatGainCalib(void) {
   TCanvas *cgr0 =new TCanvas("cgr0", "cgr0",600,600);
   TCanvas *cgr1 =new TCanvas("cgr1", "cgr1",600,600);
   // TCanvas *cgrQNvar =new TCanvas("cgrQNvar", "cgrQNvar",600,600);
   Init();
   std::vector<Double_t> vecSum;
   std::vector<std::vector<Double_t>> vecNoise;
   std::vector<std::vector<Double_t>> vecInterference;

   Double_t sum;

   Int_t npheRange[2] = {RangeMin,RangeMax};

   Double_t logstep = (TMath::Log(npheRange[1]) - TMath::Log(npheRange[0]))/ (Nstep -1);


   TClonesArray* cagrQNvar = new TClonesArray("TGraphErrors",noiselist.size());
   TClonesArray* cagrQNQ   = new TClonesArray("TGraph",noiselist.size());
   // TClonesArray* cagrInt   = new TClonesArray("TGraph",noiselist.size());
   TClonesArray* cahQNQ    = new TClonesArray("TH2D",noiselist.size());
   TClonesArray* cafpol1    = new TClonesArray("TF1",noiselist.size());
   // TF1* fpol1= new TF1("fpol1","[0]*x+[1]");
   for (int i = 0; i < noiselist.size(); i++) {
      new ((TGraphErrors*)(*cagrQNvar)[i]) TGraphErrors();
      new ((TGraph*)(*cagrQNQ)[i])   TGraph();
      // new ((TGraph*)(*cagrInt)[i])   TGraph();
      // new ((TH2D*)(*cahQNQ)[i])   TH2D(Form("hQNQ%d",i),Form("hQNQ%d",i),100,0,5000,100,0,5);
      new ((TF1*)(*cafpol1)[i])   TF1(Form("fpol1%d",i),"[0]*x+[1]");
   }


   for (int istep = 0; istep < Nstep; istep++) {
      // int iphe = (int)TMath::Exp(TMath::Log(npheRange[0]) + logstep * istep);
      int iphe = istep*(npheRange[1]-npheRange[0])/Nstep;
      for(int irep=0;irep<Nrep;irep++){
         std::vector<Double_t> noisevar;
         std::vector<Double_t> inter;
         // cout<<istep<<" "<<iphe<<" sum : "<<sum<<endl;

         for (int i = 0; i < noiselist.size(); i++) {
            Waveform* wf= new Waveform(gNbin,-100,1000);
            wf->MakeEvent(iphe);

            wf->SetNoiseLevel(noiselist[i]);
            Double_t charge=wf->GetCharge();
            wf->Differentiate(Ndiff);
            Double_t noisevar=wf->GetTotalVariance();
            // std::cout<<"variance: "<<noisevar<<std::endl;
            // noisevar.push_back(wf->GetTotalVariance());
            // inter.push_back(wf->GetInterference());
            Int_t index=istep*Nrep+irep;
            ((TGraphErrors*)(*cagrQNvar)[i])->SetPoint(     index,charge,noisevar);
            ((TGraphErrors*)(*cagrQNvar)[i])->SetPointError(index,0     ,noisevar*0.2);
            delete wf;
         }
         vecNoise.push_back(noisevar);
         vecInterference.push_back(inter);
      }
   }

   TString fname="gain1_noise5e-2.pdf";

   TF1* func= new TF1("fitfunc",fitfunc,-10,5000,3);

   cgr0->cd();
   for (int i = 0; i < noiselist.size(); i++) {

      Double_t par[3];
      Double_t err[3];
      // Double_t YACut=YCut(((TGraph*)(*cagrQNvar)[i]));
      TF1* pol1= new TF1("fitfunc",fitfunc,-100,5000,3);
      // pol1->SetParLimits(2,0,1);
      // pol1->SetParLimits(0,-1,20);
      ((TGraphErrors*)(*cagrQNvar)[i])->Fit("fitfunc","","");
      // Double_t par[3];
      for (int i = 0; i < 3; i++) {
         par[i]=pol1->GetParameter(i);
         err[i]=pol1->GetParError(i);
      }

      Double_t MesGain=par[1]/(1-par[1]*par[2]);
      Double_t pmul=par[1]*par[2];
      Double_t pmulerr=MultipleError(par[1],err[1],par[2],err[2]);
      Double_t MesGainerr=DivisionError(par[1],err[1],(double)1-pmul,pmulerr);
      std::cout<<"Gain: "<<MesGain<<"+-"<<MesGainerr<<std::endl;

      std::cout<<"Pmul: "<<pmul<<"+-"<<pmulerr<<std::endl;
      // ((TGraph*)(*cagrQNvar)[i])->Fit("fitfunc","","");
      ((TGraphErrors*)(*cagrQNvar)[i])->SetTitle("Q_{dint}vs Q^{2}_{drms};Q_{dint};Q^{2}_{drms}");
      ((TGraph*)(*cagrQNvar)[i])->SetMaximum(500);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMinimum(0);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerStyle(20);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerColor(2+i);
      if (i==0) {
         ((TGraphErrors*)(*cagrQNvar)[i])->Draw("ap");
      }else{
         ((TGraphErrors*)(*cagrQNvar)[i])->Draw("p same");
      }
   }
}

Double_t fitfunc(Double_t *x, Double_t *par){
   //par[0]: offset;
   //par[1]:
   Double_t  xx        = x[0];
   // Double_t PMultiple = TMath::Sqrt(TMath::Power(par[2]*xx,2)/1+TMath::Power(par[2]*xx,2));
   // Double_t PMultiple = par[2];
   // Double_t Gain = 1;
   return par[0]+par[1]*(xx+par[2]*xx*xx);
}

Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err){
   /*calculate error of value 1/ value 2*/
   return TMath::Sqrt(TMath::Power(value1err/value2,2)+TMath::Power(value1*value2err/(value2*value2),2));
}

Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err){
   /*calculate error of value 1* value 2*/
   return TMath::Sqrt(TMath::Power(value1err* value2,2)+TMath::Power(value1*value2err,2));
}

void LoadSimConfig(){
   ifstream fsimconf;
   TString simconfname = "./SimConfig.dat";
   std::cout<<"Loading Configuration..."<<std::endl;
   fsimconf.open(simconfname.Data(),std::ios::in);
   string linestr;
   while (!fsimconf.eof()) {
      fsimconf>>linestr;
      std::cout<<linestr<<std::endl;
      // string str( "abcdefghijk" );
      char key = ':';
      Int_t colonpos=linestr.find(key);
      string strvalue = linestr.substr(colonpos+1);
      // std::cout << colonpos << std::endl;
      // std::cout <<  << std::endl;
      if (linestr.find(strLambda)!=string::npos) {
         lambda = std::stod(strvalue);
      }
      if (linestr.find(strAlpha)!=string::npos) {
         alpha = std::stod(strvalue);
      }
      if (linestr.find(strScintDecay)!=string::npos) {
         ScintDecay = std::stod(strvalue);
      }
      if (linestr.find(strSPwidth)!=string::npos) {
         SPwidth = std::stod(strvalue);
      }
      if (linestr.find(strAPtimeconstant)!=string::npos) {
         APtimeconstant = std::stod(strvalue);
      }
      if (linestr.find(strGain)!=string::npos) {
         Gain = std::stod(strvalue);
      }
      if (linestr.find(strRangeMin)!=string::npos) {
         RangeMin = std::stoi(strvalue);
      }
      if (linestr.find(strRangeMax)!=string::npos) {
         RangeMax = std::stoi(strvalue);
      }
      if (linestr.find(strNstep)!=string::npos) {
         Nstep = std::stoi(strvalue);
      }
      if (linestr.find(strNdiff)!=string::npos) {
         Ndiff = std::stoi(strvalue);
      }
      if (linestr.find(strNevent)!=string::npos) {
         Nevent = std::stoi(strvalue);
      }
      if (linestr.find(strNoiseLevel)!=string::npos) {
         noiselist.push_back(std::stod(strvalue));
      }
   }
}

