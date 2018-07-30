#include "Waveform.h"

void Init(void);
void Event(Int_t nPhe, Waveform* wf);
// Double_t YCut(TGraph* gr);
void YCut(TGraph* gr,Double_t *par);
Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err);
Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err);
// void Differentiate(Double_t *wf,Int_t ndiff);
Double_t fitfunc(Double_t *x, Double_t *par);


TGraph* grQNQWF=new TGraph();
Double_t grQNQange[2] = {-100, 1000};

Double_t Gain=1;
Int_t gNNoise=1;
Int_t Nrep=5;
std::vector<Double_t>noiselist={0,0.05,0.1};
// Double_t noiselevel=0;
// Double_t WFRaw[gNbin] = {};
// Double_t WFDiff[gNbin] = {};
// Double_t WFNoise[gNNoise][gNbin] = {};
// Double_t WFSum[gNbin] = {};

// Double_t rawWF[gNbin]  = {};
Double_t WFT[gNbin]    = {};

void Init(void) {
   gFSctime->SetParameter(0,45);
   gFSingle->SetParameter(0,20);
   gFSingle->SetParameter(1,0);
   gFSingle->FixParameter(2,Gain);
   // for (int iBin = 0; iBin < gNbin; iBin++) WFT[iBin] = (grQNQange[1] - grQNQange[0]) /gNbin * iBin + grQNQange[0];
}

void StatGainCalib(void) {
   // gNNoise=randnoise.size();
   TCanvas *cgr0 =new TCanvas("cgr0", "cgr0",600,600);
   TCanvas *cgr1 =new TCanvas("cgr1", "cgr1",600,600);
   // TCanvas *cgrQNvar =new TCanvas("cgrQNvar", "cgrQNvar",600,600);
   Init();
   std::vector<Double_t> vecSum;
   std::vector<std::vector<Double_t>> vecNoise;
   std::vector<std::vector<Double_t>> vecInterference;

   Double_t sum;
   // Double_t noise;

   Int_t npheRange[2] = {1,3000};
   Int_t npheStep = 100;

   Double_t logstep = (TMath::Log(npheRange[1]) - TMath::Log(npheRange[0]))/ (npheStep -1);


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


   for (int istep = 0; istep < npheStep; istep++) {
      // int iphe = (int)TMath::Exp(TMath::Log(npheRange[0]) + logstep * istep);
      int iphe = istep*(npheRange[1]-npheRange[0])/npheStep;
      for(int irep=0;irep<Nrep;irep++){
         std::vector<Double_t> noisevar;
         std::vector<Double_t> inter;
         cout<<istep<<" "<<iphe<<endl;

         for (int i = 0; i < noiselist.size(); i++) {
            Waveform* wf= new Waveform(gNbin,-100,1000);
            wf->MakeEvent(iphe);

            // if (i==0) {
            // vecSum.push_back();
            // // std::cout<<"iphe: "<<iphe<<" charge : "<<wf->GetCharge()<<std::endl;
            // }

            wf->SetNoiseLevel(noiselist[i]);
            Double_t charge=wf->GetCharge();
            wf->Differentiate(2);
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


   // for(int i=0;i<vecSum.size();i++){
   //    // for(int iNoise=0;iNoise<gNNoise;iNoise++){
   //    for (int j = 0; j < noiselist.size(); j++) {
   //       ((TGraphErrors*)(*cagrQNvar)[j])->SetPoint(i,vecSum[i],vecNoise[i][j]);
   //               ((TGraphErrors*)(*cagrQNvar)[j])->SetPointError(i,0,vecNoise[i][j]*0.2);
   //       // ((TGraph*)(*cagrInt)[j])->SetPoint(i,vecSum[i],vecInterference[i][j]);
   //       // grQNQDiff[iNoise]->SetPoint(i,vecSum[i],vecDiff[iNoise][i]);
   //       // }
   //    }
   // }

   TString fname="gain1_noise5e-2.pdf";

   TF1* func= new TF1("fitfunc",fitfunc,-10,5000,3);

   cgr0->cd();
   // std::vector<Double_t> noff;
   // std::vector<Double_t> scale;
   for (int i = 0; i < noiselist.size(); i++) {

      // ((TGraph*)(*cagrInt)[i])->SetMaximum(1000);
      // ((TGraph*)(*cagrInt)[i])->SetMinimum(0);
      // ((TGraph*)(*cagrInt)[i])->SetMarkerStyle(20);
      // ((TGraph*)(*cagrInt)[i])->SetMarkerColor(2+i);
      // ((TGraph*)(*cagrInt)[i])->Draw("p same");
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
      // return pol1->GetParameter(0);
      // YCut(((TGraph*)(*cagrQNvar)[i]),par);
      // noff.push_back(par[0]);
      // scale.push_back(par[1]);
      // std::cout<<"Ycut: "<<par[0]<<std::endl;
      // std::cout<<"scale: "<<par[1]<<std::endl;
      // Double_t coeff=par[2];
      // Double_t pm=coeff*Gain/(1+coeff*Gain);
      // // Double_t Measuredpm=2.4e-4;
      // Double_t Measuredpm=1.6e-4;
      // Double_t MeasuredGain=Measuredpm/coeff/(1-Measuredpm);
      // std::cout<<"p_{multiple}: "<<pm<<std::endl;
      // std::cout<<"Gain: "<<MeasuredGain<<std::endl;
      Double_t MesGain=par[1]/(1-par[1]*par[2]);
      Double_t pmul=par[1]*par[2];
      Double_t pmulerr=MultipleError(par[1],err[1],par[2],err[2]);
      Double_t MesGainerr=DivisionError(par[1],err[1],(double)1-pmul,pmulerr);
      std::cout<<"Gain: "<<MesGain<<"+-"<<MesGainerr<<std::endl;

      std::cout<<"Pmul: "<<pmul<<"+-"<<pmulerr<<std::endl;
      // ((TGraph*)(*cagrQNvar)[i])->Fit("fitfunc","","");
      ((TGraphErrors*)(*cagrQNvar)[i])->SetTitle("Q_{dint}vs Q^{2}_{drms};Q_{dint};Q^{2}_{drms}");
      ((TGraph*)(*cagrQNvar)[i])->SetMaximum(100);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMinimum(0);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerStyle(20);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerColor(2+i);
      if (i==0) {
         ((TGraphErrors*)(*cagrQNvar)[i])->Draw("ap");
      }else{
         ((TGraphErrors*)(*cagrQNvar)[i])->Draw("p same");
      }
   }

   // Double_t pm=func->GetParameter(2);
   //
   // for(int i=0;i<vecSum.size();i++){
   //    for (int j = 0; j < noiselist.size(); j++) {
   //       // std::cout<<"point set"<<std::endl;
   //       if(vecSum[i]>0){
   //          // ((TH2D*)(*cahQNQ)[j])->Fill(vecSum[i],(vecNoise[i][j]-noff[j])/vecSum[i]);
   //          ((TGraph*)(*cagrQNQ)[j])->SetPoint(((TGraph*)(*cagrQNQ)[j])->GetN(),vecSum[i],(vecNoise[i][j]-noff[j])/vecSum[i]);
   //       }
   //    }
   // }
   //
   // cgr1->cd();
   // cgr1->Divide(2,2);
   // for (int i = 0; i < noiselist.size(); i++) {
   //
   //    // cgr1->cd(i+1);
   //    // ((TH2D*)(*cahQNQ)[i])->Draw("colz");
   //    // ((TH2D*)(*cahQNQ)[i])->FitSlicesY();
   //    //
   //    // TH1D* hFSQNQ=(TH1D*)gDirectory->Get(Form("hQNQ%d_1",i));
   //    // hFSQNQ->SetMarkerStyle(20);
   //    // hFSQNQ->SetMarkerColor(2);
   //    // hFSQNQ->Draw("same p");
   //    // hFSQNQ->Fit("fpol1");
   //    ((TGraph*)(*cagrQNQ)[i])->SetTitle("Q_{dint} vs Q^{2}_{drms}/Q_{dint};Q_{dint};Q^{2}_{drms}/Q_{dint}");
   //    ((TGraph*)(*cagrQNQ)[i])->SetMaximum(5);
   //    ((TGraph*)(*cagrQNQ)[i])->SetMinimum(0);
   //    ((TGraph*)(*cagrQNQ)[i])->SetMarkerStyle(20);
   //    ((TGraph*)(*cagrQNQ)[i])->SetMarkerColor(2+i);
   //    // ((TGraph*)(*cagrQNQ)[i])->Draw("ap");
   //    ((TGraph*)(*cagrQNQ)[i])->Fit(Form("fpol1%d",i),"Q","",0,10000);
   //    if (i==0) {
   //       ((TGraph*)(*cagrQNQ)[i])->Draw("ap");
   //    }else{
   //       ((TGraph*)(*cagrQNQ)[i])->Draw("p same");
   //    }
   //    // fpol1->SetLineColor()
   //
   //    // fpol1->Draw("same");
   //
   // }

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

void YCut(TGraph* gr,Double_t *par){
   // TF1* pol1= new TF1("fitfunc",fitfunc,-100,5000,3);
   // // pol1->SetParLimits(2,0,1);
   // // pol1->SetParLimits(0,-1,20);
   // gr->Fit("fitfunc","","");
   // for (int i = 0; i < 3; i++) {
   //    par[i]=pol1->GetParameter(i);
   // }
   // // return pol1->GetParameter(0);
}

Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err){
   /*calculate error of value 1/ value 2*/
   return TMath::Sqrt(TMath::Power(value1err/value2,2)+TMath::Power(value1*value2err/(value2*value2),2));
}

Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err){
   /*calculate error of value 1* value 2*/
   return TMath::Sqrt(TMath::Power(value1err* value2,2)+TMath::Power(value1*value2err,2));
}