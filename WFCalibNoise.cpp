//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.cpp"

void WFCalibNoise() {
   const Int_t N = 256;
   Double_t canvX = 1150., canvY = 650;
   Int_t nBin = 100;
   
   TRandom3 rndm(0);
   
   Double_t lambda = 0;
   Double_t alpha = 0;
   Double_t alphaCT = 0;
   Double_t dcr = 0;
   Double_t noise = 0;
   
   Int_t nEvent = 100;
   Int_t nPhoton = 100;
   Int_t nRepeat = 1;
//   Int_t nRepeat = 5;
   const Int_t nOV = 9;
   Int_t nPar = 4;

/*   Double_t maxMean[nOV] = {8,8,8,8,8,
                            8,8,8,8};
   Double_t maxVar[nOV] = {8e-3,8e-3,8e-3,8e-3,8e-3,
                           1e-2,1e-2,2e-2,2e-2};
//   Double_t maxVar[nOV] = {2e-3, 5e-3, 1e-2, 2e-2, 3e-2,
//                           5e-2, 1e-1, 1e-1, 2e-1};
*/
   TCanvas *canvDiffVar = new TCanvas("canvDiffVar", "DiffVar", canvX, canvY);
   canvDiffVar->Divide(3,3);
   
   TClonesArray *caPar = new TClonesArray("TGraphErrors", 4);
   for(Int_t iPar=0; iPar<nPar; iPar++) {
      new ((TGraphErrors*)(*caPar)[iPar]) TGraphErrors;
   }
   
   for(Int_t iOV=0; iOV<nOV; iOV++) {
      for(Int_t iRepeat=0; iRepeat<nRepeat; iRepeat++) {
      Double_t trueGain = .03;
      noise = .0001*iOV;
      cout << "noise: " << noise << endl;
      
      Double_t charge[nEvent*nPhoton], diffVar[nEvent*nPhoton];
      Double_t maxCharge=-1e9, maxVar = -1e9;

      for(Int_t iPhoton=0; iPhoton<nPhoton; iPhoton++) {
         for(int iEvent=0; iEvent<nEvent; iEvent++){
            Waveform* wf = new Waveform;
            wf->SetTrueGain(trueGain);
            wf->SetTauTiming(10);
            wf->MakeTemplate();
            wf->SetNoises(lambda, alpha, alphaCT, dcr, noise);
            wf->MakeEvent(iPhoton);
            charge[iPhoton*nEvent+iEvent] = wf->GetCharge();
            wf->Differentiate(2);
            diffVar[iPhoton*nEvent+iEvent] = wf->GetVariance();
            if(charge[iPhoton*nEvent+iEvent]>maxCharge)
               maxCharge = charge[iPhoton*nEvent+iEvent];
            if(diffVar[iPhoton*nEvent+iEvent]>maxVar)
               maxVar = diffVar[iPhoton*nEvent+iEvent];
         }
      }

      TH2D *histDiffVar = 
            new TH2D("histDiffVar", Form("noise: %.4f", noise),
                     nBin, 0, maxCharge, nBin, 0, maxVar);
//                     nBin, 0, maxMean[iOV], nBin, 0, maxVar[iOV]);

      for(Int_t iPhoton=0; iPhoton<nPhoton; iPhoton++) {
         for(int iEvent=0; iEvent<nEvent; iEvent++){
            histDiffVar->Fill(charge[iPhoton*nEvent+iEvent],
                              diffVar[iPhoton*nEvent+iEvent]);
         }
      }
            
      canvDiffVar->cd(iOV+1)->SetLogz();
      histDiffVar->Draw("colz");
      histDiffVar->FitSlicesY();
      TH1D *histMean;
      histMean = (TH1D*)gDirectory->Get("histDiffVar_1");
      histMean->Draw("same");

      TF1* funcQuad = new TF1("funcQuad", "pol2");
      funcQuad->SetLineColor(kRed);
      funcQuad->SetParameters(0, .001, .1);
      histMean->Fit(funcQuad, "Q");
      
      Double_t parDiffVar[3];
      Double_t errDiffVar[3];
      for (int iPar = 0; iPar < 3; iPar++) {
         parDiffVar[iPar]=funcQuad->GetParameter(iPar);
         errDiffVar[iPar]=funcQuad->GetParError(iPar);
         
         ((TGraphErrors*)(*caPar)[iPar])->
               SetPoint(iOV*nRepeat+iRepeat, noise, parDiffVar[iPar]);
         ((TGraphErrors*)(*caPar)[iPar])->
               SetPointError(iOV*nRepeat+iRepeat, 0, errDiffVar[iPar]);
      }
/*
      Double_t gainDiffVar = parDiffVar[1]/(1-parDiffVar[2]);
      Double_t gainDiffVarErr = 
            TMath::Sqrt(TMath::Power(errDiffVar[1]/(1-parDiffVar[2]), 2)
                        + TMath::Power(parDiffVar[1]*errDiffVar[2]
                                       /(1-parDiffVar[2])/(1-parDiffVar[2])
                                       , 2));
      ((TGraphErrors*)(*caPar)[3])->SetPoint(iOV, lambda, gainDiffVar);
      ((TGraphErrors*)(*caPar)[3])->SetPointError(iOV, 0, gainDiffVarErr);
*/
   }
   }
   
   TCanvas *canvPar = new TCanvas("canvPar", "Parameters", canvX, canvY);
   canvPar->Divide(2,2);
   
   ((TGraphErrors*)(*caPar)[0])->SetTitle("C0;noise;C0");
   ((TGraphErrors*)(*caPar)[1])->SetTitle("C1;noise;C1");
   ((TGraphErrors*)(*caPar)[2])->SetTitle("C2;noise;C2");
   ((TGraphErrors*)(*caPar)[3])->
         SetTitle("Calculated Gain;Gain;Calculated Gain");
   
   for(Int_t iPar=0; iPar<3; iPar++) {
      canvPar->cd(iPar+1);
      ((TGraphErrors*)(*caPar)[iPar])->Draw("ap");
      ((TGraphErrors*)(*caPar)[iPar])->SetMarkerStyle(8);
      ((TGraphErrors*)(*caPar)[iPar])->SetMarkerSize(.5);
      ((TGraphErrors*)(*caPar)[iPar])->GetXaxis()->SetLimits(0,.001);
      if(iPar==1) ((TGraphErrors*)(*caPar)[iPar])->SetMinimum(0);
   }
}
