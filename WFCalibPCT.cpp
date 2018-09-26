//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.cpp"

void WFCalibPCT() {
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
   const Int_t nOV = 9;
   Int_t nPar = 4;

   Double_t maxMean[nOV] = {1.5, 1.5, 1.5, 2, 2,
                            2, 2, 3, 3};
   Double_t maxVar[nOV] = {1.5e-3, 1.5e-3, 2e-3, 3e-3, 4e-3,
                           5e-3, 1e-2, 1e-2, 1.5e-2};

   TCanvas *canvDiffVar = new TCanvas("canvDiffVar", "DiffVar", canvX, canvY);
   canvDiffVar->Divide(3,3);
   
   TClonesArray *caMean = new TClonesArray("TGraphErrors", nOV);
   TClonesArray *histNCT = new TClonesArray("TH1D", 4);
   TClonesArray *caPar = new TClonesArray("TGraphErrors", 4);
   for(Int_t iPar=0; iPar<nPar; iPar++) {
      new ((TGraphErrors*)(*caPar)[iPar]) TGraphErrors;
      new ((TH1D*)(*histNCT)[iPar])
            TH1D(Form("hist%d", iPar),
                 Form("%.3f < height squared < %.3f", .001*iPar, .001*(iPar+1)),
                 30, 0, 30);
   }
   
   for(Int_t iOV=0; iOV<nOV; iOV++) {
//      if(iOV<8) continue;
      new ((TGraphErrors*)(*caMean)[iOV]) TGraphErrors;
      for(Int_t iRepeat=0; iRepeat<nRepeat; iRepeat++) {
      Double_t trueGain = .03;
//      noise = .0005*iOV;
      lambda = .05*iOV;
      cout << "lambda: " << lambda << endl;
      
      TH2D *histDiffVar = 
            new TH2D("histDiffVar", Form("lambda: %.2f", lambda),
                     nBin, 0, maxMean[iOV], nBin, 0, maxVar[iOV]);

      for(Int_t iPhoton=0; iPhoton<nPhoton; iPhoton++) {
         for(int iEvent=0; iEvent<nEvent; iEvent++){
            Waveform* wf = new Waveform;
            wf->SetTrueGain(trueGain);
            wf->SetTauTiming(10);
            wf->MakeTemplate();
            wf->SetNoises(lambda, alpha, alphaCT, dcr, noise);
            Int_t nHit = wf->MakeEvent(iPhoton);
            Int_t nCT = wf->GetNCT();
            Double_t charge = wf->GetCharge();
            wf->Differentiate(2);
            Double_t diffVar = wf->GetVariance();
            histDiffVar->Fill(charge, diffVar);
            if(iOV==8 && charge>.5 && charge<1.5) {
               for(Int_t iDV=0; iDV<4; iDV++) {
                  if(diffVar>iDV*0.001 && diffVar<(iDV+1)*0.001) {
                     ((TH1*)(*histNCT)[iDV])->Fill(nCT);
//                     cout << "nHit: " << nHit << ", nCT: " << nCT << endl;
                  }
               }
            }
         }
      }

      canvDiffVar->cd(iOV+1)->SetLogz();
      histDiffVar->Draw("colz");
//      histDiffVar->FitSlicesY();
//      histMean = (TH1D*)gDirectory->Get("histDiffVar_1");
      for(Int_t iBin=0; iBin<nBin; iBin++) {
         TH1D *histMean = histDiffVar->ProjectionY("histMean", iBin+1, iBin+1);
         Double_t mean = histMean->GetXaxis()->
               GetBinCenter(histMean->GetMaximumBin());
         Double_t stdDev = histMean->GetStdDev();
         TF1* funcGauss = 
               new TF1("funcGauss", "gausn", mean-stdDev, mean+stdDev);
         histMean->Fit(funcGauss, "NQ", "", mean-stdDev, mean+stdDev);
         mean = funcGauss->GetParameter(1);
         stdDev = TMath::Abs(funcGauss->GetParameter(2));
         Double_t meanErr = funcGauss->GetParError(1);
         Double_t stdDevErr = funcGauss->GetParError(2);
         ((TGraphErrors*)(*caMean)[iOV])->SetPoint(iBin, mean, stdDev*stdDev);
         ((TGraphErrors*)(*caMean)[iOV])->
               SetPointError(iBin, meanErr, 2*stdDev*stdDevErr);
      }

      TF1* funcQuad = new TF1("funcQuad", "pol2");
      funcQuad->SetLineColor(kRed);
      ((TGraphErrors*)(*caMean)[iOV])->Fit(funcQuad, "Q");
      
      Double_t parDiffVar[3];
      Double_t errDiffVar[3];
      for (int iPar = 0; iPar < 3; iPar++) {
         parDiffVar[iPar]=funcQuad->GetParameter(iPar);
         errDiffVar[iPar]=funcQuad->GetParError(iPar);
         
         ((TGraphErrors*)(*caPar)[iPar])->
               SetPoint(iOV*nRepeat+iRepeat, lambda, parDiffVar[iPar]);
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

   TCanvas* canvNCT = new TCanvas("canvNCT", "# of CT", canvX, canvY);
      
   for(Int_t iDV=0; iDV<4; iDV++) {
      if(iDV==0) {
         ((TH1*)(*histNCT)[iDV])->Draw();
      } else {
         ((TH1*)(*histNCT)[iDV])->Draw("same");
      }
      ((TH1*)(*histNCT)[iDV])->SetLineColor(iDV+2);
   }
   
   TCanvas *canvPar = new TCanvas("canvPar", "Parameters", canvX, canvY);
   canvPar->Divide(2,2);
   
   ((TGraphErrors*)(*caPar)[0])->SetTitle("C0;lambda;C0");
   ((TGraphErrors*)(*caPar)[1])->SetTitle("C1;lambda;C1");
   ((TGraphErrors*)(*caPar)[2])->SetTitle("C2;lambda;C2");
   ((TGraphErrors*)(*caPar)[3])->
         SetTitle("Calculated Gain;Gain;Calculated Gain");
   
   for(Int_t iPar=0; iPar<3; iPar++) {
      canvPar->cd(iPar+1);
      ((TGraphErrors*)(*caPar)[iPar])->Draw("ap");
      ((TGraphErrors*)(*caPar)[iPar])->SetMarkerStyle(8);
      ((TGraphErrors*)(*caPar)[iPar])->SetMarkerSize(.5);
      ((TGraphErrors*)(*caPar)[iPar])->GetXaxis()->SetLimits(0,.5);
      if(iPar==1) ((TGraphErrors*)(*caPar)[iPar])->SetMinimum(0);
   }
}
