//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.cpp"

void ENF() {
   const Int_t N = 256;
   Double_t canvX = 1150., canvY = 650;
   Int_t nBin = 1000;
   Int_t nDCR = 4;
   
   TRandom3 rndm(0);

   Double_t trueGain = .03;
   Double_t lambda = 0.2;
   Double_t alpha = 0;
   Double_t alphaCT = 0;
   Double_t dcr = 0;
   Double_t noise = 0;
   
   const Int_t nEvent = 10000;
   Int_t nMu = 9;
   Int_t nRepeat = 1;

   TClonesArray* grENF = new TClonesArray("TGraphErrors", nDCR);
   TClonesArray* canvHist = new TClonesArray("TCanvas", nDCR);
   TLegend* leg = new TLegend(0.8, 0.68, 0.99, 0.78);
   leg->SetFillColor(0);

   for(Int_t iDCR=0; iDCR<nDCR; iDCR++) {
      new ((TGraphErrors*)(*grENF)[iDCR]) TGraphErrors;
      new ((TCanvas*)(*canvHist)[iDCR])
            TCanvas(Form("%d",iDCR), Form("%d",iDCR), canvX, canvY);
      ((TCanvas*)(*canvHist)[iDCR])->Divide(3,3);

      dcr = 1000*iDCR;
      cout << "DCR: " << dcr << endl;
      leg->AddEntry(((TGraphErrors*)(*grENF)[iDCR]),
                    Form("%.0f kHz", dcr), "pl");

      for(Int_t iRepeat=0; iRepeat<nRepeat; iRepeat++) {
         for(Int_t iMu=0; iMu<nMu; iMu++) {
            Double_t mu = 0.5*(iMu+1);
            Double_t charge[nEvent];
            Double_t max=-1e9, min=1e9;
            Int_t nPed = 0;
            for(int iEvent=0; iEvent<nEvent; iEvent++){
               Int_t nPhoton = rndm.Poisson(mu);
               Waveform* wf = new Waveform;
               wf->SetTrueGain(trueGain);
               wf->MakeTemplate();
               wf->SetNoises(lambda, alpha, alphaCT, dcr, noise);
               wf->MakeEvent(nPhoton);
               charge[iEvent] = wf->GetCharge();
               if(max<charge[iEvent]) max = charge[iEvent];
               if(min>charge[iEvent]) min = charge[iEvent];
               if(charge[iEvent]<.5*trueGain) nPed++;
            }
            TH1D* hist = new TH1D("hist", Form("Mu: %.1f", mu), nBin, min, max);
            for(Int_t iEvent=0; iEvent<nEvent; iEvent++) {
               hist->Fill(charge[iEvent]);
            }
            ((TCanvas*)(*canvHist)[iDCR])->cd(iMu+1);
            hist->Draw();
            Double_t mean = hist->GetMean();
            Double_t meanErr = hist->GetMean(11);
            Double_t stdDev = hist->GetStdDev();
            Double_t stdDevErr = hist->GetStdDev(11);
            Double_t f0 = (Double_t)nPed/nEvent;
            Double_t nPedErr = TMath::Sqrt(nEvent*f0*(1-f0));
            Double_t f0Err = nPedErr/nEvent;
            Double_t muCalc = -TMath::Log(f0);
            cout << "muCalc: " << muCalc << ", muTrue: " << mu << endl;
            Double_t muErr = f0Err/f0;
            Double_t enf = mu*stdDev*stdDev/mean/mean;
            Double_t enfErr = enf *
                  TMath::Sqrt(muErr*muErr/mu/mu +
                              stdDevErr*stdDevErr/stdDev/stdDev +
                              meanErr*meanErr/mean/mean);
            ((TGraphErrors*)(*grENF)[iDCR])->
                  SetPoint(iRepeat*nMu+iMu, mu, enf);
            ((TGraphErrors*)(*grENF)[iDCR])->
                  SetPointError(iRepeat*nMu+iMu, muErr, enfErr);
         }
      }
   }   
   TCanvas* canv = new TCanvas("canv", "canv", canvX, canvY);
   ((TGraphErrors*)(*grENF)[0])->Draw("apl");
   ((TGraphErrors*)(*grENF)[0])->SetTitle("ENF;mu;ENF");
   ((TGraphErrors*)(*grENF)[0])->SetMarkerStyle(20);
   
   for(Int_t iDCR=0; iDCR<nDCR; iDCR++) {
      if(iDCR==0) continue;
      ((TGraphErrors*)(*grENF)[iDCR])->Draw("same pl");
      ((TGraphErrors*)(*grENF)[iDCR])->SetMarkerStyle(20);
      ((TGraphErrors*)(*grENF)[iDCR])->SetMarkerColor(iDCR+2);
      ((TGraphErrors*)(*grENF)[iDCR])->SetLineColor(iDCR+2);
      ((TGraphErrors*)(*grENF)[iDCR])->SetMinimum(0);
   }
   leg->Draw();
}
