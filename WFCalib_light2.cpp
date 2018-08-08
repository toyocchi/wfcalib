//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.cpp"

void WFCalib_light2() {
	const Int_t N = 256;
	Double_t canvX = 1150., canvY = 650;
	
	TRandom3 rndm(0);

	Double_t integRange = 250;
	Double_t trueGain=0.04485;

	Double_t lambda = 0;
	Double_t alpha = 0;
	Double_t alphaCT = 0;
	Double_t dcr = 0;
	Double_t noise = 0;

	Int_t nEvent = 40000;
	Int_t nWidth = 3;
	Int_t nMu = 10;

	TClonesArray *caCanv = new TClonesArray("TCanvas", nWidth);
	TClonesArray *caCanvPar = new TClonesArray("TCanvas", nWidth);

	for(Int_t iWidth=0; iWidth<nWidth; iWidth++) {
		Double_t width = 10 + 40*iWidth;
		cout << endl <<  "width: " << width << endl;
		Char_t canvName[N];
		sprintf(canvName, "DiffVar (width: %.0f)", width);
		new ((TCanvas*)(*caCanv)[iWidth]) 
			TCanvas(canvName, canvName, canvX, canvY);
		((TCanvas*)(*caCanv)[iWidth])->Divide(4,3);
		sprintf(canvName, "Parameters(width: %.0f)", width);
		new ((TCanvas*)(*caCanvPar)[iWidth]) 
			TCanvas(canvName, canvName, canvX, canvY);
		((TCanvas*)(*caCanvPar)[iWidth])->Divide(2,2);

		TClonesArray *caPar = new TClonesArray("TGraphErrors", 4);
		for(Int_t iPar=0; iPar<4; iPar++) {
			new ((TGraphErrors*)(*caPar)[iPar]) TGraphErrors;
		}
		
		for(Int_t iMu=0; iMu<nMu; iMu++) {		
			Int_t mu = iMu+1;
			cout << "mu: " << mu << endl;
			TGraph *grDiffVar = new TGraph;
			Int_t pnt = 0;
			for(int iEvent=0; iEvent<nEvent; iEvent++){
				Waveform* wf = new Waveform(integRange);
				wf->SetNoises(lambda, alpha, alphaCT, dcr, noise);
				wf->SetTauTiming(width);
				Int_t nPhoton = rndm.Poisson(mu);
				wf->MakeEvent(nPhoton);
				Double_t charge = wf->GetCharge();
				if(charge<.5*trueGain) {
					continue;
				}
				pnt++;
				wf->Differentiate(2);
				Double_t diffVar = wf->GetVariance();
				grDiffVar->SetPoint(pnt, charge, diffVar);
				delete wf;
				pnt++;
			}
			
			grDiffVar->SetMinimum(0);
			grDiffVar->SetMarkerStyle(8);
			grDiffVar->SetMarkerSize(.1);
			((TCanvas*)(*caCanv)[iWidth])->cd(iWidth+1);
			grDiffVar->Draw("ap");
			
			TF1* funcQuad = new TF1("funcQuad", "pol2");
			funcQuad->SetLineColor(kRed);
			grDiffVar->Fit(funcQuad, "Q");
			
			
			Double_t parDiffVar[3];
			Double_t errDiffVar[3];
			for (int iPar = 0; iPar < 3; iPar++) {
				parDiffVar[iPar]=funcQuad->GetParameter(iPar);
				errDiffVar[iPar]=funcQuad->GetParError(iPar);
				((TGraphErrors*)(*caPar)[iPar])->
					SetPoint(iMu, mu, parDiffVar[iPar]);
				((TGraphErrors*)(*caPar)[iPar])->
					SetPointError(iMu, 0, errDiffVar[iPar]);
			}
			
			Double_t gainDiffVar = parDiffVar[1]/(1-parDiffVar[2]);
			Double_t gainDiffVarErr = 
				TMath::Sqrt(TMath::Power(errDiffVar[1]/(1-parDiffVar[2]), 2)
								+ TMath::Power(parDiffVar[1]*errDiffVar[2]
													/(1-parDiffVar[2])/(1-parDiffVar[2])
													, 2));
			
			((TGraphErrors*)(*caPar)[3])->
				SetPoint(iMu, mu, gainDiffVar);
			((TGraphErrors*)(*caPar)[3])->SetPointError(iMu, 0, gainDiffVarErr);
		}
		
		for(Int_t iPar=0; iPar<4; iPar++) {
			((TGraphErrors*)(*caPar)[iPar])->SetMarkerStyle(8);
			((TGraphErrors*)(*caPar)[iPar])->SetMarkerSize(.5);
			((TCanvas*)(*caCanvPar)[iWidth])->cd(iPar+1);
			((TGraphErrors*)(*caPar)[iPar])->Draw("ap");
		}
	}
}
