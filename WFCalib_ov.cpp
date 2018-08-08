//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.cpp"

void WFCalib_ov() {
	const Int_t N = 256;
	Double_t canvX = 1150., canvY = 650;
	
	TRandom3 rndm(0);

	Double_t lambda = 0;
	Double_t alpha = 0;
	Double_t alphaCT = 0;
	Double_t dcr = 0;
	Double_t noise = 0;

	Int_t nEvent = 100;
	Int_t nPhoton = 50;
	Int_t nOV = 8;
	Int_t nPar = 4;

	TCanvas *canvDiffVar = new TCanvas("canvDiffVar", "DiffVar", canvX, canvY);
	canvDiffVar->Divide(3,3);

	TClonesArray *caPar = new TClonesArray("TGraphErrors", 4);
	for(Int_t iPar=0; iPar<nPar; iPar++) {
		new ((TGraphErrors*)(*caPar)[iPar]) TGraphErrors;
	}

	for(Int_t iOV=0; iOV<nOV; iOV++) {
		Double_t trueGain = .02 + .01*iOV;
		cout << "trueGain: " << trueGain << endl;

		TGraph *grDiffVar = new TGraph;
		Int_t pnt = 0;		
		for(Int_t iPhoton=0; iPhoton<nPhoton; iPhoton++) {
			for(int iEvent=0; iEvent<nEvent; iEvent++){
				Waveform* wf = new Waveform(trueGain);
				wf->SetNoises(lambda, alpha, alphaCT, dcr, noise);
				wf->MakeEvent(iPhoton);
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
		}

		grDiffVar->SetMinimum(0);
		grDiffVar->SetMarkerStyle(8);
		grDiffVar->SetMarkerSize(.3);
		Char_t nameGrDiffVar[N];
		sprintf(nameGrDiffVar, "DiffVar (Gain = %.2f);Charge;DiffVar", trueGain);
		grDiffVar->SetTitle(nameGrDiffVar);
		canvDiffVar->cd(iOV+1);
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
			 SetPoint(iOV, trueGain, parDiffVar[iPar]);
			((TGraphErrors*)(*caPar)[iPar])->
			 SetPointError(iOV, 0, errDiffVar[iPar]);
		}
		Double_t gainDiffVar = parDiffVar[1]/(1-parDiffVar[2]);
		Double_t gainDiffVarErr = 
			TMath::Sqrt(TMath::Power(errDiffVar[1]/(1-parDiffVar[2]), 2)
							+ TMath::Power(parDiffVar[1]*errDiffVar[2]
												/(1-parDiffVar[2])/(1-parDiffVar[2])
												, 2));
		((TGraphErrors*)(*caPar)[3])->SetPoint(iOV, trueGain, gainDiffVar);
		((TGraphErrors*)(*caPar)[3])->SetPointError(iOV, 0, gainDiffVarErr);
	}

	TCanvas *canvPar = new TCanvas("canvPar", "Parameters", canvX, canvY);
	canvPar->Divide(2,2);

	for(Int_t iPar=0; iPar<4; iPar++) {
		((TGraphErrors*)(*caPar)[iPar])->SetMarkerStyle(8);
		((TGraphErrors*)(*caPar)[iPar])->SetMarkerSize(.5);
		canvPar->cd(iPar+1);
		((TGraphErrors*)(*caPar)[iPar])->Draw("ap");
	}
}
