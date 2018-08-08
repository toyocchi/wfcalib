//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.cpp"

void WFCalib_ov_linear() {
	const Int_t N = 256;
	Double_t canvX = 1150., canvY = 650;
	
	TRandom3 rndm(0);

	Double_t samplingRate = 16;
	Double_t integRange = 1500;

	Double_t tauTiming = 1000;

	Double_t lambda = 0;
	Double_t alpha = 0;
	Double_t alphaCT = 0;
	Double_t dcr = 0;
	Double_t noise = 0;

	Int_t nEvent = 100;
	Int_t nPhoton = 5;
	Int_t nOV = 8;
	Int_t nPar = 2;

	TCanvas *canvDiffVar = new TCanvas("canvDiffVar", "DiffVar", canvX, canvY);
	canvDiffVar->Divide(3,3);
	TCanvas *canvMonitor = new TCanvas("canvMonitor", "Monitor", canvX, canvY);
	canvMonitor->Divide(4,4);

	TClonesArray *caPar = new TClonesArray("TGraphErrors", nPar);
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
				Waveform* wf = new Waveform;
				wf->SetSamplingRate(samplingRate);
				wf->SetIntegRange(integRange);
				wf->SetTauTiming(tauTiming);
				wf->SetTrueGain(trueGain);

				wf->MakeTemplate();
				wf->SetNoises(lambda, alpha, alphaCT, dcr, noise);
				wf->MakeEvent(iPhoton);
				Double_t charge = wf->GetCharge();
				if(charge<.5*trueGain) {
					continue;
				}
				pnt++;
				if(iOV==0 && iEvent==0) {
					canvMonitor->cd(2*iPhoton+1);
					wf->MakeGraph();
				}

				wf->Differentiate(2);
				Double_t diffVar = wf->GetVariance();
				grDiffVar->SetPoint(pnt, charge, diffVar);
				if(iOV==0 && iEvent==0) {
					canvMonitor->cd(2*iPhoton+2);
					wf->MakeGraph();
				}
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
		
		TF1* funcLinear = new TF1("funcLinear", "pol1");
		funcLinear->SetLineColor(kRed);
		grDiffVar->Fit(funcLinear, "Q");
		
		Double_t parDiffVar[3];
		Double_t errDiffVar[3];
		for (int iPar = 0; iPar < nPar; iPar++) {
			parDiffVar[iPar]=funcLinear->GetParameter(iPar);
			errDiffVar[iPar]=funcLinear->GetParError(iPar);

			((TGraphErrors*)(*caPar)[iPar])->
			 SetPoint(iOV, trueGain, parDiffVar[iPar]);
			((TGraphErrors*)(*caPar)[iPar])->
			 SetPointError(iOV, 0, errDiffVar[iPar]);
		}
	}

	TCanvas *canvPar = new TCanvas("canvPar", "Parameters", canvX, canvY);
	canvPar->Divide(2,2);

	((TGraphErrors*)(*caPar)[0])->SetTitle("C0;Gain;C0");
	((TGraphErrors*)(*caPar)[1])->SetTitle("C1;Gain;C1");

	for(Int_t iPar=0; iPar<nPar; iPar++) {
		((TGraphErrors*)(*caPar)[iPar])->SetMarkerStyle(8);
		((TGraphErrors*)(*caPar)[iPar])->SetMarkerSize(.5);
		canvPar->cd(iPar+1);
		((TGraphErrors*)(*caPar)[iPar])->Draw("ap");
	}
}
