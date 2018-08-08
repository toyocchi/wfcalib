//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#ifndef WAVEFORM_LED_H
#define WAVEFORM_LED_H

class Waveform{
protected:
	Double_t fSamplingRate = 1.6; // GHz
	Double_t fIntegRange = 150;
   Double_t fAnalysisStart = 0;
   Double_t fAnalysisEnd;
   Int_t    fNPnt;
	Double_t fPntSize;
	Double_t fTriggerTime = fAnalysisStart + 10;

	Double_t fTrueGain = .04485;
	Double_t fTauTiming = 45;
	Double_t fTauDecay = 30;
	Double_t fTauAP = 100;
	Bool_t fLED = kTRUE;

	Double_t fLambda = .067;
	Double_t fAlpha = .32;
	Double_t fAlphaCT = 0;
	Double_t fDCR = 500; // kHz
	Double_t fNoise = .00529;

   Double_t *fSignalAmp;
   Double_t *fNoiseAmp;
   Double_t *fTime;

	Int_t fHit;

	TF1 *fFuncSingle;
	TF1 *fFuncTiming;
	TF1 *fFuncAPTiming;
	TGraph *fGraph;

	void SetNoise(Double_t noise);

public:
	Waveform();
   ~Waveform();

	void MakeTemplate();
   Int_t MakeEvent(Int_t nPhoton);
	void Hit(Double_t pulseTime);
	void Differentiate(Int_t nDiff);
	void MakeGraph();
	void DeleteGraph();

	void SetSamplingRate(Double_t samplingRate);
	void SetIntegRange(Double_t integRange);
	void SetTrueGain(Double_t trueGain);
	void SetTauTiming(Double_t tauTiming);
	void SetNoises(Double_t lambda, Double_t alpha, Double_t alphaCT,
						Double_t dcr, Double_t noise);
   Double_t GetCharge();
   Double_t GetVariance();
};

#endif // WAVEFORM_LED_H
