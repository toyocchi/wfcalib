//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
class Waveform{
private:
	Double_t fSamplingRate = 1.6; // GHz

protected:
   Int_t    fNPnt;
	Double_t fPntSize;
   Double_t fAnalysisStart = 0;
   Double_t fAnalysisEnd;
	Double_t fIntegRange;
	Double_t fTriggerTime = fAnalysisStart + 10;

	Double_t fTauTiming = 45;
	Double_t fTauDecay = 30;
	Double_t fTauAP = 100;
	Bool_t fLED;

	Double_t fTrueGain = .04485;
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

	void SetNoise(Double_t noise);

public:
	Waveform(Double_t integRange=150, Bool_t led=kTRUE);
   ~Waveform();

   Int_t MakeEvent(Int_t nPhoton);
	void Hit(Double_t pulseTime);
	void Differentiate(Int_t nDiff);


	void SetNoises(Double_t lambda, Double_t alpha, Double_t alphaCT,
						Double_t dcr, Double_t noise);
	void SetTauTiming(Double_t tauTiming);
   Double_t GetCharge();
   Double_t GetVariance();
};


inline void Waveform::SetNoise(Double_t noise) {
	TRandom3 rndm(0);
	for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      fNoiseAmp[iPnt] = rndm.Gaus(0, noise);
   }
}

inline void Waveform::SetTauTiming(Double_t tauTiming) {
	fTauTiming = tauTiming;
}

inline void Waveform::SetNoises(Double_t lambda, Double_t alpha,
										  Double_t alphaCT, Double_t dcr, 
										  Double_t noise) {
	fLambda = lambda;
	fAlpha = alpha;
	fAlphaCT = alphaCT;
	fDCR = dcr;
	fNoise = noise;
}

inline Double_t Waveform::GetCharge() {
   Double_t charge = 0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      charge += (fSignalAmp[iPnt]+fNoiseAmp[iPnt])*fPntSize;
   }
   return charge;
}

inline Double_t Waveform::GetVariance() {
   Double_t variance=0;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      variance += TMath::Power((fSignalAmp[iPnt]+fNoiseAmp[iPnt]),2) * fPntSize;
   }
   return variance;
}
