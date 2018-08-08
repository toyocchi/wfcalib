//--/----|----/----|----/----|----/----|----/----|----/----|----/----|----/----|
#include "Waveform_LED.h"

Waveform::Waveform(Double_t integRange, Bool_t led) {
   fNPnt = fSamplingRate * integRange;
   fAnalysisEnd = fAnalysisStart + integRange;
   fIntegRange = integRange;
   fLED = led;

   fTime = new Double_t[fNPnt];
   fSignalAmp = new Double_t[fNPnt];
   fNoiseAmp  = new Double_t[fNPnt];

   fFuncSingle = new TF1("funcSingle", "[0]/[1]*exp(-x/[1])",
								 fAnalysisStart, fAnalysisEnd);
   fFuncSingle->SetParameters(fTrueGain, fTauDecay);
   fFuncTiming = new TF1("funcTiming", "exp(-x/[0])", 0, 99999);
   fFuncTiming->SetParameters(0, fTauTiming);
   fFuncAPTiming = new TF1("funcAPTiming", "exp(-x/[0])", 0, 99999);
   fFuncAPTiming->SetParameters(0, fTauAP);

   memset(fSignalAmp, 0, sizeof(Double_t)*fNPnt);
   memset(fNoiseAmp,  0, sizeof(Double_t)*fNPnt);

   fPntSize = (fAnalysisEnd-fAnalysisStart)/fNPnt;
   for (int iPnt = 0; iPnt < fNPnt; iPnt++) {
      fTime[iPnt] = fAnalysisStart + fPntSize*iPnt;
   }
}


Waveform::~Waveform(){
   delete[] fTime;
   delete[] fSignalAmp;
   delete[] fNoiseAmp;
   delete fFuncSingle;
}


void Waveform::Hit(Double_t pulseTime) {
	TRandom3 rndm(0);

   fHit++;
   for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
      fSignalAmp[iPnt] += fFuncSingle->Eval(fTime[iPnt]-pulseTime);
   }

   // prompt crosstalk
   if(fLambda>0) {
      Int_t numCT = rndm.Poisson(fLambda);
      for(Int_t iCT=0; iCT<numCT; iCT++) {
         Hit(pulseTime);
      }
   }


   // delayed cross talk                                                      
   if(fAlphaCT>0) {
      if(rndm.Uniform() < fAlphaCT) {
         Double_t dCTTime = fFuncAPTiming->GetRandom();
         Double_t dCTPulseTime = pulseTime + dCTTime;

         if(pulseTime < fAnalysisEnd) {
            for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
               fSignalAmp[iPnt] += fFuncSingle->Eval(fTime[iPnt]-dCTPulseTime);
            }

            // Crosstalk of delayed crosstalk      
            if(fLambda>0) {
               Int_t numCT = rndm.Poisson(fLambda);
               for(Int_t iCT=0; iCT<numCT; iCT++) {
                  Hit(dCTPulseTime);
               }
            }
         }
      }
   }

   // Afterpulsing
   if(fAlpha>0) {
      if(rndm.Uniform() < fAlpha) {
         Double_t apTime= fFuncAPTiming->GetRandom();
         Double_t apPulseTime = pulseTime + apTime;

         if(pulseTime < fAnalysisEnd) {
            for(Int_t iPnt=0; iPnt<fNPnt; iPnt++) {
               fSignalAmp[iPnt] +=
                  (1-exp(-apTime/fTauDecay)) *
                  fFuncSingle->Eval(fTime[iPnt]-apPulseTime);
            }

            // Cross Talk of After Pulse
            Int_t numCT = rndm.Poisson(fLambda);
            for(Int_t iCT=0; iCT<numCT; iCT++) {
               Hit(pulseTime);
            }
         }
      }
   }
}


Int_t Waveform::MakeEvent(Int_t nPhoton) {
	TRandom3 rndm(0);
	
   fHit = 0;
   for (size_t iPhoton = 0; iPhoton < nPhoton; iPhoton++) {
      Double_t pulseTime;
      if(fLED) {
         pulseTime = rndm.Uniform(fTauTiming);
      } else {
         pulseTime = fFuncTiming->GetRandom();
      }
      pulseTime += fTriggerTime;

      if(pulseTime < fAnalysisEnd) {
         Hit(pulseTime);
      }
   }

   // dark
   if(fDCR>0) {
      Double_t nExp = fDCR*pow(10,-6) * fIntegRange;
      Int_t nDC = rndm.Poisson(nExp);
      for(Int_t iDC=0; iDC<nDC; iDC++) {
         Double_t pulseTime = rndm.Uniform(fAnalysisStart, fAnalysisEnd);
         Hit(pulseTime);
      }
   }
	
	SetNoise(fNoise);

   return fHit;
}


void Waveform::Differentiate(Int_t nDiff) {
   if (nDiff>0) {
      Double_t* diffSignalAmp= new Double_t[fNPnt];
      Double_t* diffNoiseAmp= new Double_t[fNPnt];

      for(int iDiff=0; iDiff<nDiff; iDiff++){
         for (int iPnt = 1; iPnt < fNPnt-1; iPnt++) {
            diffSignalAmp[iPnt] =
               (fSignalAmp[iPnt+1]-fSignalAmp[iPnt-1])/fPntSize;
            diffNoiseAmp[iPnt] =
               (fNoiseAmp[iPnt+1]-fNoiseAmp[iPnt-1])/fPntSize;
         }
         diffSignalAmp[0+iDiff] = 0;
         diffNoiseAmp[0+iDiff] = 0;
         diffSignalAmp[fNPnt-1-iDiff] = 0;
         diffNoiseAmp[fNPnt-1-iDiff]  = 0;

         memcpy(fSignalAmp,diffSignalAmp,fNPnt*sizeof(Double_t) );
         memcpy(fNoiseAmp,diffNoiseAmp,fNPnt*sizeof(Double_t) );
      }
   }
}

