#include "phFunc.h"
Double_t ExpectedSpectrum(Double_t *x,Double_t *par);
void phrec(){

	Double_t kmax  = 50;
	Double_t alpha = 0.2;
	Double_t mu    = 2;
	Double_t lambda= 0.1;
	Double_t sigma_0=0.1;
	Double_t sigma_1=0.0;
	Double_t beta   =100;
	// std::cout<<"debug"<<std::endl;
	TF1* fexp= new TF1("fexp",ExpectedSpectrum,-0.5,20,7);
	// std::cout<<"debug"<<std::endl;
	fexp->SetParameters(kmax,alpha,mu,lambda,sigma_0,sigma_1,beta);
	// std::cout<<"debug"<<std::endl;
	// TF1* fborel= new TF1("fborel",Borel);
	// fborel->Draw();
	// fexp->Eval(0);
	// std::cout<<fexp->Eval(0)<<std::endl;
	fexp->Draw("pl");
	fexp->SetMarkerStyle(3);
}

Double_t ExpectedSpectrum(Double_t *x,Double_t *par){
	Double_t xx    = x[0];
	Double_t kmax  = par[0];
	Double_t alpha = par[1];
	Double_t mu    = par[2];
	Double_t lambda= par[3];
	Double_t sigma_0=par[4];
	Double_t sigma_1=par[5];
	Double_t beta   =par[6];
	Double_t gain   =1.;
	Double_t ped    =0.;
	// std::cout<<"debug"<<std::endl;
	TF1* fgpoi = new TF1("fgpoi",GeneralizedPoisson,-100,100,2);
	TF1* fgph  = new TF1("fgph",GaussPH,-100,100,4);
	TF1* fborel= new TF1("fborel",Borel,-100,100,2);
	TF1* fdpdph_1= new TF1("fdpdph_1",SingleAPProb,-100,100,5);
	TF1* fdpdph_i= new TF1("fdpdph_i",DiffProb,-100,100,5);
	TF1* fsigmak= new TF1("fsigmak",SigmaK,0,100,2);
	// std::cout<<"debug"<<std::endl;
	fgpoi->SetParameters(mu,lambda);
	fgph-> SetParameters(0,sigma_0,ped,gain);
	// std::cout<<"ped: "<<ped<<std::endl;
	Double_t vgpoi_ped=fgpoi-> Eval(0);
	Double_t vgph_ped=fgph-> Eval(xx);
	Double_t term_ped = vgpoi_ped*vgph_ped;
	Double_t term_ill = 0;
	// std::cout<<"xx: "<<xx<<std::endl;
	// std::cout<<"term_ped"<<term_ped<<std::endl;
	// std::cout<<"gpoi"<<vgpoi_ped<<std::endl;
	// std::cout<<"gph"<<vgph_ped<<std::endl;
	fsigmak->SetParameters(sigma_0,sigma_1);
	for(Int_t k=1;k<kmax+1;k++){
		// std::cout<<"k: "<<k<<std::endl;
		Double_t kk=(Double_t)k;
		Double_t sigma_k= fsigmak->Eval(kk);
		// std::cout<<"sigma_k: "<<sigma_k<<std::endl;
		fgpoi  -> SetParameters(mu,lambda);
		Double_t GP= fgpoi->Eval(kk);
		// 0 term
		fborel -> SetParameters(0,alpha);
		fgph   -> SetParameters(kk,sigma_k,ped,gain);
		Double_t vborel_0 = fborel->Eval(kk);
		Double_t vgrph_0  = fgph->Eval(xx);
		Double_t term_0 = vborel_0*vgrph_0;
		// 1 term
		fborel -> SetParameters(1,alpha);
		fdpdph_1 -> SetParameters(kk,beta,sigma_k,ped,gain);
		Double_t vborel_1= fborel->Eval(kk);
		Double_t vdpdph_1  = fdpdph_1->Eval(xx);
		Double_t term_1=vborel_1* vdpdph_1;
		// i > 2
		Double_t term_i=0;
		for(Int_t i=2;i<k+1;i++){
			// std::cout<<"i: "<<i<<std::endl;
			fborel -> SetParameters(i,alpha);
			fdpdph_i -> SetParameters(kk,i,beta,ped,gain);
			term_i+=fborel->Eval(kk)*fdpdph_i->Eval(xx);
		}
		term_ill+=GP*(term_0+term_1+term_i);
		// term_ill+=GP*term_0;
		// std::cout<<"term_0"<<term_0<<std::endl;
		// std::cout<<"term_1"<<term_1<<std::endl;
		// std::cout<<"term_i"<<term_i<<std::endl;
		// std::cout<<"GP    "<<GP<<std::endl;

	}
	// std::cout<<"term_ill"<<term_ill<<std::endl;
	return term_ped+ term_ill;
}