/* **************************************************
 *                                                  *
 *  Authors: Yuanjing Ji                            *
 *           Guannan Xie <guannanxie@lbl.gov>       *
 *           Mustafa Mustafa <mmustafa@lbl.gov>     *
 *                                                  *
 * **************************************************
 */
 /* ******************************************************************************************
  * read PicoDst document about FixedTarget AuAu 3.85GeV for produciton within TOF acceptance*
  * ******************************************************************************************
 */
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoETofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoDstarMixedMaker.h"
#include "StAnaCuts.h"
#include "StMemStat.h"//?
//#include "calmean.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"

#ifndef DEBUG
#define DEBUG 1
#endif

const Int_t kMagBins = 1;// 磁场的取向（+/-）
const Int_t kCenBins = anaCuts::nCenBins;//中心度分成的bin数
const Int_t kVzBins = 1;// Vz分成的bin数
Int_t magBufferIndex, cenBufferIndex, vzBufferIndex;// 索引到每个bin
const Int_t kMaxEventsInBuffer = 500;// 每个buffer中事件数最大值
const Int_t kMaxElectrons = 50;// 每个event中正、负电子数的上界

Int_t current_nPositron;// 当前事件中通过筛选的正电子（e+）数
Int_t current_nElectron;// 当前事件中通过筛选的电子（e-）数
TLorentzVector current_positron[kMaxElectrons];// 当前事件中通过筛选的正电子的四动量
TLorentzVector current_electron[kMaxElectrons];// 当前事件中通过筛选的电子的四动量

Int_t nEventsInBuffer[kMagBins][kCenBins][kVzBins];// 每个buffer每个bin中事件数
Bool_t bufferFullFlag[kMagBins][kCenBins][kVzBins];// 每个buffer每个bin是否满员
Int_t buffer_nEPlus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer];// 每个buffer每个bin每个事件中正电子数
Int_t buffer_nEMinus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer];// 每个buffer每个bin每个事件中电子数
TLorentzVector buffer_ePlus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer][kMaxElectrons];// 每个buffer每个bin每个事件每个正电子的四动量
TLorentzVector buffer_eMinus[kMagBins][kCenBins][kVzBins][kMaxEventsInBuffer][kMaxElectrons];// 每个buffer每个bin每个事件每个电子的四动量

ClassImp(StPicoDstarMixedMaker)

// 成员初始化函数的类外实现
StPicoDstarMixedMaker::StPicoDstarMixedMaker(Char_t const *name, TString const inputFilesList, TString const outFileBaseName, StPicoDstMaker* picoDstMaker):
    StMaker(name),// char*实际上是指向存储字符串的字符数组第一个元素的地址
	mPicoDstMaker(picoDstMaker),
    mInputFilesList(inputFilesList),
	mOutFileBaseName(outFileBaseName)
{}
StPicoDstarMixedMaker::~StPicoDstarMixedMaker() {}

Int_t StPicoDstarMixedMaker::Init()
{
  mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
  mFile = new TFile(mOutFileBaseName+".root", "RECREATE");

  // initialize histograms and trees
  initHists();
  return kStOK;
}

void StPicoDstarMixedMaker::initHists()
{
	ifstream readnum;
	readnum.open(mRunNumList);// mRunNumList由readPico.cxx通过setRunNumList函数赋值，即runnumber.list
	Int_t runum = 0;
	Int_t totalNum = 0;
	//if (DEBUG) cout << "start initial run number..." << endl;
	while (readnum >> runum) {
		mrunnum[runum] = totalNum;
		totalNum++;
		//if (DEBUG) cout << "run number : " << runum << " id :" << mrunnum[runum] << endl;
	}
	readnum.close();


	// 把一段连续的内存 每个字节 都设为0,例如：Int_t buf[1024]; memset(buf, 0, sizeof(buf));
	memset(nEventsInBuffer, 0, sizeof(nEventsInBuffer));
	memset(bufferFullFlag, 0, sizeof(bufferFullFlag));
	memset(buffer_nEPlus, 0, sizeof(buffer_nEPlus));
	memset(buffer_nEMinus, 0, sizeof(buffer_nEMinus));

	// Event histograms
	fphiVcut = new TF1("fphiVcut", "0.84326*exp(-49.4819*x)-0.996609*x+0.19801", 0.0, 1.0);
	h_RunNum = new TH1F("h_RunNum", "h_RunNum", totalNum, -0.5, totalNum-0.5);
	h_RunNum__goodTrigger = new TH1F("h_RunNum__goodTrigger", "h_RunNum__goodTrigger", totalNum, -0.5, totalNum - 0.5);
	h_cen9 = new TH1F("h_cen9", "h_cen9", 9, 0., 9.);
	h_cen9__gE = new TH1F("h_cen9__gE", "h_cen9__gE", 9, 0., 9.);
	h_cen9__gT = new TH1F("h_cen9__gT", "h_cen9__gT", 9, 0., 9.);
	h_cen9->GetXaxis()->SetBinLabel(9, "0~5%");
	h_cen9->GetXaxis()->SetBinLabel(8, "5~10%");
	h_cen9->GetXaxis()->SetBinLabel(7, "10~20%");
	h_cen9->GetXaxis()->SetBinLabel(6, "20~30%");
	h_cen9->GetXaxis()->SetBinLabel(5, "30~40%");
	h_cen9->GetXaxis()->SetBinLabel(4, "40~50%");
	h_cen9->GetXaxis()->SetBinLabel(3, "50~60%");
	h_cen9->GetXaxis()->SetBinLabel(2, "60~70%");
	h_cen9->GetXaxis()->SetBinLabel(1, "70~80%");
	h_fxtMult = new TH1F("h_fxtMult", "h_fxtMult", 250, 0, 250);// 参考多重数
	h_nTofMat_fxtMult = new TH2F("h_nTofMat_fxtMult", "fxtMult vs nTofmatch;nTofMatch;fxtMult", 250, 0, 250, 250, 0, 250);// 与TOF匹配的径迹数vs多重数关系
	h_passEvtcut = new TH1F("h_passEvtcut", "pass event cut", 7, -0.5, 6.5);
	h_passEvtcut->GetXaxis()->SetBinLabel(1, "All");
	h_passEvtcut->GetXaxis()->SetBinLabel(2, "Good Runs");
	h_passEvtcut->GetXaxis()->SetBinLabel(3, "passVz");
	h_passEvtcut->GetXaxis()->SetBinLabel(4, "passVr");
	h_passEvtcut->GetXaxis()->SetBinLabel(5, "passVerror");
	h_passEvtcut->GetXaxis()->SetBinLabel(6, "notPileUp");
	h_passEvtcut->GetXaxis()->SetBinLabel(7, "0-80%");

	h_passTrkcut = new TH1F("h_passTrkcut", "tracks in different conditions", 3, -0.5, 2.5);
	h_passTrkcut->GetXaxis()->SetBinLabel(1, "All");
	h_passTrkcut->GetXaxis()->SetBinLabel(2, "Primary Tracks");
	h_passTrkcut->GetXaxis()->SetBinLabel(3, "Good Tracks");

	//顶点位置
	h_Vx_Vy_Vz = new TH3F("h_Vx_Vy_Vz", "Vz vs Vy vs Vx;Vx(cm);Vy(cm);Vz(cm)", 250, -5, 5, 250, -5, 5, 100, 150, 250);
	h_Vx_Vy = new TH2F("h_Vx_Vy", "Vy vs Vx;Vx(cm);Vy(cm)", 1400, -7, 7, 1400, -7, 7);
	h_Vr = new TH1F("h_Vr", "Vr;Vr(cm);Counts", 400, 0, 4);
	h_Vz = new TH1F("h_Vz", "Vz;Vz(cm);Counts", 400, 198, 202);

	// 径迹信息&主径迹信息
	h_nHitsFit = new TH1F("h_nHitsFit", "nHitsFit;nHitsFit", 160, -80., 80.);
	h_nHitsPoss = new TH1F("h_nHitsPoss", "nHitsPoss;nHitsPoss", 160, -80., 80.);
	h_nHitsDEdx = new TH1F("h_nHitsDEdx", "nHitsDedx;nHitsDedx", 160, -80., 80.);
	h_nHitsFit_Pt_Eta = new TH3F("h_nHitsFit_Pt_Eta", "nHitsFit vs p_{T} vs #eta;p_{T} (GeV/c);#eta;nHitsFit", 100, 0., 10., 35, -3., 0.5, 80, 0., 80.);
	h_nHitsDEdx_Pt_Eta = new TH3F("h_nHitsDEdx_Pt_Eta", "nHitsDedx vs p_{T} vs #eta;p_{T} (GeV/c);#eta;nHitsDedx", 100, 0., 10., 35, -3., 0.5, 80, 0., 80.);

	h_pDca = new TH1F("h_pDca", "pDca;DCA;counts", 50, 0., 5.);// p代表primary，pDCA指主径迹与重建顶点的最小距离
	h_ppT = new TH1F("h_ppT", "primary p_{T};p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_pP = new TH1F("h_pP", "primary p;p (GeV/c);counts", 1000, 0., 10.);
	h_pP_ppT = new TH2F("h_pP_ppT", "primary p vs primary p_{T};p (GeV/c);p_{T} (GeV/c)", 1000, 0., 10., 1000, 0., 10.);
	h_gPt = new TH1F("h_gPt", "global p_{T} of track;global p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_pEta = new TH1F("h_pEta", "primary #eta;#eta;counts", 300, -3., 0.5);
	h_pPhi = new TH1F("h_pPhi", "primary #phi;#phi;counts", 80, -4.0, 4.0);
	h_ppTc_pEta = new TH2F("h_ppTc_pEta", "p_{T}*q vs #eta;p_{T}*q (GeV/c);#eta", 400, -10, 10, 300, -3., 0.5);
	h_ppTc_pPhi = new TH2F("h_ppTc_pPhi", "p_{T}*q vs #phi;p_{T}*q (GeV/c);#phi", 400, -10, 10, 800, -4, 4);
	h_pDca_Eta_NHitsFit = new TH3F("h_pDca_Eta_NHitsFit", "NHitsFit vs #eta vs DCA;DCA;#eta;NHitsFit", 50, 0., 5., 35, -3., 0.5, 90, 0, 90);
	h_pDca_Pt_Eta = new TH3F("h_pDca_Pt_Eta", "#eta vs p_{T} vs DCA;DCA;p_{T};#eta", 50, 0., 5., 1000, 0., 10., 35, -3., 0.5);
	//GoodTrack径迹信息
	h_nSigmaElectron_P = new TH2F("h_nSigmaElectron_P", "n#sigma_{e} vs p;p (GeV/c);n#sigma_{e}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaPion_P = new TH2F("h_nSigmaPion_P", "n#sigma_{pi} vs p;p (GeV/c);n#sigma_{pi}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaKaon_P = new TH2F("h_nSigmaKaon_P", "n#sigma_{K} vs p;p (GeV/c);n#sigma_{K}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaProton_P = new TH2F("h_nSigmaProton_P", "n#sigma_{P} vs p;p (GeV/c);n#sigma_{P}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaEcorr_P = new TH2F("h_nSigmaEcorr_P", "corrected n#sigma_{e} vs p;p (GeV/c);corrected n#sigma_{e}", 500, 0, 5, 4000, -20, 20);
	h_nSigmaPicorr_P = new TH2F("h_nSigmaPicorr_P", "corrected n#sigma_{pi} vs p;p (GeV/c);corrected n#sigma_{pi}", 500, 0, 5, 4000, -20, 20);
	h_dEdx_Pc = new TH2F("h_dEdx_Pc", "dE/dx vs p*q;p*q(GeV/c);#frac{dE}{dx} (GeV cm^{2}/g)", 1000, -5, 5, 400, 0, 25);
	h_m2    = new TH1F("h_m2", "m^{2};m^{2};counts", 8000, -0.5, 15.5);
	h_m2_Pc = new TH2F("h_m2_Pc", "m^{2} vs p*q;p*q (GeV/c);m^{2} (GeV/c^{2})^{2}", 1000, -5, 5, 1600, -0.5, 15.5);
	//nSigmaE
	h_Pt_Cen_nSigmaE__PureE = new TH3F("h_Pt_Cen_nSigmaE__PureE", "n#sigma_{e} vs p_{T} vs Cen;p_{T} (GeV/c);Cen;n#sigma_{e}", 50,0,5, 9,0,9, 200, -10, 10);
	h_Eta_Cen_nSigmaE__PureE = new TH3F("h_Eta_Cen_nSigmaE__PureE", "n#sigma_{e} vs #eta vs Cen;#eta;Cen;n#sigma_{e}", 30,-2.5,0.5, 9,0,9, 200, -10, 10);
	h_Phi_Cen_nSigmaE__PureE = new TH3F("h_Phi_Cen_nSigmaE__PureE", "n#sigma_{e} vs #phi vs Cen;#phi;Cen;n#sigma_{e}", 64,-3.2,3.2, 9,0,9, 200, -10, 10);
	//h_Pt_Cen_nSigmaEcorr__PureE = new TH3F("h_Pt_Cen_nSigmaEcorr__PureE", "n#sigma_{e} vs p_{T} vs Cen;p_{T} (GeV/c);Cen;n#sigma_{e}", 500, 0, 5, 9, 0, 9, 1000, -10, 10);
	//h_Eta_Cen_nSigmaEcorr__PureE = new TH3F("h_Eta_Cen_nSigmaEcorr__PureE", "n#sigma_{e} vs #eta vs Cen;#eta;Cen;n#sigma_{e}", 300, -2.5, 0.5, 9, 0, 9, 1000, -10, 10);
	//h_Phi_Cen_nSigmaEcorr__PureE = new TH3F("h_Phi_Cen_nSigmaEcorr__PureE", "n#sigma_{e} vs #phi vs Cen;#phi;Cen;n#sigma_{e}", 640, -3.2, 3.2, 9, 0, 9, 1000, -10, 10);
	//nSigmaPi
	h_Pt_Cen_nSigmaPion = new TH3F("h_Pt_Cen_nSigmaPion", "n#sigma_{pi} vs p_{T} vs Cen;p_{T} (GeV/c);Cen;n#sigma_{pi}", 500, 0, 5, 9, 0, 9, 1000, -10, 10);
	h_Eta_Cen_nSigmaPion = new TH3F("h_Eta_Cen_nSigmaPion", "n#sigma_{pi} vs #eta vs Cen;#eta;Cen;n#sigma_{pi}", 300, -2.5, 0.5, 9, 0, 9, 1000, -10, 10);
	h_Phi_Cen_nSigmaPion = new TH3F("h_Phi_Cen_nSigmaPion", "n#sigma_{pi} vs #phi vs Cen;#phi;Cen;n#sigma_{pi}", 640, -3.2, 3.2, 9, 0, 9, 1000, -10, 10);

	// group 1
	h_invBeta_P__TOFMatch = new TH2F("h_invBeta_P__TOFMatch", "1/#beta vs pc;pc (GeV/c);1/#beta", 500, 0., 5., 5000, 0, 5);
	h_pT__TOFMatch = new TH1F("h_pT__TOFMatch", "p_{T} of TOF matched tracks;p_{T} (GeV/c);counts", 1000, 0., 10.);
	h_Eta__TOFMatch = new TH1F("h_Eta__TOFMatch", "#eta of TOF matched tracks;#eta;counts", 300, -3., 0.5);
	h_Phi__TOFMatch = new TH1F("h_Phi__TOFMatch", "#phi of TOF matched tracks;#phi;counts", 80, -4, 4);
	h_nSigmaElectron_P__1 = new TH2F("h_nSigmaElectron_P__1", "n#sigma_{e} vs p (p_{T}>0.2,-1.4<#eta<0);p (GeV/c);n#sigma_{e}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaProton_P__1 = new TH2F("h_nSigmaProton_P__1", "n#sigma_{p} vs p (p_{T}>0.2,-1.4<#eta<0);p (GeV/c);n#sigma_{p}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaElectron_P__TOFMatch = new TH2F("h_nSigmaElectron_P__TOFMatch", "n#sigma_{e} vs p (p_{T}>0.2,-1.4<#eta<0);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_P__PIDcut_1 = new TH2F("h_nSigmaElectron_P__PIDcut_1", "n#sigma_{e} vs p (p_{T}>0.2,-1.4<#eta<0);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	// group 2
	h_nSigmaElectron_P__2 = new TH2F("h_nSigmaElectron_P__2", "n#sigma_{e} vs p (p_{T}>0.2,-2.5<#eta<-1.4);p (GeV/c);n#sigma_{e}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaPion_P__2 = new TH2F("h_nSigmaPion_P__2", "n#sigma_{pi} vs p (p_{T}>0.2,-2.5<#eta<-1.4);p (GeV/c);n#sigma_{pi}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaKaon_P__2 = new TH2F("h_nSigmaKaon_P__2", "n#sigma_{K} vs p (p_{T}>0.2,-2.5<#eta<-1.4);p (GeV/c);n#sigma_{K}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaProton_P__2 = new TH2F("h_nSigmaProton_P__2", "n#sigma_{p} vs p (p_{T}>0.2,-2.5<#eta<-1.4);p (GeV/c);n#sigma_{p}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaElectron_P__PIDcut_2 = new TH2F("h_nSigmaElectron_P__PIDcut_2", "n#sigma_{e} vs p (p_{T}>0.2,-2.5<#eta<-1.4);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	// group 3
	h_nSigmaElectron_P__3 = new TH2F("h_nSigmaElectron_P__3", "n#sigma_{e} vs p (p_{T}<0.2,#eta<0);p (GeV/c);n#sigma_{e}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaPion_P__3 = new TH2F("h_nSigmaPion_P__3", "n#sigma_{pi} vs p (p_{T}<0.2,#eta<0);p (GeV/c);n#sigma_{pi}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaKaon_P__3 = new TH2F("h_nSigmaKaon_P__3", "n#sigma_{K} vs p (p_{T}<0.2,#eta<);p (GeV/c);n#sigma_{K}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaProton_P__3 = new TH2F("h_nSigmaProton_P__3", "n#sigma_{p} vs p (p_{T}<0.2,#eta<0);p (GeV/c);n#sigma_{p}", 500, 0, 5, 1000, -10, 10);
	h_nSigmaElectron_P__PIDcut_3 = new TH2F("h_nSigmaElectron_P__PIDcut_3", "n#sigma_{e} vs p (p_{T}>0.2,-2.5<#eta<-1.4);p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	h_nSigmaElectron_Eta__EIDcut_3 = new TH2F("h_nSigmaElectron_Eta__EIDcut_3", "n#sigma_{e} vs p (p_{T}>0.2,-2.5<#eta<-1.4);#eta;n#sigma_{e}", 300, -2.5, 0.5, 2000, -10, 10);
	// total
	h_nSigmaElectron_P__EIDcut_total = new TH2F("h_nSigmaElectron_P__EIDcut_total", "n#sigma_{e} vs p;p (GeV/c);n#sigma_{e}", 500, 0, 5, 2000, -10, 10);
	// 经TPC和TOF判选后的电子信息
	h_eNumber_Cen = new TH2F("h_eNumber_Cen", ";Num.;Cen", 50, 0., 50., 9, 0., 9.);
	h_pT__electrons = new TH1F("h_pT__electrons", "p_{T} of electrons;p_{T} (GeV/c)", 50, 0, 1.);
	h_eta__electrons = new TH1F("h_eta__electrons", "#eta of electrons;#eta", 35, -3., 0.5);
	h_phi__electrons = new TH1F("h_phi__electrons", "#phi of electrons;#phi", 80, -4, 4);
	h_pT__electrons_w_PhiV_Cut = new TH1F("h_pT__electrons_w_PhiV_Cut", "p_{T} of electrons with #phi_{V} cut;p_{T} (GeV/c)", 50, 0., 1.);//"w"表示with, "wo"表示without
	h_eta__electrons_w_PhiV_Cut = new TH1F("h_eta__electrons_w_PhiV_Cut", "#eta of electrons with #phi_{V} cut;#eta", 35, -3., 0.5);
	h_phi__electrons_w_PhiV_Cut = new TH1F("h_phi__electrons_w_PhiV_Cut", "RFF #phi of electrons with #phi_{V} cut;#phi", 80, -4, 4);
	h_pT__positrons = new TH1F("h_pT__positrons", "p_{T} of positrons;p_{T} (GeV/c)", 50, 0, 1.);
	h_eta__positrons = new TH1F("h_eta__positrons", "#eta of positrons;#eta", 35, -3., 0.5);
	h_phi__positrons = new TH1F("h_phi__positrons", "#phi of positrons;#phi", 80, -4, 4);
	h_pT__positrons_w_PhiV_Cut = new TH1F("h_pT__positrons_w_PhiV_Cut", "p_{T} of positrons with #phi_{V} cut;p_{T} (GeV/c)", 50, 0., 1.);//"w"表示with, "wo"表示without
	h_eta__positrons_w_PhiV_Cut = new TH1F("h_eta__positrons_w_PhiV_Cut", "#eta of positrons with #phi_{V} cut;#eta", 35, -3., 0.5);
	h_phi__positrons_w_PhiV_Cut = new TH1F("h_phi__positrons_w_PhiV_Cut", "RFF #phi of positrons with #phi_{V} cut;#phi", 80, -4, 4);

	h_Rapidity__unlikeSame = new TH1F("h_Rapidity__unlikeSame", "y distribution of e^{+}e^{-};y;counts", 300, -2.5, 0.5);
	h_Mee_PhiV__unlikeSame = new TH2F("h_Mee_PhiV__unlikeSame", "Mee vs #phi_{V};Mee(GeV/c^{2});#phi_{V}", 800, 0, 4, 100, 0, 1);
	h_Mee__unlikeSame = new TH1F("h_Mee__unlikeSame", "Mee without #phi_{V} cut;Mee(GeV/c^{2})", 800, 0, 4);
	h_Mee__unlikeSame__w_PhiV_Cut = new TH1F("h_Mee__unlikeSame__w_PhiV_Cut", "Mee with #phi_{V} cut;Mee(GeV/c^{2})", 800, 0, 4);

	h_Mee_Pt_Cen__unlikeSame = new TH3F("h_Mee_Pt_Cen__unlikeSame", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 9, 0., 9.);
	h_Mee_Pt_Cen__likemm = new TH3F("h_Mee_Pt_Cen__likemm", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 9, 0., 9.);
	h_Mee_Pt_Cen__likepp = new TH3F("h_Mee_Pt_Cen__likepp", "Mee vs p_{T} vs Cen;Mee (GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 9, 0., 9.);

	h_Mee_Pt_Cen__unlikeMixed = new TH3F("h_Mee_Pt_Cen__unlikeMixed", "Mee vs p_{T} vs Cen;Mee(GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 9, 0., 9.);
	h_Mee_Pt_Cen__likemmMixed = new TH3F("h_Mee_Pt_Cen__likemmMixed", "Mee vs p_{T} vs Cen;Mee(GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 9, 0., 9.);
	h_Mee_Pt_Cen__likeppMixed = new TH3F("h_Mee_Pt_Cen__likeppMixed", "Mee vs p_{T} vs Cen;Mee(GeV/c^{2});p_{T} (GeV/c);Cen", 800, 0, 4, 100, 0, 5, 9, 0., 9.);

}// 

Int_t StPicoDstarMixedMaker::Make()
{
  ParticleInfo particleinfo;
  vector<ParticleInfo> electroninfo;
  vector<ParticleInfo> positroninfo;

  // StMemStat mem;
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoDstarMixedMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  if (!picoDst)
  {
    LOG_WARN << "StPicoDstarMixedMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  // -------------- USER ANALYSIS -------------------------


  StPicoEvent const *picoEvent = picoDst->event();

  mRunId = picoEvent->runId();
  h_RunNum->Fill(mrunnum[mRunId]);// 填充不同run号的事例数；
  if (!isGoodTrigger(picoEvent)) return kStOK;
  h_RunNum__goodTrigger->Fill(mrunnum[mRunId]);// 填充不同run号的事例数；
  // 清空正负电子信息缓存
  electroninfo.clear();
  positroninfo.clear();


  TVector3 pVtx = picoEvent->primaryVertex();// 获得TPC重建的该事例顶点的三维坐标
  mVx = pVtx.x();
  mVy = pVtx.y();
  mVz = pVtx.z();
  mVr = sqrt(pow(mVx+0.5,2)+pow(mVy+2,2));
  
  // event and track level QA
  h_passEvtcut->Fill(0);// 原始事例数+1

  if (!isBadrun(mRunId))// bad run list 
  {
	h_passEvtcut->Fill(1);// 通过bad run cut的事例数+1
	h_Vx_Vy_Vz->Fill(mVx, mVy, mVz);
	h_Vz->Fill(mVz);
	h_Vx_Vy->Fill(mVx, mVy);
	h_Vr->Fill(mVr);


	// 获取固定靶多重数以及计算中心度
	Double_t fxtMult = picoEvent->fxtMult();
	Int_t mCen9;
	if (fxtMult > 183) {
		mCen9 = 8;//0~5%
	}
	else if (fxtMult > 155) {
		mCen9 = 7;//5~10%
	}
	else if (fxtMult > 111) {
		mCen9 = 6;//10~20 %
	}
	else if (fxtMult > 77) {
		mCen9 = 5;//20~30 %
	}
	else if (fxtMult > 51) {
		mCen9 = 4;//30~40 %
	}
	else if (fxtMult > 33) {
		mCen9 = 3;//40~50 %
	}
	else if (fxtMult > 19) {
		mCen9 = 2;//50~60 %
	}
	else if (fxtMult > 10) {
		mCen9 = 1;//60~70 %
	}
	else if (fxtMult > 5) {
		mCen9 = 0;//70~80%
	}
	else {
		mCen9 = -1;
		return kStOk;
	}
	h_cen9->Fill(mCen9);// 填充中心度
	h_fxtMult->Fill(fxtMult);
	h_nTofMat_fxtMult->Fill(picoEvent->nBTOFMatch(), fxtMult);

	// 不同条件cut后的事例数统计
    Bool_t vzcut = mVz < anaCuts::Vz_up && mVz > anaCuts::Vz_low;
	Bool_t vrcut = mVr < anaCuts::Vr;
    //Bool_t verrcut = fabs(mVx) > anaCuts::Verror && fabs(mVy) > anaCuts::Verror && fabs(mVz) > anaCuts::Verror;// |Vx|,|Vy,|Vz|>1.0e-5 cm, why? too small that better than resolution?
	Bool_t notPileUp = kTRUE;// !mRefMultCorrUtil->isPileUpEvent(mRefmult6, picoEvent->nBTOFMatch(), mVz, mTotnMIP);
    Bool_t cen0280cut = mCen9 > -1;

    if (vzcut) h_passEvtcut->Fill(2);
    if (vzcut && vrcut) h_passEvtcut->Fill(3);
    if (vzcut && vrcut && notPileUp ) h_passEvtcut->Fill(5);
    if (vzcut && vrcut && notPileUp && cen0280cut)
    {
      // ******************以下分析均基于good event, notPileup, 0~80%中心度***********************
	  h_passEvtcut->Fill(6);
	  h_cen9__gE->Fill(mCen9);
	  mBfield = picoEvent->bField();// 获取磁场


      magBufferIndex = 0;
      cenBufferIndex = mCen9;
	  vzBufferIndex = 0;

	  current_nPositron = 0;
      current_nElectron = 0;

	  // ************************开始径迹分析*****************************
      Int_t nTracks = picoDst->numberOfTracks();
	  for (Int_t itrack = 0; itrack < nTracks; itrack++)
	  {
		  h_passTrkcut->Fill(0);
		  StPicoTrack* trk = picoDst->track(itrack);
		  // ******************以下分析均基于主径迹**************************
		  Bool_t isPrimaryTrack = trk->isPrimary();
		  if (!isPrimaryTrack) continue;

		  h_passTrkcut->Fill(1);
		  TVector3 mom = trk->pMom();// 提取当前径迹在“主顶点处”的三维动量矢量？
		  Float_t mgDCAs = trk->gDCA(picoEvent->primaryVertex()).Mag();// 返回该径迹与主顶点的最短距离

		  // 获得nSigmaElectron等图像
		  Double_t nSigmaE = trk->nSigmaElectron();
		  Double_t nSigmaPi = trk->nSigmaPion();
		  Double_t nSigmaK = trk->nSigmaKaon();
		  Double_t nSigmaP = trk->nSigmaProton();

		  Double_t beta = getTofBeta(trk);
		  Double_t m2 = pow(mom.Mag()*1.0 / beta, 2)*(1 - beta * beta);
		  Bool_t isTOFMatch = !std::isnan(beta) && beta > 0;
		  Bool_t isBemcMatch = trk->isBemcMatchedTrack();

		  // 填充直方图
		  if (fabs(mgDCAs + 999.0) > 1e-2) h_pDca->Fill(mgDCAs);// STAR约定：当 DCA 无法计算（例如径迹缺 hit、拟合失败）时，把 mgDCAs 设为 - 999 cm 作为无效标志。
		  h_pP_ppT->Fill(mom.Mag(), mom.Perp());
		  h_ppT->Fill(mom.Perp());
		  h_pP->Fill(mom.Mag());
		  h_gPt->Fill(trk->gPt());
		  h_pEta->Fill(mom.Eta());
		  h_pPhi->Fill(mom.Phi());
		  h_ppTc_pEta->Fill(mom.Perp()*trk->charge(), mom.Eta());// mom.Perp()指横动量;charge(){return (mNHitsFit > 0) ? 1 : -1;}只能返回正负1
		  h_ppTc_pPhi->Fill(mom.Perp()*trk->charge(), mom.Phi());

		  h_nHitsFit->Fill(trk->nHitsFit()*trk->charge());
		  h_nHitsPoss->Fill(trk->nHitsPoss()*trk->charge());
		  h_nHitsDEdx->Fill(trk->nHitsDedx()*trk->charge());
		  h_nHitsFit_Pt_Eta->Fill(mom.Perp(), mom.Eta(), trk->nHitsFit());
		  h_nHitsDEdx_Pt_Eta->Fill(mom.Perp(), mom.Eta(), trk->nHitsDedx());
		  h_pDca_Eta_NHitsFit->Fill(mgDCAs, mom.Eta(), trk->nHitsFit());
		  h_pDca_Pt_Eta->Fill(mgDCAs, mom.Perp(), mom.Eta());
		  // ******************以下分析均基于满足goodTrackCuts的主径迹**************************
		  if (!isGoodTrack(trk, picoEvent)) continue;
		  h_passTrkcut->Fill(2);
		  h_cen9__gT->Fill(mCen9);
		  h_nSigmaElectron_P->Fill(mom.Mag(), nSigmaE);
		  h_nSigmaPion_P->Fill(mom.Mag(), nSigmaPi);
		  h_nSigmaKaon_P->Fill(mom.Mag(), nSigmaK);
		  h_nSigmaProton_P->Fill(mom.Mag(), nSigmaP);

		  h_dEdx_Pc->Fill(mom.Mag()*trk->charge(), trk->dEdx());// mom.Mag()指动量模
		  h_m2->Fill(m2);
		  h_m2_Pc->Fill(mom.Mag()*trk->charge(), m2);

		  // Recalibrate nSigmaElectron nSigmaPion
		  Double_t nSigmaE_corrfactor = getNSigmaECorr(mom);
		  Double_t nSigmaE_corr = nSigmaE - nSigmaE_corrfactor;
		  Double_t nSigmaPi_corrfactor = getNSigmaPiKPCorr(trk->charge() > 0 ? 1 : -1, mom);
		  Double_t nSigmaPi_corr = nSigmaPi - nSigmaPi_corrfactor;
		  h_nSigmaEcorr_P->Fill(mom.Mag(), nSigmaE_corr);
		  h_nSigmaPicorr_P->Fill(mom.Mag(), nSigmaPi_corr);

		  if (nSigmaPi > -3.5 && nSigmaPi < 3.5)
		  {
			  h_Pt_Cen_nSigmaPion->Fill(mom.Perp(), nSigmaPi, mCen9);
			  h_Eta_Cen_nSigmaPion->Fill(mom.Eta(), nSigmaPi, mCen9);
			  h_Phi_Cen_nSigmaPion->Fill(mom.Phi(), nSigmaPi, mCen9);
		  }

		  // Electrons IDentification
		  Bool_t isTPCElectron__1 = kFALSE;
		  Bool_t isTPCProton__1 = kFALSE;
		  Bool_t isBemcElectron__1 = kFALSE;
		  Bool_t isElectronGroup__1 = kFALSE;

		  Bool_t isTPCElectron__2 = kFALSE;
		  Bool_t isTPCPion__2 = kFALSE;
		  Bool_t isTPCKaon__2 = kFALSE;
		  Bool_t isTPCProton__2 = kFALSE;
		  Bool_t isTPCX__2 = kFALSE;
		  Bool_t isTPCDeuteron__2 = kFALSE;
		  Bool_t isElectronGroup__2 = kFALSE;

		  Bool_t isTPCElectron__3 = kFALSE;
		  Bool_t isTPCPion__3 = kFALSE;
		  Bool_t isTPCKaon__3 = kFALSE;
		  Bool_t isTPCProton__3 = kFALSE;
		  Bool_t isElectronGroup__3 = kFALSE;

		  // group 1
		  isTPCProton__1 = trk->nSigmaProton() > -2 && trk->nSigmaProton() < -1*mom.Mag()+5;
		  isTPCElectron__1 = ((mom.Mag() <= 1.0) ? (nSigmaE < 3. && nSigmaE > (3.5 * mom.Mag() - 4.5)) : (nSigmaE<3. && nSigmaE > -1.)) && mom.Mag() < 1.6;
		  if (isTOFMatch)
		  {
			  StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(trk->bTofPidTraitsIndex());
			  Int_t tofid = tofPid->btofCellId();

			  if (isGoodBTofCell(tofid))// && fabs(tofPid->btofYLocal()) < 1.8)//详见note_2013_JieZhao
			  {
				  h_pT__TOFMatch->Fill(mom.Perp());
				  h_Eta__TOFMatch->Fill(mom.Eta());
				  h_Phi__TOFMatch->Fill(mom.Phi());
				  if (1)//mom.Perp() > 0.2 && mom.Eta() > -1.5)
				  {
					  h_invBeta_P__TOFMatch->Fill(mom.Mag(), 1. / beta);
					  h_nSigmaElectron_P__1->Fill(mom.Mag(), nSigmaE);
					  h_nSigmaProton_P__1->Fill(mom.Mag(), nSigmaP);

					  if (fabs(1.0/beta - 1.0) < 0.025)
				      {
						h_nSigmaElectron_P__TOFMatch->Fill(mom.Mag(), nSigmaE);
						if (!isTPCProton__1)
						{
							h_nSigmaElectron_P__PIDcut_1->Fill(mom.Mag(), nSigmaE);
							isElectronGroup__1 = isTPCElectron__1;
						}
					  }
				  }
			  }
		  }
		  // EMC cut(有问题，还需要修改)(dilepton一般用不到BEMC)
		  //if (isBemcMatch && 0)
		  //{
			  //float E0, bemcE, bemcZ, bemcPhi, bemcTowerPhi, bemcTowerEta, bemcSmdN_eta, bemcSmdN_phi;
			  //int bemcId = trk->bemcPidTraitsIndex();
			  //if (bemcId >= 0)
			  //{
				 // StPicoBEmcPidTraits * bemcTrait = picoDst->bemcPidTraits(bemcId);
				 // E0 = bemcTrait->bemcE0();
				 // bemcE = bemcTrait->bemcE();
				 // bemcZ = bemcTrait->bemcZDist();
				 // bemcPhi = bemcTrait->bemcPhiDist();
				 // bemcTowerPhi = bemcTrait->btowPhiDist();
				 // bemcTowerEta = bemcTrait->btowEtaDist();
				 // bemcSmdN_eta = bemcTrait->bemcSmdNEta();
				 // bemcSmdN_phi = bemcTrait->bemcSmdNPhi();
				 // if (mom.Perp() >= 1.0)
				 // {
					//  isBemcElectron__1 = E0/mom.Mag() < 1.5 && E0/mom.Mag() > 0.5;
				 // }
			  //}
		  //}

		  // group 2
		  isTPCPion__2 = nSigmaPi>-4.0 && nSigmaPi<3.0;
		  isTPCKaon__2 = nSigmaK > -2.0 && nSigmaK < 2.0;;
		  isTPCProton__2 = nSigmaP>-3 && nSigmaP<-1.2*mom.Mag()+6.0;
		  isTPCX__2 = kFALSE;
		  isTPCDeuteron__2 = kFALSE;
		  isTPCElectron__2 = nSigmaE>-2.0 && nSigmaE<2.0 && mom.Mag()<1.5;
		  if (mom.Perp() > 0.2 && mom.Eta() > -1.4)
		  {
			  h_nSigmaElectron_P__2->Fill(mom.Mag(), nSigmaE);
			  h_nSigmaPion_P__2->Fill(mom.Mag(), nSigmaPi);
			  h_nSigmaKaon_P__2->Fill(mom.Mag(), nSigmaK);
			  h_nSigmaProton_P__2->Fill(mom.Mag(), nSigmaP);
			  if (!isTPCPion__2 && !isTPCKaon__2 && !isTPCProton__2)
			  {
				  h_nSigmaElectron_P__PIDcut_2->Fill(mom.Mag(), nSigmaE);
				  isElectronGroup__2 = isTPCElectron__2;
			  }
		  }

		  // group 3
		  isTPCPion__3 = nSigmaPi > -4.0 && nSigmaPi < 3.0;
		  isTPCKaon__3 = kFALSE;
		  isTPCProton__3 = kFALSE;
		  isTPCElectron__3 = nSigmaE > -4.0 && nSigmaE < 2.5 && mom.Mag() < 0.7;
		  if (mom.Perp() < 0.2 && mom.Eta() < 0) 
		  {
			  h_nSigmaElectron_P__3->Fill(mom.Mag(), nSigmaE);
			  h_nSigmaPion_P__3->Fill(mom.Mag(), nSigmaPi);
			  h_nSigmaKaon_P__3->Fill(mom.Mag(), nSigmaK);
			  h_nSigmaProton_P__3->Fill(mom.Mag(), nSigmaP);
			  if (!isTPCPion__3 && !isTPCKaon__3 && !isTPCProton__3)
			  {
				  h_nSigmaElectron_P__PIDcut_3->Fill(mom.Mag(), nSigmaE);
				  if (isTPCElectron__3)
				  {
					  isElectronGroup__3 = kTRUE;
					  h_nSigmaElectron_Eta__EIDcut_3->Fill(mom.Eta(),nSigmaE);
				  }
			  }
		  }

		  //填充EID得到的电子信息
		  if (isElectronGroup__1)
		  {
			h_nSigmaElectron_P__EIDcut_total->Fill(mom.Mag(), nSigmaE);
			if (trk->charge() < 0) // electron
			{
				particleinfo.charge = trk->charge();
				particleinfo.pt = mom.Perp();
				particleinfo.eta = mom.Eta();
				particleinfo.phi = mom.Phi();
				particleinfo.p = mom.Mag();
				particleinfo.nSigmaE = nSigmaE;
				particleinfo.beta = beta;
				particleinfo.energy = sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0));
				particleinfo.p1 = mom.X();
				particleinfo.p2 = mom.Y();
				particleinfo.p3 = mom.Z();
				particleinfo.isPhotonicE = kFALSE;
				particleinfo.isPureE = kFALSE;
				electroninfo.push_back(particleinfo);

				//current_electron[current_nElectron].SetPx(mom.X());
				//current_electron[current_nElectron].SetPy(mom.Y());
				//current_electron[current_nElectron].SetPz(mom.Z());
				//current_electron[current_nElectron].SetE(sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0)));
				//current_nElectron++;
			}
			if (trk->charge() > 0) // positron
			{
				particleinfo.charge = trk->charge();
				particleinfo.pt = mom.Perp();
				particleinfo.eta = mom.Eta();
				particleinfo.phi = mom.Phi();
				particleinfo.p = mom.Mag();
				particleinfo.nSigmaE = nSigmaE;
				particleinfo.beta = beta;
				particleinfo.energy = sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0));
				particleinfo.p1 = mom.X();
				particleinfo.p2 = mom.Y();
				particleinfo.p3 = mom.Z();
				particleinfo.isPhotonicE = kFALSE;
				particleinfo.isPureE = kFALSE;
				positroninfo.push_back(particleinfo);

				//current_positron[current_nPositron].SetPx(mom.X());
				//current_positron[current_nPositron].SetPy(mom.Y());
				//current_positron[current_nPositron].SetPz(mom.Z());
				//current_positron[current_nPositron].SetE(sqrt(pow(M_electron, 2.0) + pow(mom.Mag(), 2.0)));
				//current_nPositron++;
			}
          }// 填充单径迹的正/负电子信息
      }// 填充单事例所有径迹的正负电子信息

      Int_t x=0;
      Int_t y=0;
      Int_t num_electron = electroninfo.size();
      Int_t num_positron = positroninfo.size();
	  h_eNumber_Cen->Fill(num_electron + num_positron, mCen9);
      TLorentzVector eepair(0,0,0,0);
      TLorentzVector particle1_4V(0,0,0,0);
      TLorentzVector particle2_4V(0,0,0,0);
	  // \phi_v cut && pure electron sample selection
	  for (x = 0; x < num_positron; x++)
	  {
		  particle1_4V.SetPx(positroninfo[x].p1);
		  particle1_4V.SetPy(positroninfo[x].p2);
		  particle1_4V.SetPz(positroninfo[x].p3);
		  particle1_4V.SetE(positroninfo[x].energy);
		  for (y = 0; y < num_electron; y++)
		  {
			  particle2_4V.SetPx(electroninfo[y].p1);
			  particle2_4V.SetPy(electroninfo[y].p2);
			  particle2_4V.SetPz(electroninfo[y].p3);
			  particle2_4V.SetE(electroninfo[y].energy);

			  eepair = particle1_4V + particle2_4V;
			  if (eepair.M() < 0.015)
			  {
				  positroninfo[x].isPureE = kTRUE;
				  electroninfo[y].isPureE = kTRUE;
			  }
			  Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);// 注意参数1、-1的选取要求
			  Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			  h_Mee_PhiV__unlikeSame->Fill(eepair.M(), angleV);
			  if (eepair.M() < anaCuts::PhiVCutMRange && angleV < angleVcut)
			  {
				  positroninfo[x].isPhotonicE = kTRUE;
				  electroninfo[y].isPhotonicE = kTRUE;
			  }
		  }
	  }
	  for (x = 0; x < num_positron; x++)//填充current_positron，只填充isPhotonicE=False的正电子
	  {
		  if (positroninfo[x].isPureE)
		  {
			  h_Pt_Cen_nSigmaE__PureE->Fill(positroninfo[x].pt, mCen9, positroninfo[x].nSigmaE);
			  h_Eta_Cen_nSigmaE__PureE->Fill(positroninfo[x].eta, mCen9, positroninfo[x].nSigmaE);
			  h_Phi_Cen_nSigmaE__PureE->Fill(positroninfo[x].phi, mCen9, positroninfo[x].nSigmaE);
		  }
		  if (!positroninfo[x].isPhotonicE)
		  {
			  current_positron[current_nPositron].SetPx(positroninfo[x].p1);
			  current_positron[current_nPositron].SetPy(positroninfo[x].p2);
			  current_positron[current_nPositron].SetPz(positroninfo[x].p3);
			  current_positron[current_nPositron].SetE(sqrt(pow(M_electron, 2.0) + pow(positroninfo[x].p, 2.0)));
			  current_nPositron++;
		  }
	  }
	  for (x = 0; x < num_electron; x++)//填充current_positron，只填充isPhotonicE=False的电子
	  {
		  if (electroninfo[x].isPureE)
		  {
			  h_Pt_Cen_nSigmaE__PureE->Fill(positroninfo[x].pt, mCen9, positroninfo[x].nSigmaE);
			  h_Eta_Cen_nSigmaE__PureE->Fill(positroninfo[x].eta, mCen9, positroninfo[x].nSigmaE);
			  h_Phi_Cen_nSigmaE__PureE->Fill(positroninfo[x].phi, mCen9, positroninfo[x].nSigmaE);
		  }
		  if (!electroninfo[x].isPhotonicE)
		  {
			  current_electron[current_nElectron].SetPx(electroninfo[x].p1);
			  current_electron[current_nElectron].SetPy(electroninfo[x].p2);
			  current_electron[current_nElectron].SetPz(electroninfo[x].p3);
			  current_electron[current_nElectron].SetE(sqrt(pow(M_electron, 2.0) + pow(electroninfo[x].p, 2.0)));
			  current_nElectron++;
		  }
	  }

	  // 通过随机组合++, --, +-电子对重建信号
	  // +-
	  for (x = 0; x < num_positron; x++)
	  {
		  particle1_4V.SetPx(positroninfo[x].p1);
		  particle1_4V.SetPy(positroninfo[x].p2);
		  particle1_4V.SetPz(positroninfo[x].p3);
		  particle1_4V.SetE(positroninfo[x].energy);
		  for (y = 0; y < num_electron; y++)
		  {
			  particle2_4V.SetPx(electroninfo[y].p1);
			  particle2_4V.SetPy(electroninfo[y].p2);
			  particle2_4V.SetPz(electroninfo[y].p3);
			  particle2_4V.SetE(electroninfo[y].energy);
			  eepair = particle1_4V + particle2_4V;
			  Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);
			  Double_t angleVcut = fphiVcut->Eval(eepair.M());
			  //if (!positroninfo[x].isPhotonicE && !electroninfo[y].isPhotonicE)
			  if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			  {
				  h_Rapidity__unlikeSame->Fill(eepair.Rapidity());
				  h_Mee_Pt_Cen__unlikeSame->Fill(eepair.M(), eepair.Perp(), mCen9);
			  }
		  }
	  }// end: +-
	  // --
      for(x=0;x<num_electron;x++)
      {
		particle1_4V.SetPx(electroninfo[x].p1);
		particle1_4V.SetPy(electroninfo[x].p2);
		particle1_4V.SetPz(electroninfo[x].p3);
		particle1_4V.SetE(electroninfo[x].energy);
		for (y=x+1; y<num_electron; y++)//从x+1开始，避免自组合和重复组合
		{
			particle2_4V.SetPx(electroninfo[y].p1);
			particle2_4V.SetPy(electroninfo[y].p2);
			particle2_4V.SetPz(electroninfo[y].p3);
			particle2_4V.SetE(electroninfo[y].energy);
			eepair = particle1_4V + particle2_4V;
			Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);
			Double_t angleVcut = fphiVcut->Eval(eepair.M());
			//if (!electroninfo[x].isPhotonicE && !electroninfo[y].isPhotonicE)
			if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			{
				h_Mee_Pt_Cen__likemm->Fill(eepair.M(), eepair.Perp(), mCen9);
			}
		}
      }// end: for(x=0;x<num_electron;x++)
	  // ++
      for(x=0;x<num_positron;x++)
      {
        particle1_4V.SetPx(positroninfo[x].p1);
        particle1_4V.SetPy(positroninfo[x].p2);
        particle1_4V.SetPz(positroninfo[x].p3);
        particle1_4V.SetE(positroninfo[x].energy);
        for(y=x+1;y<num_positron;y++)
        {
			particle2_4V.SetPx(positroninfo[y].p1);
			particle2_4V.SetPy(positroninfo[y].p2);
			particle2_4V.SetPz(positroninfo[y].p3);
			particle2_4V.SetE(positroninfo[y].energy);
			eepair = particle1_4V + particle2_4V;
			Double_t angleV = getPhiVAngle(particle1_4V, particle2_4V, 1, -1);
			Double_t angleVcut = fphiVcut->Eval(eepair.M());
			//if (!positroninfo[x].isPhotonicE && !positroninfo[y].isPhotonicE)
			if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			{
				h_Mee_Pt_Cen__likepp->Fill(eepair.M(), eepair.Perp(), mCen9);
			}
        }
      }// end: for(x=0;x<num_positron;x++)

	  // 读取经过\phi_v cut前、后的正、负电子的横动量，赝快度，方位角
      for(x=0;x<num_positron;x++)
      {        
        h_pT__positrons->Fill(positroninfo[x].pt);
        h_eta__positrons->Fill(positroninfo[x].eta);
        h_phi__positrons->Fill(positroninfo[x].phi);
        if(!positroninfo[x].isPhotonicE)
        {
		  h_pT__positrons_w_PhiV_Cut->Fill(positroninfo[x].pt);
          h_eta__positrons_w_PhiV_Cut->Fill(positroninfo[x].eta);
          h_phi__positrons_w_PhiV_Cut->Fill(positroninfo[x].phi);
        }
      }// end: for(x=0;x<num_positron;x++)
      for(y=0;y<num_electron;y++)
      {
        h_pT__electrons->Fill(electroninfo[y].pt);
        h_eta__electrons->Fill(electroninfo[y].eta);
        h_phi__electrons->Fill(electroninfo[y].phi);
        if(!electroninfo[y].isPhotonicE)
        {        
		  h_pT__electrons_w_PhiV_Cut->Fill(electroninfo[y].pt);
          h_eta__electrons_w_PhiV_Cut->Fill(electroninfo[y].eta);
          h_phi__electrons_w_PhiV_Cut->Fill(electroninfo[y].phi);
        }
      }//

      // Event Mixing
      for(Int_t iBufferEvent=0;iBufferEvent<nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex];iBufferEvent++)
      {
        for(x=0;x<current_nPositron;x++)// unlike-sign +-
        {
			for (y = 0; y < buffer_nEMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent]; y++)
			{
				eepair = current_positron[x] + buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
				//Double_t angleV = getPhiVAngle(current_positron[x], buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
				//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
				////if(fabs(eepair.Rapidity()) <= 1.5)
				//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
				//{
					h_Mee_Pt_Cen__unlikeMixed->Fill(eepair.M(), eepair.Perp(), mCen9);
				//}
			}
        }        
        for(x=0;x<current_nElectron;x++)// unlike-sign -+
        {
          for(y=0;y<buffer_nEPlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent];y++)
          {
            eepair = current_electron[x] + buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
			//Double_t angleV = getPhiVAngle(current_electron[x], buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
			//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			////if(fabs(eepair.Rapidity()) <= 1.5)
			//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			//{
				h_Mee_Pt_Cen__unlikeMixed->Fill(eepair.M(), eepair.Perp(), mCen9);
			//}
          }
        }          
        for(x=0;x<current_nPositron;x++)// like-sign ++
        {
          for(y=0;y<buffer_nEPlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent];y++)
          {
            eepair = current_positron[x] + buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
			//Double_t angleV = getPhiVAngle(current_positron[x], buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
			//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			////if(fabs(eepair.Rapidity()) <= 1.5)
			//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			//{
				h_Mee_Pt_Cen__likeppMixed->Fill(eepair.M(), eepair.Perp(), mCen9);
			//}
          }
        }
        for(x=0;x<current_nElectron;x++)// like-sign --
        {
          for(y=0;y<buffer_nEMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent];y++)
          {
            eepair = current_electron[x] + buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y];
			//Double_t angleV = getPhiVAngle(current_electron[x], buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][iBufferEvent][y], 1, -1);// 注意参数1、-1的选取要求
			//Double_t angleVcut = fphiVcut->Eval(eepair.M()); // 根据fphiVcut关于pair-M的函数取值
			////if(fabs(eepair.Rapidity()) <= 1.5)
			//if (eepair.M() > anaCuts::PhiVCutMRange || angleV > angleVcut)
			//{
				h_Mee_Pt_Cen__likemmMixed->Fill(eepair.M(), eepair.Perp(), mCen9);
			//}
          }
        }
      }// End Mixed Event
      copyCurrentToBuffer();
    }// Good Event
  }// Good Run

  return kStOK;
}// end Make
Int_t StPicoDstarMixedMaker::Finish()
{
	mFile->cd();// 将当前 ROOT 的输出目录切换到 mFile 所指向的 ROOT 文件

	// 写入直方图
	h_RunNum->Write();
	h_RunNum__goodTrigger->Write();
	h_passEvtcut->Write();
	h_passTrkcut->Write();
	// event QA
	h_cen9->Write();
	h_cen9__gE->Write();
	h_cen9__gT->Write();
	h_Vx_Vy_Vz->Write();
	h_Vz->Write();
	h_Vr->Write();
	h_Vx_Vy->Write();
	h_nTofMat_fxtMult->Write();
	h_fxtMult->Write();
	h_pDca_Eta_NHitsFit->Write();
	h_pDca_Pt_Eta->Write();
	// track QA
	h_nHitsFit->Write();
	h_nHitsPoss->Write();
	h_nHitsDEdx->Write();
	h_nHitsFit_Pt_Eta->Write();
	h_nHitsDEdx_Pt_Eta->Write();
	// TPC TOF track

	h_dEdx_Pc->Write();
	h_m2->Write();
	h_m2_Pc->Write();
	h_pDca->Write();
	h_ppT->Write();
	h_pP->Write();
	h_pP_ppT->Write();
	h_gPt->Write();
	h_pEta->Write();
	h_pPhi->Write();
	h_ppTc_pEta->Write();
	h_ppTc_pPhi->Write();

	h_nSigmaElectron_P->Write();
	h_nSigmaPion_P->Write();
	h_nSigmaKaon_P->Write();
	h_nSigmaProton_P->Write();
	//h_nSigmaEcorr_P->Write();
	//h_nSigmaPicorr_P->Write();

	h_Pt_Cen_nSigmaE__PureE->Write();
	h_Eta_Cen_nSigmaE__PureE->Write();
	h_Phi_Cen_nSigmaE__PureE->Write();
	//h_Pt_Cen_nSigmaPion->Write();
	//h_Eta_Cen_nSigmaPion->Write();
	//h_Phi_Cen_nSigmaPion->Write();
	// group 1
	h_invBeta_P__TOFMatch->Write();
	h_pT__TOFMatch->Write();
	h_Eta__TOFMatch->Write();
	h_Phi__TOFMatch->Write();
	h_nSigmaElectron_P__1->Write();
	h_nSigmaProton_P__1->Write();
	h_nSigmaElectron_P__TOFMatch->Write();
	h_nSigmaElectron_P__PIDcut_1->Write();
	// group 2
	h_nSigmaElectron_P__2->Write();
	h_nSigmaPion_P__2->Write();
	h_nSigmaKaon_P__2->Write();
	h_nSigmaProton_P__2->Write();
	h_nSigmaElectron_P__PIDcut_2->Write();
	// group 3
	h_nSigmaElectron_P__3->Write();
	h_nSigmaPion_P__3->Write();
	h_nSigmaKaon_P__3->Write();
	h_nSigmaProton_P__3->Write();
	h_nSigmaElectron_P__PIDcut_3->Write();
	h_nSigmaElectron_Eta__EIDcut_3->Write();

	//sum
	h_nSigmaElectron_P__EIDcut_total->Write();
	h_eNumber_Cen->Write();
	h_pT__positrons->Write();
	h_eta__positrons->Write();
	h_phi__positrons->Write();
	h_pT__positrons_w_PhiV_Cut->Write();
	h_eta__positrons_w_PhiV_Cut->Write();
	h_phi__positrons_w_PhiV_Cut->Write();
	h_pT__electrons->Write();
	h_eta__electrons->Write();
	h_phi__electrons->Write();
	h_pT__electrons_w_PhiV_Cut->Write();
	h_eta__electrons_w_PhiV_Cut->Write();
	h_phi__electrons_w_PhiV_Cut->Write();

	h_Rapidity__unlikeSame->Write();
	h_Mee_PhiV__unlikeSame->Write();
	h_Mee__unlikeSame->Write();
	h_Mee__unlikeSame__w_PhiV_Cut->Write();

	h_Mee_Pt_Cen__unlikeSame->Write();
	h_Mee_Pt_Cen__likemm->Write();
	h_Mee_Pt_Cen__likepp->Write();

	h_Mee_Pt_Cen__unlikeMixed->Write();
	h_Mee_Pt_Cen__likemmMixed->Write();
	h_Mee_Pt_Cen__likeppMixed->Write();

	mFile->Close();
	return kStOK;
}

Bool_t StPicoDstarMixedMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
  for(auto trg : anaCuts::trigNumber)
  {
    if (picoEvent->isTrigger(trg)) return kTRUE;
  }
  return kFALSE;
}

Bool_t StPicoDstarMixedMaker::isGoodEvent(StPicoEvent const* const picoEvent) const
{
  TVector3 pVtx = picoEvent->primaryVertex();
  return pVtx.z() < anaCuts::Vz_up &&
		 pVtx.z() > anaCuts::Vz_low &&
		 //fabs(pVtx.x()) > anaCuts::Verror &&
		 //fabs(pVtx.y()) > anaCuts::Verror &&
		 //fabs(pVtx.z()) > anaCuts::Verror &&
		 sqrt(pow(pVtx.x() + 0.5, 2) + pow(pVtx.y() + 2, 2)) < anaCuts::Vr;
}

Bool_t StPicoDstarMixedMaker::isGoodTrack(StPicoTrack const* trk, StPicoEvent const* const picoEvent) const
{
	TVector3 mom = trk->pMom();
	return 
	trk->gPt() > anaCuts::GPt &&
    trk->gMom().Eta()> anaCuts::Eta &&
	fabs(trk->nHitsFit()) >= anaCuts::NHitsFit &&
	fabs(trk->nHitsDedx()) >= anaCuts::NHitsDedx;// &&
	trk->gDCA(picoEvent->primaryVertex()).Mag() <= anaCuts::Dca &&
	fabs(trk->nHitsFit()*1.0 / trk->nHitsMax()) >= anaCuts::NHitsFitRatio;
}

Float_t StPicoDstarMixedMaker::getTofBeta(StPicoTrack const* const trk) const
{
  Int_t index2tof = trk->bTofPidTraitsIndex();
  Float_t beta = std::numeric_limits<Float_t>::quiet_NaN();
  if (index2tof >= 0)
  {
    StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);
    if (tofPid)
    {
      beta = tofPid->btofBeta();
      if (beta < 1e-4)
      {
        TVector3 const vtx3 = mPicoDstMaker->picoDst()->event()->primaryVertex();
        StThreeVectorF vtx(vtx3.x(),vtx3.y(),vtx3.z());
        TVector3 const btofHitPos3 = tofPid->btofHitPos();
        StThreeVectorF btofHitPos(btofHitPos3.x(),btofHitPos3.y(),btofHitPos3.z());
        StPicoPhysicalHelix helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());
        Float_t L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
        Float_t tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<Float_t>::quiet_NaN();
      }
    }
  } 
  return beta;
}

Double_t StPicoDstarMixedMaker::getPhiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2) const
{
	Double_t pt1 = e1.Pt();
	Double_t eta1 = e1.Eta();
	Double_t phi1 = e1.Phi();

	Double_t pt2 = e2.Pt();
	Double_t eta2 = e2.Eta();
	Double_t phi2 = e2.Phi();

	TVector3 e1Mom,e2Mom;
	if(q1>0&&q2<0){
		e2Mom.SetPtEtaPhi(pt1,eta1,phi1);// e+
		e1Mom.SetPtEtaPhi(pt2,eta2,phi2);// e-
	}else if(q1<0&&q2>0){
		e2Mom.SetPtEtaPhi(pt2,eta2,phi2);// e+
		e1Mom.SetPtEtaPhi(pt1,eta1,phi1);// e-
	}else if(q1==q2&&TMath::Abs(q1)==1){
		Double_t ran = gRandom->Uniform(-1,1);
		if(ran>0){
			e2Mom.SetPtEtaPhi(pt1,eta1,phi1);
			e1Mom.SetPtEtaPhi(pt2,eta2,phi2);
		}
		else{
			e2Mom.SetPtEtaPhi(pt2,eta2,phi2);
			e1Mom.SetPtEtaPhi(pt1,eta1,phi1);
		}
	}else return -1;
	Double_t mN = 0.;
	if(mBfield<0.) mN = -1.;
	if(mBfield>0.) mN = 1.;

	TVector3 pu=e1Mom+e2Mom;
	TVector3 pv=e1Mom.Cross(e2Mom);
	TVector3 pw=pu.Cross(pv);
	TVector3 pnz(0.,0.,mN);
	TVector3 pwc=pu.Cross(pnz);

	Double_t angleV = pw.Angle(pwc);

	return angleV;
}


Double_t StPicoDstarMixedMaker::getNSigmaECorr(TVector3 mom) const // from zihan
{
	return 0.;
}
Double_t StPicoDstarMixedMaker::getNSigmaPiKPCorr(Int_t num_variable, TVector3 mom) const //from Wendi
{
	return 0.;
}

void StPicoDstarMixedMaker::copyCurrentToBuffer()
{
	//判断是否已满
	if (nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex] >= kMaxEventsInBuffer)
		bufferFullFlag[magBufferIndex][cenBufferIndex][vzBufferIndex] = kTRUE;

	Int_t eventIndex = -1;
	//如果已满随机替换，未满则按顺序添加event
	if (bufferFullFlag[magBufferIndex][cenBufferIndex][vzBufferIndex])
		eventIndex = gRandom->Integer(kMaxEventsInBuffer);
	else
	{
		eventIndex = nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex];
		nEventsInBuffer[magBufferIndex][cenBufferIndex][vzBufferIndex]++;
	}

	//填充该事例中e+
	buffer_nEPlus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex] = current_nPositron;
	for (Int_t i = 0; i < current_nPositron; i++)
		buffer_ePlus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex][i] = current_positron[i];

	//填充该事例中e-
	buffer_nEMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex] = current_nElectron;
	for (Int_t i = 0; i < current_nElectron; i++)
		buffer_eMinus[magBufferIndex][cenBufferIndex][vzBufferIndex][eventIndex][i] = current_electron[i];		
}