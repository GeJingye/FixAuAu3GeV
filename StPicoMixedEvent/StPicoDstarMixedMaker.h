#ifndef StPicoDstarMixedMaker_H
#define StPicoDstarMixedMaker_H

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoDstarEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */
 /* ******************************************************************************************
  * read PicoDst document about FixedTarget AuAu 3.85GeV for produciton within TOF acceptance*
  * ******************************************************************************************
  */
#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "StPicoEvent/StPicoDst.h"
#include "../StRefMultCorr/StRefMultCorr.h"
#include "../StRefMultCorr/CentralityMaker.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLorentzVector.h"
class TString;
class TFile;
class TNtuple;
class StPicoTrack;
class StPicoDstMaker;
class StPicoEvent;
class ParticleInfo;

class ParticleInfo
{
  public:
	Int_t charge;
	Float_t pt;
	Float_t eta;
	Float_t phi;
	Float_t p;
	Float_t nSigmaE;
	// Float_t nSigmaPi;
	// Float_t nSigmaK;
	// Float_t nSigmaP;
	Float_t beta;

	Float_t energy;
	Float_t p1;
	Float_t p2;
	Float_t p3;
	Bool_t isPhotonicE;//isPhotonicE
	Bool_t isPureE;
};

class StPicoDstarMixedMaker : public StMaker
{
  public:
    StPicoDstarMixedMaker(Char_t const *name, TString const inputFilesList, TString const outBaseName, StPicoDstMaker* picoDstMaker);
    virtual ~StPicoDstarMixedMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

  private:
    StPicoDstMaker* mPicoDstMaker;
	TString mInputFilesList;
	TString mOutFileBaseName;

    StPicoDstarMixedMaker() {}
    void initHists();
    Bool_t isGoodTrigger(StPicoEvent const*) const;
    Bool_t isGoodEvent(StPicoEvent const* const picoEvent) const;
    Bool_t isGoodTrack(StPicoTrack const* trk, StPicoEvent const* const picoEvent) const;

    Float_t getTofBeta(StPicoTrack const* const trk) const;
	Double_t getPhiVAngle(TLorentzVector e1, TLorentzVector e2, Int_t q1, Int_t q2) const;
	Double_t getNSigmaECorr(TVector3 mom) const;
	Double_t getNSigmaPiKPCorr(Int_t num_variable, TVector3 mom) const;
	void copyCurrentToBuffer();
    
    // -------------- USER variables -------------------------
    // add your member variables here. 
    // Remember that ntuples size can be really big, use histograms where appropriate 

  public:
    void setRunNumList(string list){
      mRunNumList = list;
    }
    void getBadruns(string inputFileName);
    void getBadTofCells(string inputFileName);
    Bool_t isBadrun(Int_t runId);
    Bool_t isGoodBTofCell(Int_t const tofid);

  private: 
    TFile* mFile;
	// std::map<Key, T> 是 C++ 标准库提供的有序关联容器。建立“键→值”的一一映射，支持对数复杂度（O(log n)）的插入、查找、删除，且按键自动排序。
	std::map<Int_t, Int_t> mrunnum;
	string mRunNumList;
	// std::vector<T> 是 C++ 标准库提供的动态顺序容器。“可变长数组”，在尾部高效 O(1) 地插入 / 删除，支持随机访问 O(1)，内存连续。
    vector<Int_t> mBadRun;
    vector<Int_t> mBadTofCell;
    // Event level
    Int_t  mRunId;
    // primaryVertex vertex
    Float_t mVx;
    Float_t mVy;
    Float_t mVr;
    Float_t mVz;

	Float_t mBfield;// 磁场强度
	Int_t mCen9;
	TF1* fphiVcut;

	// Histograms
    TH1F* h_RunNum;
	TH1F* h_RunNum__goodTrigger;
    TH1F* h_passEvtcut;
	TH1F* h_passTrkcut;
 
    // event QA
	TH1F* h_cen9;
	TH1F* h_cen9__gE;
	TH1F* h_cen9__gT;
	TH2F* h_nTofMat_fxtMult;
	TH1F* h_fxtMult;
    TH3F* h_Vx_Vy_Vz;
	TH2F* h_Vx_Vy;
	TH1F* h_Vr;
    TH1F* h_Vz;
    TH3F* h_pDca_Eta_NHitsFit;
	TH3F* h_pDca_Pt_Eta;

    // track level QA
    TH1F* h_nHitsFit;
    TH1F* h_nHitsPoss;
    TH1F* h_nHitsDEdx;
	TH3F* h_nHitsFit_Pt_Eta;
    TH3F* h_nHitsDEdx_Pt_Eta;

	TH1F* h_m2;
    TH2F* h_m2_Pc;
    TH1F* h_pDca;
    TH1F* h_ppT;
	TH1F* h_pP;
	TH2F* h_pP_ppT;
	TH1F* h_gPt;
    TH1F* h_pEta;
    TH1F* h_pPhi;
    TH2F* h_ppTc_pEta;
    TH2F* h_ppTc_pPhi;

	TH2F* h_dEdx_Pc;
	TH2F* h_nSigmaElectron_P;
	TH2F* h_nSigmaPion_P;
	TH2F* h_nSigmaKaon_P;
	TH2F* h_nSigmaProton_P;
	TH2F* h_nSigmaEcorr_P;
	TH2F* h_nSigmaPicorr_P;
    
	TH3F* h_Pt_Cen_nSigmaE__PureE;
	TH3F* h_Eta_Cen_nSigmaE__PureE;
	TH3F* h_Phi_Cen_nSigmaE__PureE;
	TH3F* h_Pt_Cen_nSigmaPion;
	TH3F* h_Eta_Cen_nSigmaPion;
	TH3F* h_Phi_Cen_nSigmaPion;

	// group 1
	TH2F* h_invBeta_P__TOFMatch;
	TH1F* h_pT__TOFMatch;
	TH1F* h_Eta__TOFMatch;
	TH1F* h_Phi__TOFMatch;
	TH2F* h_nSigmaElectron_P__1;
	TH2F* h_nSigmaProton_P__1;
	TH2F* h_nSigmaElectron_P__TOFMatch;
	TH2F* h_nSigmaElectron_P__PIDcut_1;
	// group 2
	TH2F* h_nSigmaElectron_P__2;
	TH2F* h_nSigmaPion_P__2;
	TH2F* h_nSigmaKaon_P__2;
	TH2F* h_nSigmaProton_P__2;
	TH2F* h_nSigmaElectron_P__PIDcut_2;
	// group 3
	TH2F* h_nSigmaElectron_P__3;
	TH2F* h_nSigmaPion_P__3;
	TH2F* h_nSigmaKaon_P__3;
	TH2F* h_nSigmaProton_P__3;
	TH2F* h_nSigmaElectron_P__PIDcut_3;
	TH2F* h_nSigmaElectron_Eta__EIDcut_3;
	TH2F* h_nSigmaElectron_P__EIDcut_total;
	// phiV cut
	TH2F* h_eNumber_Cen;
	TH1F* h_pT__positrons;
	TH1F* h_eta__positrons;
	TH1F* h_phi__positrons;
	TH1F* h_pT__positrons_w_PhiV_Cut;
	TH1F* h_eta__positrons_w_PhiV_Cut;
	TH1F* h_phi__positrons_w_PhiV_Cut;
	TH1F* h_pT__electrons;
	TH1F* h_eta__electrons;
	TH1F* h_phi__electrons;
	TH1F* h_pT__electrons_w_PhiV_Cut;
	TH1F* h_eta__electrons_w_PhiV_Cut;
	TH1F* h_phi__electrons_w_PhiV_Cut;

	TH1F* h_Rapidity__unlikeSame;
    TH2F* h_Mee_PhiV__unlikeSame;
    TH1F* h_Mee__unlikeSame;
    TH1F* h_Mee__unlikeSame__w_PhiV_Cut;

    TH3F* h_Mee_Pt_Cen__unlikeSame;
    TH3F* h_Mee_Pt_Cen__likemm;
    TH3F* h_Mee_Pt_Cen__likepp;

    TH3F* h_Mee_Pt_Cen__unlikeMixed;
    TH3F* h_Mee_Pt_Cen__likemmMixed;
    TH3F* h_Mee_Pt_Cen__likeppMixed;

    const Float_t M_electron=0.000511;// GeV
	const Float_t M_proton = 0.938;

    ClassDef(StPicoDstarMixedMaker, 1)
};

inline void StPicoDstarMixedMaker::getBadruns(string inputFileName)
{
    ifstream fin(inputFileName.c_str());
    if(!fin)
    {
      cout <<"no Bad runs list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId )
    {
      mBadRun.push_back(runId);
    }
    cout << "get Bad runs list: [OK]" << endl;
}

inline void StPicoDstarMixedMaker::getBadTofCells(string inputFileName)
{
    ifstream fin(inputFileName.c_str());
    if(!fin)
    {
      cout << "no Bad BTOF Cells list" << endl;
      return;
    }
    cout << "  " << inputFileName << flush;

    Int_t runId = 0 ;
    while( fin >> runId )
    {
      mBadTofCell.push_back(runId);
    }
    cout << "get Bad BTOF Cells list: [OK]" << endl;
}

inline Bool_t StPicoDstarMixedMaker::isBadrun(Int_t runId)
{
    vector<Int_t>::iterator iter = std::find(mBadRun.begin(), mBadRun.end(), runId);
    return ( iter != mBadRun.end() ) ;
}

inline Bool_t StPicoDstarMixedMaker::isGoodBTofCell(Int_t const tofid)
{
  if(std::count(mBadTofCell.begin(),mBadTofCell.end(),tofid)>0) return kFALSE;
  else return kTRUE;
}

#endif
