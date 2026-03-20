//---------------------从AuAuFixed.root中提取直方图，并进行设置更改和元素添加----------------------
#include "someFunction.h"
//1_20260201_Fixtest_forEID_29605
void DrawNSigmaEPiKP(TString inFilename = "roots/11_20260311_FixAu_group1_Mixed_1_500_50_cutPhiVByPair.root", Int_t number = 11)//6_20260227_Fix_3group_onlyMb
{
    // 从root文件中导入待拟合的直方图
	TFile *inFile = new TFile(inFilename); if (!inFile)
    {
	    cout << "The input file is not found! Exiting..." << endl;
	    return;
    }
    TH1F* h_passEvtcut = (TH1F*)inFile->Get("h_passEvtcut");
    TH1F* h_passTrkcut = (TH1F*)inFile->Get("h_passTrkcut");
    TH2F* h_Vx_Vy = (TH2F*)inFile->Get("h_Vx_Vy");

    TH3F* h_nHitsFit_Pt_Eta = (TH3F*)inFile->Get("h_nHitsFit_Pt_Eta");
    TH2F* h_nHitsFit_Pt = (TH2F*)h_nHitsFit_Pt_Eta->Project3D("zx");
    TH2F* h_nHitsFit_Eta = (TH2F*)h_nHitsFit_Pt_Eta->Project3D("zy");

    TH3F* h_nHitsDEdx_Pt_Eta = (TH3F*)inFile->Get("h_nHitsDEdx_Pt_Eta");
    TH2F* h_nHitsDEdx_Pt = (TH2F*)h_nHitsDEdx_Pt_Eta->Project3D("zx");
    TH2F* h_nHitsDEdx_Eta = (TH2F*)h_nHitsDEdx_Pt_Eta->Project3D("zy");

    TH3F* h_pDca_Pt_Eta = (TH3F*)inFile->Get("h_pDca_Pt_Eta");
    TH2F* h_pDca_Pt = (TH2F*)h_pDca_Pt_Eta->Project3D("xy");
    TH2F* h_pDca_Eta = (TH2F*)h_pDca_Pt_Eta->Project3D("xz");

    TH1F* h_ppT  = (TH1F*)inFile->Get("h_ppT");
    TH1F* h_pEta = (TH1F*)inFile->Get("h_pEta");
    TH1F* h_pPhi = (TH1F*)inFile->Get("h_pPhi");
    TH2F* h_pP_ppT = (TH2F*)inFile->Get("h_pP_ppT");
    TH2F* h_ppTc_pEta = (TH2F*)inFile->Get("h_ppTc_pEta");
	TH2F* h_dEdx_Pc = (TH2F*)inFile->Get("h_dEdx_Pc");
	TH2F* h_m2 = (TH2F*)inFile->Get("h_m2");

    TH2F* h_nSigmaElectron_P = (TH2F*)inFile->Get("h_nSigmaElectron_P");
    TH2F* h_nSigmaEcorr_P = (TH2F*)inFile->Get("h_nSigmaEcorr_P");
    TH2F* h_nSigmaPion_P = (TH2F*)inFile->Get("h_nSigmaPion_P");
    TH2F* h_nSigmaPicorr_P = (TH2F*)inFile->Get("h_nSigmaPicorr_P");

	// group 1
	TH1F* h_pT__TOFMatch = (TH1F*)inFile->Get("h_pT__TOFMatch");
	TH1F* h_Eta__TOFMatch = (TH1F*)inFile->Get("h_Eta__TOFMatch");
	TH1F* h_Phi__TOFMatch = (TH1F*)inFile->Get("h_Phi__TOFMatch");
	TH2F* h_invBeta_P__TOFMatch = (TH2F*)inFile->Get("h_invBeta_P__TOFMatch");
	TH2F* h_nSigmaElectron_P__1 = (TH2F*)inFile->Get("h_nSigmaElectron_P__1");
	TH2F* h_nSigmaProton_P__1 = (TH2F*)inFile->Get("h_nSigmaProton_P__1");
	TH2F* h_nSigmaElectron_P__TOFMatch = (TH2F*)inFile->Get("h_nSigmaElectron_P__TOFMatch");
	TH2F* h_nSigmaElectron_P__PIDcut_1 = (TH2F*)inFile->Get("h_nSigmaElectron_P__PIDcut_1");
	// group 2
	TH2F* h_nSigmaElectron_P__2 = (TH2F*)inFile->Get("h_nSigmaElectron_P__2");
	TH2F* h_nSigmaPion_P__2 = (TH2F*)inFile->Get("h_nSigmaPion_P__2");
	TH2F* h_nSigmaKaon_P__2 = (TH2F*)inFile->Get("h_nSigmaKaon_P__2");
	TH2F* h_nSigmaProton_P__2 = (TH2F*)inFile->Get("h_nSigmaProton_P__2");
	TH2F* h_nSigmaElectron_P__PIDcut_2 = (TH2F*)inFile->Get("h_nSigmaElectron_P__PIDcut_2");
	// group 3
	TH2F* h_nSigmaElectron_P__3 = (TH2F*)inFile->Get("h_nSigmaElectron_P__3");
	TH2F* h_nSigmaPion_P__3 = (TH2F*)inFile->Get("h_nSigmaPion_P__3");
	TH2F* h_nSigmaKaon_P__3 = (TH2F*)inFile->Get("h_nSigmaKaon_P__3");
	TH2F* h_nSigmaProton_P__3 = (TH2F*)inFile->Get("h_nSigmaProton_P__3");
	TH2F* h_nSigmaElectron_P__PIDcut_3 = (TH2F*)inFile->Get("h_nSigmaElectron_P__PIDcut_3");
	TH2F* h_nSigmaElectron_Eta__EIDcut_3 = (TH2F*)inFile->Get("h_nSigmaElectron_Eta__EIDcut_3");

	//sum
	TH1F* h_nSigmaElectron_P__EIDcut_total = (TH1F*)inFile->Get("h_nSigmaElectron_P__EIDcut_total");
	TH2F* h_eNumber_Cen = (TH2F*)inFile->Get("h_eNumber_Cen");
	TH1F* h_pT__electrons = (TH1F*)inFile->Get("h_pT__electrons");
	TH1F* h_eta__electrons = (TH1F*)inFile->Get("h_eta__electrons");
	TH1F* h_phi__electrons = (TH1F*)inFile->Get("h_phi__electrons");
	TH1F* h_pT__positrons = (TH1F*)inFile->Get("h_pT__positrons");
	TH1F* h_eta__positrons = (TH1F*)inFile->Get("h_eta__positrons");
	TH1F* h_phi__positrons = (TH1F*)inFile->Get("h_phi__positrons");
	TH1F* h_pT_ratio_EperP = (TH1F*)h_pT__electrons->Clone("h_pT_ratio_EperP");	h_pT_ratio_EperP->SetTitle("e^{-}/e^{+} p_{T} ratio;p_{T} [GeV/c];Ratio");	h_pT_ratio_EperP->Divide(h_pT__positrons);
	TH1F* h_eta_ratio_EperP = (TH1F*)h_eta__electrons->Clone("h_eta_ratio_EperP");	h_eta_ratio_EperP->SetTitle("e^{-}/e^{+} #eta ratio;#eta [GeV/c];Ratio");	h_eta_ratio_EperP->Divide(h_eta__positrons);
	TH1F* h_phi_ratio_EperP = (TH1F*)h_phi__electrons->Clone("h_phi_ratio_EperP");	h_phi_ratio_EperP->SetTitle("e^{-}/e^{+} #phi ratio;#phi [GeV/c];Ratio");	h_phi_ratio_EperP->Divide(h_phi__positrons);

    if (1)//nsigmaElectron__EID, h_pP_ppT, h_passEvtcut, h_passTrkcut
    {
	    //误差条，设置Marker形状颜色，设置线条颜色，设置图例，设置坐标轴标题，设置对数Y坐标）
	    h_Vx_Vy->SetTitle("V_{z} vs V_{x};V_{x} (cm);V_{z} (cm)");
	    TCanvas *c_QA = new TCanvas("c_QA", "c_QA", 200, 100, 1300, 1200);
	    c_QA->Divide(4, 4);

	    c_QA->cd(1);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
		h_dEdx_Pc->DrawClone("col z");

	    c_QA->cd(2);
	    gPad->SetLogy(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
		h_m2->DrawClone("col z");

	    c_QA->cd(3);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);

	    c_QA->cd(4);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);

	    c_QA->cd(5);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_Vx_Vy->DrawClone("col z");

	    c_QA->cd(6);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);


	    c_QA->cd(7);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_passEvtcut->DrawClone("TEXT");
	    h_passEvtcut->DrawClone("same");

	    c_QA->cd(8);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_passTrkcut->DrawClone("TEXT0");
	    h_passTrkcut->DrawClone("same");

	    c_QA->cd(9);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nHitsFit_Pt->SetTitle("nHitsFit vs p_{T}");
	    h_nHitsFit_Pt->GetXaxis()->SetRangeUser(0.0,3.0);
	    h_nHitsFit_Pt->DrawClone("col z");

	    c_QA->cd(10);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nHitsFit_Eta->SetTitle("nHitsFit vs #eta");
	    h_nHitsFit_Eta->DrawClone("col z");

	    c_QA->cd(11);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nHitsDEdx_Pt->SetTitle("nHitsDEdx vs p_{T}");
	    h_nHitsDEdx_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
	    h_nHitsDEdx_Pt->DrawClone("col z");

	    c_QA->cd(12);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nHitsDEdx_Eta->SetTitle("nHitsDEdx vs #eta");
	    h_nHitsDEdx_Eta->DrawClone("col z");

	    c_QA->cd(13);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);


	    c_QA->cd(14);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_pDca_Pt->SetTitle("DCA vs p_{T}");
	    h_pDca_Pt->GetXaxis()->SetRangeUser(0.0, 3.0);
	    h_pDca_Pt->DrawClone("col z");

	    c_QA->cd(15);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_pDca_Eta->SetTitle("DCA vs #eta");
	    h_pDca_Eta->DrawClone("col z");

	    c_QA->cd(16);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);

	    c_QA->SaveAs(Form("roots/%d_Track_and_Event_Check.png", number));
    }
	if (0) 
	{
		TCanvas *c_sum = new TCanvas("c_sum", "c_sum", 200, 100, 1000, 900);
		c_sum->Divide(2, 2);

		c_sum->cd(1);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P->GetYaxis()->SetRangeUser(-10, 10);
		h_nSigmaElectron_P->Draw("col z");

		c_sum->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__EIDcut_total->Draw("col z");

		c_sum->SaveAs(Form("roots/%d_EIDsum.png", number));
	}
    if (1)//track QA and TOF track check
    {
	    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800);
	    c1->Divide(3, 2);
	    c1->cd(1);
	    gPad->SetLogy(1);
	    h_ppT->GetXaxis()->SetRangeUser(0., 3.0);
	    //h_ppT->GetYaxis()->SetRangeUser(1., 1e8);
	    h_ppT->Draw();

	    c1->cd(2);
	    gPad->SetLogy(1);
	    //h_pEta->GetYaxis()->SetRangeUser(1., 1e8);
	    h_pEta->Draw();

	    c1->cd(3);
	    gPad->SetLogy(1);
	    h_pPhi->GetXaxis()->SetRangeUser(-TMath::Pi() * 4 / 3, TMath::Pi() * 4 / 3);
	    //h_pPhi->GetYaxis()->SetRangeUser(1., 1e8);
	    h_pPhi->GetXaxis()->SetNdivisions(-108);    // 负号 = 用户自定义，202 = π 刻度
	    h_pPhi->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "-4#pi/3"); // 让 ROOT 自动写 kπ
	    h_pPhi->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "-#pi");
	    h_pPhi->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "-2#pi/3");
	    h_pPhi->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "-#pi/3");
	    h_pPhi->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "0");
	    h_pPhi->GetXaxis()->ChangeLabel(-4, -1, -1, -1, -1, -1, "#pi/3");
	    h_pPhi->GetXaxis()->ChangeLabel(-3, -1, -1, -1, -1, -1, "2#pi/3");
	    h_pPhi->GetXaxis()->ChangeLabel(-2, -1, -1, -1, -1, -1, "#pi");
	    h_pPhi->GetXaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, "4#pi/3"); // 让 ROOT 自动写 kπ
	    h_pPhi->Draw();

	    c1->cd(4);
	    gPad->SetLogy(1);
	    h_pT__TOFMatch->GetXaxis()->SetRangeUser(0., 3.0);
	    //h_pT__TOFMatch->GetYaxis()->SetRangeUser(1., 1e8);
	    h_pT__TOFMatch->Draw();

	    c1->cd(5);
	    gPad->SetLogy(1);
	    //h_Eta__TOFMatch->GetYaxis()->SetRangeUser(1., 1e8);
	    h_Eta__TOFMatch->Draw();

	    c1->cd(6);
	    gPad->SetLogy(1);
	    h_Phi__TOFMatch->GetXaxis()->SetRangeUser(-TMath::Pi() * 4 / 3, TMath::Pi() * 4 / 3);
	    //h_Phi__TOFMatch->GetYaxis()->SetRangeUser(1., 1e8);
	    h_Phi__TOFMatch->GetXaxis()->SetNdivisions(-108);    // 负号 = 用户自定义，202 = π 刻度
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, "-4#pi/3"); // 让 ROOT 自动写 kπ
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(2, -1, -1, -1, -1, -1, "-#pi");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(3, -1, -1, -1, -1, -1, "-2#pi/3");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(4, -1, -1, -1, -1, -1, "-#pi/3");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(5, -1, -1, -1, -1, -1, "0");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-4, -1, -1, -1, -1, -1, "#pi/3");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-3, -1, -1, -1, -1, -1, "2#pi/3");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-2, -1, -1, -1, -1, -1, "#pi");
	    h_Phi__TOFMatch->GetXaxis()->ChangeLabel(-1, -1, -1, -1, -1, -1, "4#pi/3"); // 让 ROOT 自动写 kπ
	    h_Phi__TOFMatch->Draw();

	    c1->SaveAs(Form("roots/%d_TrackQA_and_TrackTOFMatch.png", number));
    }
	if (1)//e+ e- track QA
	{
		h_pT__positrons->SetLineColor(kRed);
		h_eta__positrons->SetLineColor(kRed);
		h_phi__positrons->SetLineColor(kRed);

		TCanvas *c_eTrack = new TCanvas("c_eTrack", "c_eTrack", 1200, 800);
		c_eTrack->Divide(3, 2);

		c_eTrack->cd(1);
		gPad->SetLogy(1);
		TLegend *leg1 = new TLegend(0.65, 0.75, 0.88, 0.88);
		leg1->AddEntry(h_pT__electrons, "e^{-} (electrons)", "l");
		leg1->AddEntry(h_pT__positrons, "e^{+} (positrons)", "l");
		leg1->SetBorderSize(0);
		leg1->SetFillStyle(0);
		h_pT__electrons->GetXaxis()->SetRangeUser(0., 3.0);
		h_pT__electrons->Draw();
		h_pT__positrons->Draw("same");
		leg1->Draw("same");

		c_eTrack->cd(2);
		gPad->SetLogy(1);
		h_eta__electrons->GetXaxis()->SetRangeUser(-2.0, 0.0);
		h_eta__electrons->Draw();
		h_eta__positrons->Draw("same");

		c_eTrack->cd(3);
		gPad->SetLogy(1);
		h_phi__electrons->GetXaxis()->SetRangeUser(-3.2, 3.2);
		h_phi__electrons->Draw();
		h_phi__positrons->Draw("same");

		c_eTrack->cd(4);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_pT_ratio_EperP->SetMaximum(2);
		h_pT_ratio_EperP->SetMinimum(0.5);
		h_pT_ratio_EperP->Draw("colz");

		c_eTrack->cd(5);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_eta_ratio_EperP->SetMaximum(2);
		h_eta_ratio_EperP->SetMinimum(0.5);
		h_eta_ratio_EperP->Draw("colz");

		c_eTrack->cd(6);
		gPad->SetLogy(0);
		gStyle->SetOptStat(0);
		h_phi_ratio_EperP->SetMaximum(2);
		h_phi_ratio_EperP->SetMinimum(0.5);
		h_phi_ratio_EperP->Draw("colz");

		c_eTrack->SaveAs(Form("roots/%d_ep_TrackQA.png", number));
	}

    if (0)// EID in group1(pT>0.2, -1.4<eta<0)
    {
	    // 设置直方图格式
	    // 去除误差条，设置Marker形状颜色，设置线条颜色，设置图例，设置坐标轴标题，设置对数Y坐标）

	    //画图
	    TCanvas *c_g1 = new TCanvas("c_g1", "c_g1", 1200, 800);
	    c_g1->Divide(3, 2);
	    c_g1->cd(1);//nSigmaE
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nSigmaElectron_P__1->Draw("col z");

	    c_g1->cd(2);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_invBeta_P__TOFMatch->Draw("col z");
		TLine *line1_beta1 = new TLine(0, 1 - 0.025, 5, 1 - 0.025);	line1_beta1->SetLineColor(kRed); line1_beta1->Draw("same");
		TLine *line1_beta2 = new TLine(0, 1 + 0.025, 5, 1 + 0.025);	line1_beta2->SetLineColor(kRed); line1_beta2->Draw("same");
		TF1 *invBeta_pi    = new TF1("invBeta_pi", "sqrt(x*x+[0]*[0])/x", 0., 5.);   invBeta_pi->SetParameters(0.139);   invBeta_pi->SetLineColor(kRed);   invBeta_pi->Draw("same");
		TF1 *invBeta_k     = new TF1("invBeta_k", "sqrt(x*x+[0]*[0])/x", 0., 5.);    invBeta_k->SetParameters(0.494);    invBeta_k->SetLineColor(kRed);    invBeta_k->Draw("same");
		TF1 *invBeta_p     = new TF1("invBeta_p", "sqrt(x*x+[0]*[0])/x", 0., 5.);    invBeta_p->SetParameters(0.938);    invBeta_p->SetLineColor(kRed);    invBeta_p->Draw("same");
		TF1 *invBeta_pn = new TF1("invBeta_pn", "sqrt(x*x+[0]*[0])/x", 0., 5.);   invBeta_pn->SetParameters(1.876);   invBeta_pn->SetLineColor(kRed);   invBeta_pn->Draw("same");
		TF1 *invBeta_pnn = new TF1("invBeta_pnn", "sqrt(x*x+[0]*[0])/x", 0., 5.);  invBeta_pnn->SetParameters(2.814);  invBeta_pnn->SetLineColor(kRed);  invBeta_pnn->Draw("same");
		TF1 *invBeta_ppn = new TF1("invBeta_ppn", "sqrt((x*2)*(x*2)+[0]*[0])/(x*2)", 0., 3.);  invBeta_ppn->SetParameters(2.814);  invBeta_ppn->SetLineColor(kRed);  invBeta_ppn->Draw("same");
		TF1 *invBeta_ppnn = new TF1("invBeta_ppnn", "sqrt((x*2)*(x*2)+[0]*[0])/(x*2)", 0., 3.); invBeta_ppnn->SetParameters(3.762); invBeta_ppnn->SetLineColor(kRed); invBeta_ppnn->Draw("same");

	    c_g1->cd(3);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nSigmaElectron_P__TOFMatch->Draw("col z");

	    c_g1->cd(4);
	    gPad->SetLogz(1);
	    gPad->SetLeftMargin(0.12);
	    gPad->SetRightMargin(0.12);
	    gStyle->SetOptStat(0);
	    h_nSigmaProton_P__1->Draw("col z");
		TLine *line1_P1 = new TLine(0.0, 5.0, 5.0, 0.0);	line1_P1->SetLineColor(kRed); line1_P1->Draw("same");
		TLine *line1_P2 = new TLine(0.0, -2.0, 5.0, -2.0);	line1_P2->SetLineColor(kRed); line1_P2->Draw("same");

		c_g1->cd(5);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__PIDcut_1->Draw("col z");
		TLine *line1_TOF1 = new TLine(0.0, 3.0, 1.6, 3.0); line1_TOF1->SetLineColor(kRed); line1_TOF1->Draw("same");
		TLine *line1_TOF3 = new TLine(0.0, -4.5, 1.0, -1.0); line1_TOF3->SetLineColor(kRed); line1_TOF3->Draw("same");
		TLine *line1_TOF2 = new TLine(1.0, -1.0, 1.6, -1.0); line1_TOF2->SetLineColor(kRed); line1_TOF2->Draw("same");
		TLine *line1_TOF4 = new TLine(1.6, -1.0, 1.6, 3.0); line1_TOF4->SetLineColor(kRed); line1_TOF4->Draw("same");

	    c_g1->SaveAs(Form("roots/%d_group1.png", number));
    }
	if (0)//EID in group2(pT>0.2, -2.5<eta<-1.4)
	{
		// 设置直方图格式
		// 去除误差条，设置Marker形状颜色，设置线条颜色，设置图例，设置坐标轴标题，设置对数Y坐标）

		//画图
		TCanvas *c_g2 = new TCanvas("c_g2", "c_g2", 1200, 800);
		c_g2->Divide(3, 2);
		c_g2->cd(1);//nSigmaE
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__2->Draw("col z");

		c_g2->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaPion_P__2->Draw("col z");
		TLine *line2_Pi1 = new TLine(0.0, 3.0, 5.0, 3.0);	line2_Pi1->SetLineColor(kRed); line2_Pi1->Draw("same");
		TLine *line2_Pi2 = new TLine(0.0, -4.0, 5.0, -4.0);	line2_Pi2->SetLineColor(kRed); line2_Pi2->Draw("same");

		c_g2->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaKaon_P__2->Draw("col z");
		TLine *line2_K1 = new TLine(0.0, 2.0, 5.0, 2.0);	line2_K1->SetLineColor(kRed); line2_K1->Draw("same");
		TLine *line2_K2 = new TLine(0.0, -2.0, 5.0, -2.0);	line2_K2->SetLineColor(kRed); line2_K2->Draw("same");

		c_g2->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaProton_P__2->Draw("col z");
		TLine *line2_P1 = new TLine(0.0, 6.0, 5.0, 0.0);	line2_P1->SetLineColor(kRed); line2_P1->Draw("same");
		TLine *line2_P2 = new TLine(0.0, -3.0, 5.0, -3.0);	line2_P2->SetLineColor(kRed); line2_P2->Draw("same");

		c_g2->cd(5);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__PIDcut_2->Draw("col z");
		TLine *line2_E1 = new TLine(0.0, 2.0, 1.5, 2.0);	line2_E1->SetLineColor(kRed); line2_E1->Draw("same");
		TLine *line2_E2 = new TLine(0.0, -2.0, 1.5, -2.0);	line2_E2->SetLineColor(kRed); line2_E2->Draw("same");
		TLine *line2_E3 = new TLine(1.5, -2.0, 1.5, 2.0);	line2_E3->SetLineColor(kRed); line2_E3->Draw("same");

		c_g2->SaveAs(Form("roots/%d_group2.png", number));
	}
	if (0)//EID in group3(pT<0.2, eta<0.0)
	{
		// 设置直方图格式
		// 去除误差条，设置Marker形状颜色，设置线条颜色，设置图例，设置坐标轴标题，设置对数Y坐标）
		h_nSigmaElectron_P__3->GetXaxis()->SetRangeUser(0., 1.);
		h_nSigmaPion_P__3->GetXaxis()->SetRangeUser(0., 1.);
		h_nSigmaKaon_P__3->GetXaxis()->SetRangeUser(0., 1.);
		h_nSigmaProton_P__3->GetXaxis()->SetRangeUser(0., 1.);
		h_nSigmaElectron_P__PIDcut_3->GetXaxis()->SetRangeUser(0., 1.);
		//画图
		TCanvas *c_g3 = new TCanvas("c_g3", "c_g3", 1200, 800);
		c_g3->Divide(3, 2);
		c_g3->cd(1);//nSigmaE
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__3->Draw("col z");

		c_g3->cd(2);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaPion_P__3->Draw("col z");
		TLine *line3_Pi1 = new TLine(0.0, 3.0, 1.0, 3.0);	line3_Pi1->SetLineColor(kRed); line3_Pi1->Draw("same");
		TLine *line3_Pi2 = new TLine(0.0, -4.0, 1.0, -4.0);	line3_Pi2->SetLineColor(kRed); line3_Pi2->Draw("same");

		c_g3->cd(3);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaKaon_P__3->Draw("col z");

		c_g3->cd(4);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaProton_P__3->Draw("col z");

		c_g3->cd(5);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_P__PIDcut_3->Draw("col z");
		TLine *line3_E1 = new TLine(0.0, 2.5, 0.7, 2.5);	line3_E1->SetLineColor(kRed); line3_E1->Draw("same");
		TLine *line3_E2 = new TLine(0.0, -4.0, 0.7, -4.0);	line3_E2->SetLineColor(kRed); line3_E2->Draw("same");
		TLine *line3_E3 = new TLine(0.7, -4.0, 0.7, 2.5);	line3_E3->SetLineColor(kRed); line3_E3->Draw("same");

		c_g3->cd(6);
		gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h_nSigmaElectron_Eta__EIDcut_3->Draw("col z");

		c_g3->SaveAs(Form("roots/%d_group3.png", number));
	}
}

