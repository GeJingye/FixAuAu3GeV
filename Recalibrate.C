//---------------------从AuAU.root中提取直方图，进行n\sigma_{e}重新标定的计算----------------------
#include "someFunction.h"
#include <iomanip>
void Recalibrate(TString inFileName = "roots/11_20260311_FixAu_group1_Mixed_1_500_50_cutPhiVByPair.root", Int_t number = 11)//
{
	// 从root文件中导入待拟合的直方图
	TFile *inFile = new TFile(inFileName);
	if (!inFile) { cout << "The input file is not found! Exiting..." << endl; return; }
	// 导入直方图（from PureE）
	// 选择纯度(purity)较高的电子，例如组成电子对Mee<0.015的电子。
	// TH3F* h_Pt_Eta_Phi__PureE	   = (TH3F*)inFile->Get("h_Pt_Eta_Phi__PureE");
	TH3F* h_Pt_Cen_nSigmaE__PureE = (TH3F*)inFile->Get("h_Pt_Cen_nSigmaE__PureE"); h_Pt_Cen_nSigmaE__PureE->RebinX(10);
	TH3F* h_Eta_Cen_nSigmaE__PureE = (TH3F*)inFile->Get("h_Eta_Cen_nSigmaE__PureE"); h_Eta_Cen_nSigmaE__PureE->RebinX(10);
	TH3F* h_Phi_Cen_nSigmaE__PureE = (TH3F*)inFile->Get("h_Phi_Cen_nSigmaE__PureE"); h_Phi_Cen_nSigmaE__PureE->RebinX(10);
	//TH3F* h_Pt_Cen_nSigmaPi__PurePi = (TH3F*)inFile->Get("h_Pt_Cen_nSigmaPi__PurePi");
	//TH3F* h_Eta_Cen_nSigmaPi__PurePi = (TH3F*)inFile->Get("h_Eta_Cen_nSigmaPi__PurePi");
	//TH3F* h_Phi_Cen_nSigmaPi__PurePi = (TH3F*)inFile->Get("h_Phi_Cen_nSigmaPi__PurePi");
	//TH3F* h_Pt_Eta_Phi__PurePi = (TH3F*)inFile->Get("h_Pt_Eta_Phi__PurePi");
	if (!h_Pt_Cen_nSigmaE__PureE || !h_Eta_Cen_nSigmaE__PureE || !h_Phi_Cen_nSigmaE__PureE) { cout << "Some histograms are not found! Exiting..." << endl; return; }

	//通过函数计算均值分布
	//Pt
	TH2F* h2_zmean_binCount_Pt = (TH2F*)MeanOfH3D_binCount(h_Pt_Cen_nSigmaE__PureE, "h2_zmean_binCount_Pt");
	TH2F* h2_zmean_gausFit_Pt = (TH2F*)MeanOfH3D_gausFit(h_Pt_Cen_nSigmaE__PureE, "h2_zmean_gausFit_Pt");
	TH1F* h1_xmean_cen__binCount_Pt = (TH1F*)Meanof2DAlongX(h2_zmean_binCount_Pt, "h1_xmean_cen__binCount_Pt");
	TH1F* h1_xmean_cen__gausFit_Pt = (TH1F*)Meanof2DAlongX(h2_zmean_gausFit_Pt, "h1_xmean_cen__gausFit_Pt");
	//Eta
	TH2F* h2_zmean_binCount_Eta = (TH2F*)MeanOfH3D_binCount(h_Eta_Cen_nSigmaE__PureE, "h2_zmean_binCount_Eta");
	TH2F* h2_zmean_gausFit_Eta = (TH2F*)MeanOfH3D_gausFit(h_Eta_Cen_nSigmaE__PureE, "h2_zmean_gausFit_Eta");
	TH1F* h1_xmean_cen__binCount_Eta = (TH1F*)Meanof2DAlongX(h2_zmean_binCount_Eta, "h1_xmean_cen__binCount_Eta");
	TH1F* h1_xmean_cen__gausFit_Eta = (TH1F*)Meanof2DAlongX(h2_zmean_gausFit_Eta, "h1_xmean_cen__gausFit_Eta");

	//Phi
	TH2F* h2_zmean_binCount_Phi = (TH2F*)MeanOfH3D_binCount(h_Phi_Cen_nSigmaE__PureE, "h2_zmean_binCount_Phi");
	TH2F* h2_zmean_gausFit_Phi = (TH2F*)MeanOfH3D_gausFit(h_Phi_Cen_nSigmaE__PureE, "h2_zmean_gausFit_Phi");
	TH1F* h1_xmean_cen__binCount_Phi = (TH1F*)Meanof2DAlongX(h2_zmean_binCount_Phi, "h1_xmean_cen__binCount_Phi");
	TH1F* h1_xmean_cen__gausFit_Phi = (TH1F*)Meanof2DAlongX(h2_zmean_gausFit_Phi, "h1_xmean_cen__gausFit_Phi");

	//输出2D-eta-cen平均值
	std::cout << "mean[cen][eta]:" << endl;
	for (int iy = 1; iy <= h2_zmean_gausFit_Eta->GetNbinsY(); iy++)//h2_zmean_gausFit_Eta;
	{
		std::cout << "{";
		for (int ix = 1; ix <= h2_zmean_gausFit_Eta->GetNbinsX(); ix++)
		{
			float content = h2_zmean_gausFit_Eta->GetBinContent(ix, iy);
			std::cout << std::fixed << std::setprecision(3);
			std::cout << content << ',';
		}
		std::cout << "}," << endl;
	}
	//输出2D-phi-cen平均值
	std::cout << "mean[cen][phi]:" << endl;
	for (int iy = 1; iy <= h2_zmean_gausFit_Phi->GetNbinsY(); iy++)//h2_zmean_gausFit_Phi;
	{
		std::cout << "{";
		for (int ix = 1; ix <= h2_zmean_gausFit_Phi->GetNbinsX(); ix++)
		{
			float content = h2_zmean_gausFit_Phi->GetBinContent(ix, iy);
			std::cout << std::fixed << std::setprecision(3);
			std::cout << content << ',';
		}
		std::cout << "}," << endl;
	}
	//输出1D-phi-cen平均值
	std::cout << "phimean[cen]:" << endl;
	for (int ix = 1; ix <= h1_xmean_cen__gausFit_Phi->GetNbinsX(); ix++)
	{
		float content = h1_xmean_cen__gausFit_Phi->GetBinContent(ix);
		std::cout << content << "," << endl;
		if(ix== h1_xmean_cen__gausFit_Phi->GetNbinsX())cout<<endl;
	}

	if (1)//nsigma_e vs pT
	{
		//画图
		TCanvas *c_Pt = new TCanvas("c_Pt", "c_Pt", 1400, 600);
		c_Pt->Divide(4, 2);
		c_Pt->cd(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h2_zmean_binCount_Pt->SetMaximum(1);
		h2_zmean_binCount_Pt->SetMinimum(-1.0);
		h2_zmean_binCount_Pt->DrawClone("col z");

		c_Pt->cd(2);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h2_zmean_gausFit_Pt->SetMaximum(1);
		h2_zmean_gausFit_Pt->SetMinimum(-1.0);
		h2_zmean_gausFit_Pt->DrawClone("col z");

		c_Pt->cd(3);
		//gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h1_xmean_cen__binCount_Pt->DrawClone("");

		c_Pt->cd(4);
		//gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h1_xmean_cen__gausFit_Pt->DrawClone("");

		c_Pt->cd(5);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		TLegend *leg_E = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg_E->SetBorderSize(1);
		leg_E->SetFillColor(0);
		TH1F* proj1_bC = (TH1F*)h2_zmean_binCount_Pt->ProjectionX("_1", 1, 1);
		proj1_bC->SetLineColor(1);
		proj1_bC->DrawClone();
		leg_E->AddEntry(proj1_bC, "70~80%");
		for (int iy = 2; iy <= h2_zmean_binCount_Pt->GetYaxis()->GetNbins(); iy++)
		{
			TH1F* proj_bC = (TH1F*)h2_zmean_binCount_Pt->ProjectionX(Form("_%d", iy), iy, iy);
			proj_bC->SetLineColor(iy);
			leg_E->AddEntry(proj_bC, Form("%d~%d%", 80 - 10 * iy, 90 - 10 * iy), "l");
			proj_bC->DrawClone("same");
		}
		leg_E->DrawClone("same");

		c_Pt->cd(6);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		TLegend *leg_Ecorr = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg_Ecorr->SetBorderSize(1);
		leg_Ecorr->SetFillColor(0);
		TH1F* proj1_gF = (TH1F*)h2_zmean_gausFit_Pt->ProjectionX("_1", 1, 1);
		proj1_gF->SetLineColor(1);
		proj1_gF->DrawClone();
		leg_Ecorr->AddEntry(proj1_gF, "70~80%");
		for (int iy = 2; iy <= h2_zmean_gausFit_Pt->GetYaxis()->GetNbins(); iy++)
		{
			TH1F* proj_gF = (TH1F*)h2_zmean_gausFit_Pt->ProjectionX(Form("_%d", iy), iy, iy);
			proj_gF->SetLineColor(iy);
			leg_Ecorr->AddEntry(proj_gF, Form("%d~%d%", 80 - 10 * iy, 90 - 10 * iy), "l");
			proj_gF->DrawClone("same");
		}
		leg_Ecorr->DrawClone("same");

		c_Pt->SaveAs(Form("roots/%d_nsigmaMean_Pt.png", number));
	}
	if (1)//nsigma_e vs Eta
	{
		//画图
		TCanvas *c_Eta = new TCanvas("c_Eta", "c_Eta", 1400, 600);
		c_Eta->Divide(4, 2);
		c_Eta->cd(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h2_zmean_binCount_Eta->SetMaximum(1);
		h2_zmean_binCount_Eta->SetMinimum(-1.0);
		h2_zmean_binCount_Eta->DrawClone("col z");

		c_Eta->cd(2);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h2_zmean_gausFit_Eta->SetMaximum(1);
		h2_zmean_gausFit_Eta->SetMinimum(-1.0);
		h2_zmean_gausFit_Eta->DrawClone("col z");

		c_Eta->cd(3);
		//gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h1_xmean_cen__binCount_Eta->DrawClone("");

		c_Eta->cd(4);
		//gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h1_xmean_cen__gausFit_Eta->DrawClone("");

		c_Eta->cd(5);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		TLegend *leg_E = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg_E->SetBorderSize(1);
		leg_E->SetFillColor(0);
		TH1F* proj1_bC = (TH1F*)h2_zmean_binCount_Eta->ProjectionX("_1", 1, 1);
		proj1_bC->SetLineColor(1);
		proj1_bC->DrawClone();
		leg_E->AddEntry(proj1_bC, "70~80%");
		for (int iy = 2; iy <= h2_zmean_binCount_Eta->GetYaxis()->GetNbins(); iy++)
		{
			TH1F* proj_bC = (TH1F*)h2_zmean_binCount_Eta->ProjectionX(Form("_%d", iy), iy, iy);
			proj_bC->SetLineColor(iy);
			leg_E->AddEntry(proj_bC, Form("%d~%d%", 80-10*iy,90-10*iy), "l");
			proj_bC->DrawClone("same");
		}
		leg_E->DrawClone("same");

		c_Eta->cd(6);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		TLegend *leg_Ecorr = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg_Ecorr->SetBorderSize(1);
		leg_Ecorr->SetFillColor(0);
		TH1F* proj1_gF = (TH1F*)h2_zmean_gausFit_Eta->ProjectionX("_1", 1, 1);
		proj1_gF->SetLineColor(1);
		proj1_gF->DrawClone();
		leg_Ecorr->AddEntry(proj1_gF, "70~80%");
		for (int iy = 2; iy <= h2_zmean_gausFit_Eta->GetYaxis()->GetNbins(); iy++)
		{
			TH1F* proj_gF = (TH1F*)h2_zmean_gausFit_Eta->ProjectionX(Form("_%d", iy), iy, iy);
			proj_gF->SetLineColor(iy);
			leg_Ecorr->AddEntry(proj_gF, Form("%d~%d%", 80 - 10 * iy, 90 - 10 * iy), "l");
			proj_gF->DrawClone("same");
		}
		leg_Ecorr->DrawClone("same");


		c_Eta->SaveAs(Form("roots/%d_nsigmaMean_Eta.png", number));
	}
	if (1)//nsigma_e vs Phi
	{
		//画图
		TCanvas *c_Phi = new TCanvas("c_Phi", "c_Phi", 1400, 600);
		c_Phi->Divide(4, 2);
		c_Phi->cd(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h2_zmean_binCount_Phi->SetMaximum(1);
		h2_zmean_binCount_Phi->SetMinimum(-1.0);
		h2_zmean_binCount_Phi->DrawClone("col z");

		c_Phi->cd(2);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h2_zmean_gausFit_Phi->SetMaximum(1);
		h2_zmean_gausFit_Phi->SetMinimum(-1.0);
		h2_zmean_gausFit_Phi->DrawClone("col z");

		c_Phi->cd(3);
		//gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h1_xmean_cen__binCount_Phi->DrawClone("");

		c_Phi->cd(4);
		//gPad->SetLogz(1);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);
		h1_xmean_cen__gausFit_Phi->DrawClone("");

		c_Phi->cd(5);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		TLegend *leg_E = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg_E->SetBorderSize(1);
		leg_E->SetFillColor(0);
		TH1F* proj1_bC = (TH1F*)h2_zmean_binCount_Phi->ProjectionX("_1", 1, 1);
		proj1_bC->SetLineColor(1);
		proj1_bC->DrawClone();
		leg_E->AddEntry(proj1_bC, "70~80%");
		for (int iy = 2; iy <= h2_zmean_binCount_Phi->GetYaxis()->GetNbins(); iy++)
		{
			TH1F* proj_bC = (TH1F*)h2_zmean_binCount_Phi->ProjectionX(Form("_%d", iy), iy, iy);
			proj_bC->SetLineColor(iy);
			leg_E->AddEntry(proj_bC, Form("%d~%d%", 80 - 10 * iy, 90 - 10 * iy), "l");
			proj_bC->DrawClone("same");
		}
		leg_E->DrawClone("same");

		c_Phi->cd(6);
		gPad->SetLeftMargin(0.12);
		gPad->SetRightMargin(0.12);
		gStyle->SetOptStat(0);

		TLegend *leg_Ecorr = new TLegend(0.7, 0.6, 0.9, 0.9);
		leg_Ecorr->SetBorderSize(1);
		leg_Ecorr->SetFillColor(0);
		TH1F* proj1_gF = (TH1F*)h2_zmean_gausFit_Phi->ProjectionX("_1", 1, 1);
		proj1_gF->SetLineColor(1);
		proj1_gF->DrawClone();
		leg_Ecorr->AddEntry(proj1_gF, "70~80%");
		for (int iy = 2; iy <= h2_zmean_gausFit_Phi->GetYaxis()->GetNbins(); iy++)
		{
			TH1F* proj_gF = (TH1F*)h2_zmean_gausFit_Phi->ProjectionX(Form("_%d", iy), iy, iy);
			proj_gF->SetLineColor(iy);
			leg_Ecorr->AddEntry(proj_gF, Form("%d~%d%", 80 - 10 * iy, 90 - 10 * iy), "l");
			proj_gF->DrawClone("same");
		}
		leg_Ecorr->DrawClone("same");


		c_Phi->SaveAs(Form("roots/%d_nsigmaMean_Phi.png", number));
	}	
}