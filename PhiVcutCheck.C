//---------------------从root中提取直方图，并进行设置更改和元素添加----------------------
#include "someFunction.h"
//15_20251204_nsigmaE_3regions;16_20251229_TOFElectron_and_EcorrPicorr_P24ia;17_20251221_TOFElectron_and_EcorrPicorr_P24iy;18_20251230_iTPCElectron_and_EcorrPicorr_P24iy;20_20260108_iTPCElectron_goodCut4LowHighP_rmLowP_region3;21_20260112_iTPCElectron_goodCut4LowHighP_withLowP_region3;22_20260113_iTPCElectron_goodCut4LowHighP_rmLowEta_region3
void PhiVcutCheck(TString inFilename = "roots/11_20260311_FixAu_group1_Mixed_1_500_50_cutPhiVByPair.root", Int_t number = 11)
{
   // 从root文件中导入待拟合的直方图
   TFile *inFile = new TFile(inFilename); if (!inFile)
   {
	   cout << "The input file is not found! Exiting..." << endl;
	   return;
   }
   TH2F* h_Mee_PhiV__ranComb		 = (TH2F*)inFile->Get("h_Mee_PhiV__unlikeSame");//h_Mee_PhiV__ranComb
   TH1F* h_Mee__ranComb				 = (TH1F*)inFile->Get("h_Mee__unlikeSame");//h_Mee__ranComb
   TH1F* h_Mee__ranComb__w_PhiV_Cut  = (TH1F*)inFile->Get("h_Mee__unlikeSame__w_PhiV_Cut");//h_Mee__ranComb__w_PhiV_Cut

   auto fphiVcut = new TF1("fphiVcut", "0.84326*exp(-49.4819*x)-0.996609*x+0.19801", 0, 1.0);
   fphiVcut->SetMaximum(1.0);
   fphiVcut->SetNpx(1000);
   if (1)//check PhiV
   {
	   h_Mee__ranComb->SetLineColor(kBlack);
	   h_Mee__ranComb->GetYaxis()->SetTitleOffset(1.5);
	   h_Mee__ranComb__w_PhiV_Cut->SetLineColor(kBlue);
	   h_Mee__ranComb->SetTitle("M_{ee} w/wo #phi_{V} cut;M_{ee} (GeV/c^{2});Counts");

	   TCanvas* c_temp = new TCanvas("c_temp", "c_temp", 800, 800);
	   c_temp->Divide(2, 2);
	   c_temp->cd(1);
	   gPad->SetLogz(1);
	   gPad->SetLeftMargin(0.12);
	   gPad->SetRightMargin(0.12);
	   gStyle->SetOptStat(0);
	   h_Mee_PhiV__ranComb->GetXaxis()->SetRangeUser(0.0, 0.4);
	   h_Mee_PhiV__ranComb->GetYaxis()->SetRangeUser(0.0, 1.0);
	   fphiVcut->GetXaxis()->SetRangeUser(0, 0.2);
	   h_Mee_PhiV__ranComb->DrawClone("col z");
	   fphiVcut->DrawClone("same");

	   c_temp->cd(2);
	   gPad->SetLogz(1);
	   gPad->SetLeftMargin(0.12);
	   gPad->SetRightMargin(0.12);
	   gStyle->SetOptStat(0);
	   h_Mee__ranComb->GetXaxis()->SetRangeUser(0, 0.2);
	   h_Mee__ranComb__w_PhiV_Cut->GetXaxis()->SetRangeUser(0, 0.2);
	   h_Mee__ranComb->DrawClone("");
	   h_Mee__ranComb__w_PhiV_Cut->DrawClone("same");
	   auto legend = new TLegend(0.53, 0.58, 0.78, 0.68);
	   legend->SetFillColor(0); legend->SetBorderSize(0);
	   legend->AddEntry(h_Mee__ranComb, "\t without #phi_{V} cut", "lp");
	   legend->AddEntry(h_Mee__ranComb__w_PhiV_Cut, "\t with #phi_{V} cut", "lp");
	   legend->SetMargin(0.20);
	   gStyle->SetLegendTextSize(0.04);
	   legend->DrawClone("same");

	   c_temp->cd(3);
	   gPad->SetLogz(1);
	   gPad->SetLeftMargin(0.12);
	   gPad->SetRightMargin(0.12);
	   gStyle->SetOptStat(0);
	   h_Mee__ranComb->GetXaxis()->SetRangeUser(0, 1);
	   h_Mee__ranComb->DrawClone("");

	   c_temp->cd(4);
	   gPad->SetLogz(1);
	   gPad->SetLeftMargin(0.12);
	   gPad->SetRightMargin(0.12);
	   gStyle->SetOptStat(0);

	   c_temp->SaveAs(Form("roots/%d_PhiV_Check.png", number));
   }

}

