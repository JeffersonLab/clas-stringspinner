// ninja && clas-stringspinner --seed 29877 --trig 10000 --cut-inclusive 11,211,-211|grep M_X| awk '{print $2}' > output.dat
void draw_mx_2(TString file_name = "output.both.dat") {
  gStyle->SetOptStat(0);
  auto tr = new TTree();
  tr->ReadFile(file_name, "mx_before/D:mx_after/D");
  // new TCanvas();
  // tr->Draw("mx");
  auto c = new TCanvas("c", "c", 800, 600);
  auto mx_dist_before = new TH1D("mx_dist_before", "M_{X} distribution", 60, 1.0, 1.35);
  auto mx_dist_after = new TH1D("mx_dist_after",   "M_{X} distribution", 60, 1.0, 1.35);
  tr->Project("mx_dist_before", "mx_before");
  tr->Project("mx_dist_after", "mx_after");
  mx_dist_before->SetLineColor(kBlue);
  mx_dist_after->SetLineColor(kRed);
  mx_dist_after->SetLineStyle(kDashed);
  mx_dist_after->SetLineWidth(2);
  mx_dist_before->SetLineWidth(2);
  mx_dist_before->Draw();
  mx_dist_after->Draw("same");
  auto delta = new TLine(1.232, 0, 1.232, mx_dist_before->GetMaximum());
  delta->SetLineColor(kGreen+1);
  delta->SetLineWidth(3);
  delta->Draw();
  c->SaveAs(file_name+".png");
}
