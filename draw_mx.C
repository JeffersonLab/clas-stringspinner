// ninja && clas-stringspinner --seed 29877 --trig 10000 --cut-inclusive 11,211,-211|grep M_X| awk '{print $2}' > output.dat
void draw_mx(TString file_name = "output.dat", TString title = "m_x") {
  auto tr = new TTree();
  tr->ReadFile(file_name, "mx/D");
  // new TCanvas();
  // tr->Draw("mx");
  auto c = new TCanvas("c", "c", 800, 600);
  auto mx_dist = new TH1D("mx_dist",title,100,0.9,1.5);
  tr->Project("mx_dist", "mx");
  mx_dist->Draw();
  auto peak = mx_dist->GetBinCenter(mx_dist->GetMaximumBin());
  auto dev = 0.03;
  // mx_dist->Fit("gaus","","",peak-dev,peak+dev);
  auto delta = new TLine(1.232, 0, 1.232, mx_dist->GetMaximum());
  delta->SetLineColor(kRed);
  delta->SetLineWidth(3);
  delta->Draw();
  c->SaveAs(file_name+".png");
}
