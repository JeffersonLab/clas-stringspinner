// ninja && clas-stringspinner --seed 29877 --trig 10000 --cut-inclusive 11,211,-211|grep M_X| awk '{print $2}' > output.dat
void draw_mx(TString file_name = "output.dat") {
  auto tr = new TTree();
  tr->ReadFile(file_name, "mx/D");
  // new TCanvas();
  // tr->Draw("mx");
  auto c = new TCanvas("c", "c", 800, 600);
  auto mx_dist = new TH1D("mx_dist","m_x",100,0.9,1.5);
  tr->Project("mx_dist", "mx");
  mx_dist->Draw();
  auto peak = mx_dist->GetBinCenter(mx_dist->GetMaximumBin());
  auto dev = 0.03;
  mx_dist->Fit("gaus","","",peak-dev,peak+dev);
  c->SaveAs(file_name+".png");
}
