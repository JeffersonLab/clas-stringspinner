// ninja && clas-stringspinner --seed 29877 --trig 10000 --config beam_test --cut-inclusive 11,211,-211|grep M_X| awk '{print $2}' > output.dat
void draw_mx(TString file_name = "build/output.dat") {
  auto tr = new TTree();
  tr->ReadFile(file_name, "mx/D");
  new TCanvas();
  tr->Draw("mx");
  new TCanvas();
  auto mx_dist = new TH1D("mx_dist","m_x",200,0,2);
  tr->Project("mx_dist", "mx");
  mx_dist->Draw();
  mx_dist->Fit("gaus","","",0.8,0.95);
}
