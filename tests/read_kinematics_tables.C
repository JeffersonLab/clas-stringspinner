// read kinematics table into a TTree
void read_kinematics_tables(TString lund_name = "clas-stringspinner.dat")
{
  auto out_name = lund_name + ".root";
  auto out_file = new TFile(out_name, "RECREATE");
  auto trdis = new TTree("trdis", "trdis");
  auto tr1h  = new TTree("tr1h", "tr1h");
  auto tr2h  = new TTree("tr2h", "tr2h");
  trdis->ReadFile(lund_name + ".dis.table");
  tr1h->ReadFile(lund_name + ".1h.table");
  tr2h->ReadFile(lund_name + ".2h.table");
  trdis->Write();
  tr1h->Write();
  tr2h->Write();
  out_file->Close();
  auto read_file = new TFile(out_name, "READ");
  new TBrowser();
}
