// read kinematics table into a TTree
void read_table_file(TString table_file_name = "clas-stringspinner.dat.2h.table")
{
  auto tr = new TTree("tr", "kinematics tree");
  tr->ReadFile(table_file_name);
}
