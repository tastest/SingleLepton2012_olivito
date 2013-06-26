void scan_nsig(TString dir, TString output, bool do2d = true) {

  TChain* t = new TChain("Events");
  t->Add(dir + "/merged*.root");
  t->SetAlias("sparm_values","floats_sParmMaker_sparmvalues_CMS2.obj");

  TFile* fout = new TFile(output,"RECREATE");

  // 1D scan
  if (!do2d) {

    TH1F* h_masses = new TH1F("masses","", 21, 0., 525.);
    t->Draw("sparm_values[0]>>masses");

  // 2D scan
  } else {

    TH2F* h_masses = new TH2F("masses","", 21, 0., 525., 21, 0., 525.);
    t->Draw("sparm_values[1]:sparm_values[0]>>masses");

  }

  fout->Write();
  fout->Close();

}
