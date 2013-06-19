void scan_nsig(TString dir, TString output, int xmin, int xmax, int xstep, int ymin = -1, int ymax = -1, int ystep = -1, int delta = -1) {

  TChain* t = new TChain("Events");
  t->Add(dir + "/merged*.root");
  t->SetAlias("sparm_values","floats_sParmMaker_sparmvalues_CMS2.obj");

  TFile* fout = new TFile(output,"RECREATE");

  int xbins = (xmax - xmin)/xstep + 1;
  // 1D scan: only loop over x values
  if (ymin < 0) {

    TH1F* h_masses = new TH1F("masses","", 21, 0., 525.);

    for (int ix = 0; ix < xbins; ++ix) {
      int xval = ix * xstep + xmin;
      int n = t->GetEntries(Form("sparm_values[0] == %d",xval));
      h_masses->SetBinContent(ix + xmin/xstep + 1, n);
    }

  // 2D scan: loop over x and y values
  } else {

    TH2F* h_masses = new TH2F("masses","", 21, 0., 525., 21, 0., 525.);
    int ybins = (ymax - ymin)/ystep + 1;

    for (int ix = 0; ix < xbins; ++ix) {
      int xval = ix * xstep + xmin;
      for (int iy = 0; iy < ybins; ++iy) {
	int yval = iy * ystep + ymin;
	if (xval - yval < delta) continue;
	int n = t->GetEntries(Form("sparm_values[0] == %d && sparm_values[1] == %d",xval,yval));
        h_masses->SetBinContent(ix + xmin/xstep + 1, iy + ymin/ystep + 1, n);
      }
    }

  }

  fout->Write();
  fout->Close();

}
