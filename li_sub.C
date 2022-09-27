{

  gROOT->Reset();
  TCanvas canvas("li_sub");
  //canvas.SetLogy();

  TH2S frame("frame","",10,1,10,10,.1,50000);
  frame.GetXaxis()->SetTitle("E^{*} [MeV]");
  frame.GetYaxis()->SetTitle("Counts");
  frame.SetStats(kFALSE);
  frame.SetTitle("6Li contamination");

  frame.Draw();

  TFile f7Be("/home/Carbon8/daq/corr_7Be.root");
  TH1I* hpa_7Be = (TH1I*) f7Be.Get("Be7/Erel_7Be_a3He");
  hpa_7Be->Draw("same");



  TFile f6Li("/home/Carbon8/daq/corr_6Li.root");
  TH1I* hpa_6Li = (TH1I*) f6Li.Get("Be7/Erel_7Be_a3He");


  hpa_6Li->Scale(1./7.979);
  hpa_6Li->SetLineColor(2); 
  hpa_6Li->Draw("same");

  TF1 * funct1 = new TF1("funct1","(x-[0])*[1] + pow(x-[0],3)*[2]",1.5,5.5);
  funct1->SetParameter(0,2.1);
  funct1->SetParameter(1,2300);
  funct1->SetParameter(2,-80.0);
  funct1->SetLineColor(2);
  funct1->Draw("same L");

  float sum = 0.;
  for (int i=1;i<1000;i++)
    {
      float x = hpa_7Be->GetBinCenter(i);
      if (x < 3.5) continue;
      if (x > 5.5) break;
      float y = hpa_7Be->GetBinContent(i);
      sum += y - funct1->Eval(x);
    }

  cout << sum << endl;

}
