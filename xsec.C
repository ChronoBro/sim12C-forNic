{
  gROOT->Reset();
  TCanvas canvas("xsec");
   canvas.SetLogy();

  ifstream feff("/home/Carbon8/nout/li6Beam/da/el/all/eff.dat");
  if (!feff.is_open()) 
    {
      cout << "eff file not found"<< endl;
     return;
    }
  double x[100];
  double y[100];
  int N = 0;

  double xx,yy;
  for (;;)
    {
      feff >> xx >> yy;
      if (feff.eof())break;
      if (feff.bad())break;

      x[N] = xx;
      y[N] = yy;
      N++;
    }
  feff.close();
  feff.clear();

  TH2S frame("frame","",10,0,20,10,.1,30.);
  frame.GetXaxis()->SetTitle("E^{*}(^{6}Li) [MeV]");
  frame.GetYaxis()->SetTitle("d#sigma/dE^{*} [mb/MeV]");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();

  frame.SetStats(kFALSE);
  frame.Draw();

  TFile f("/home/Carbon8/daq/corr_6Li.root");
  TH1I* hist = (TH1I*) f.Get("Li6/Erel_6Li_el");


  //TH1I* hist_pda = (TH1I*) f.Get("Be7/Erel_7Be_pda");
  //TH1I* hist_p6Li = (TH1I*) f.Get("Be7/Erel_7Be_p6Li");
  //TH1I* hist_pt3He = (TH1I*) f.Get("Be7/Erel_7Be_pt3He");



  TH1F* histN = new TH1F("histN","",1434,-3,30);
  //TH1F* histN_pda = new TH1F("histN_pda","",1000,0,60);
  //TH1F* histN_p6Li = new TH1F("histN_p6Li","",1000,0,60);
  //TH1F* histN_pt3He = new TH1F("histN_pt3He","",1000,0,60);


  float binWidth = 33./1434.;
  cout << "bin Width= " << binWidth << endl;


  

  TF1 *funct = new TF1("funct","[0]*pow((x-1.6),1)*exp(-(x-1.6)/[1]) + [2]*pow(x-1.6,2)*exp(-(x-1.6)/[3]) + [4]*pow(x-1.6,4)*exp(-(x-1.6)/[5])",1.6,20);



  //exponent == 1
  
  funct->SetParameter(0,1.31);
  funct->SetParameter(1,1.65);
  funct->SetParameter(2,.0);
  funct->SetParameter(3,4.4);
  funct->SetParameter(4,.000041);
  funct->SetParameter(5,3);
  
  //exponent = 0.5
  //funct->SetParameter(0,4.9);
  //funct->SetParameter(1,5.1);
  //funct->SetParameter(2,17);
  //funct->SetParameter(3,2.);
  //funct->Draw("same");
  double sum_3 = 0.;
  double sum_2 = 0.;
  double sum_b = 0.;
  funct->Draw("same");


  for (int i=1;i<1000;i++)
    {




      /*
      double tot = hist_pda->GetBinContent(i);
      tot /= 100.;

      histN_pda->SetBinContent(i,tot);



      tot = hist_p6Li->GetBinContent(i);
      tot /= 100.;

      histN_p6Li->SetBinContent(i,tot);


      tot = hist_pt3He->GetBinContent(i);
      tot /= 100.;

      histN_pt3He->SetBinContent(i,tot);
      */

      float e = hist->GetBinCenter(i);
      if (e-1.4743 < x[0]) continue;
      if (e > 20.) continue;
      float ET = e - 1.4743;
  



      int j = 0;
      for (;;)
	{

	  if (x[j] >=ET) break;
	  j++;
          if (j == N) break;
	}




      float eff = y[j-1] + (y[j]-y[j-1])/(x[j]-x[j-1])*(ET-x[j-1]);




      double tot = hist->GetBinContent(i);


      tot /= eff;

      tot /= binWidth;
      tot *= 7.629e-6;





      histN->SetBinContent(i,tot);

      if ( e > 3.3 && e < 5.3) sum_2 += tot - funct->Eval(e);
      else if ( e > 5.3 && e < 9.8) sum_3 += tot - funct->Eval(e);
      sum_b += funct->Eval(e);




    }

  histN->SetLineColor(4);
  histN->Draw("same");
  cout << histN->Integral()*binWidth << endl;


  return;
  histN_pda->SetLineColor(6);
  histN_pda->Draw("same");

  histN_p6Li->SetLineColor(3);
  histN_p6Li->Draw("same");

  histN_pt3He->SetLineColor(2);
  histN_pt3He->Draw("same");

  //  histNR->SetLineColor(3);
  //histNR->Draw("L same");

  TLatex tt; 
  tt.SetNDC();
  tt.DrawLatex(.4,.83,"^{3}He+#alpha breakup of ^{7}Be");


  cout << "sum in 7Be_2 = " << sum_2*binWidth << endl;
  cout << "sum in 7Be_3 = " << sum_3*binWidth << endl;
  cout << "sum back = " << sum_b*binWidth << endl;
  cout << " sum all = " << histN_a3He->Integral()*binWidth << endl;




  return;

  //***6Li line shape from Rmatrix fit
      ifstream lfile("star3/line.dat");




  double xl[300];
  double yl[300];
  int Nl = 0;
  for (;;)
    {
      lfile >> xx >> yy;
      if (lfile.eof()) break;
      if (lfile.bad()) break;
      xl[Nl] = xx;
      yl[Nl] = yy*1. + funct->Eval(xx);
      Nl++;

    }
  TGraph gl(Nl,xl,yl);
  gl.SetLineColor(6);
  gl.Draw("L");

  TLine line;
  line.SetLineStyle(2);
  line.DrawLine(4.57,0,4.57,10);
  line.DrawLine(6.73,0,6.73,10);
  line.DrawLine(7.21,0,7.21,10);
  line.DrawLine(9.27,0,9.27,10);
  line.DrawLine(9.9,0,9.9,10);


}
