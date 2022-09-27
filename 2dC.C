{

  gROOT->Reset();
TStyle * Sty = new TStyle("MyStyle","MyStyle");
Sty->SetOptTitle(0);
Sty->SetOptStat(0);
//Sty->SetPalette(8,0);
Sty->SetCanvasColor(10);
Sty->SetCanvasBorderMode(0);
Sty->SetFrameLineWidth(3);
Sty->SetFrameFillColor(10);
Sty->SetPadColor(10);
Sty->SetPadTickX(1);
Sty->SetPadTickY(1);
Sty->SetPadBottomMargin(.16);
Sty->SetPadTopMargin(.11);
Sty->SetPadLeftMargin(.19);
Sty->SetPadRightMargin(.15);
Sty->SetHistLineWidth(3);
Sty->SetHistLineColor(kRed);
Sty->SetFuncWidth(3);
Sty->SetFuncColor(kGreen);
Sty->SetLineWidth(3);
Sty->SetLabelSize(0.06,"xyz");
Sty->SetLabelOffset(0.02,"y");
Sty->SetLabelOffset(0.02,"x");
Sty->SetLabelColor(kBlack,"xyz");
Sty->SetTitleSize(0.06,"xyz");
Sty->SetTitleOffset(-1.7,"y");
Sty->SetTitleOffset(1.2,"x");
Sty->SetTitleFillColor(10);
Sty->SetTitleTextColor(kBlack);
Sty->SetTickLength(.03,"xz");
Sty->SetTickLength(.02,"y");
Sty->SetNdivisions(5,"xyz");
Sty->SetEndErrorSize(0);
gROOT->SetStyle("MyStyle");
gROOT->ForceStyle();



  TCanvas canvas("2dC7Be72","",800,800);

  TLatex tt;
  TLine ll;
  ll.SetLineWidth(3);  

double overlapy = .06;
TPad *pad1 = new TPad("pad1","",0.,0.,.5+overlapy,1.);
TPad *pad2 = new TPad("pad2","",.5-overlapy/2.,0.,1.,1.);
pad1->SetFillStyle(4000);
pad2->SetFillStyle(4000);

pad1->Draw();
pad2->Draw();

pad1->cd();
double overlap = .05;
TPad *pad11 = new TPad("pad1","",0.,0.,1.,.5+overlap);
TPad *pad21 = new TPad("pad2","",0.,.5-overlap/2.,1.,1.);
pad11->SetFillStyle(4000);
pad21->SetFillStyle(4000);

pad11->Draw();
pad21->Draw();
pad2->cd();
TPad *pad12 = new TPad("pad1","",0.,0.,1.,.5+overlap);
TPad *pad22 = new TPad("pad2","",0.,.5-overlap/2.,1.,1.);
pad12->SetFillStyle(4000);
pad22->SetFillStyle(4000);

pad12->Draw();
pad22->Draw();

  TH2S frame1("frame1","",10,-1,1,10,0,360);
  //frame1.GetXaxis()->SetTitle("\\cos(\\psi)");
  frame1.GetYaxis()->SetTitle("\\chi \\,\\, [deg]");

  frame1.GetXaxis()->CenterTitle();
  frame1.GetXaxis()->SetLabelSize(0.);
  frame1.GetYaxis()->CenterTitle();
  frame1.SetStats(kFALSE);
  frame1.GetXaxis()->SetTicks("+-");
  frame1.GetYaxis()->SetTicks("+-");
  frame1.GetYaxis()->SetNdivisions(4,kFALSE);

  TFile file("/home/Carbon8/daq/corr_7Be.root");
  TH2I* histX = (TH2I*) file.Get("Be7/theta3DeltaPhi_7Be_2_el");



  float Nx = histX->Integral();
  cout << Nx << endl;

  TFile fileS("sim.root");
  TH2I* sr = (TH2I*) fileS.Get("theta3DeltaPhi_R");
  double Nsr = sr->Integral();
  TH2I* sp = (TH2I*) fileS.Get("theta3DeltaPhi_P");
  double Nsp = sp->Integral();



  TH2F* eff = new TH2F("eff","",50,-1,1,90,0,360);
  TH2F* corr = new TH2F("corr","",50,-1,1,90,0,360);
  TH2F* weight = new TH2F("weight","",50,-1,1,90,0,360);
  TH2F* fit = new TH2F("fit","",50,-1,1,90,0,360);
  eff->SetStats(kFALSE);
  for (int i=1;i<=50;i++)
    {
      for (int j=1;j<=90;j++)
	{
	double yr  = sr->GetBinContent(i,j);
	double yp  = sp->GetBinContent(i,j);
	double effic = yr/yp/Nsr*Nsp;


	eff->SetBinContent(i,j,effic*100.);
        double yexp = histX->GetBinContent(i,j);        
	if (effic == 0. || yexp == 0.)
  	  {
           corr->SetBinContent(i,j,0.);
           weight->SetBinContent(i,j,0.);
	  }

        else 
	  {
           corr->SetBinContent(i,j,yexp/effic);
	   weight->SetBinContent(i,j,pow(effic,2)/yexp);
	  }



	}
    }


  pad22->cd();
  TH2S frame2("frame2","",10,-1,1,10,0,360);
  //frame2.GetXaxis()->SetTitle("\\cos(\\psi)");
  //frame2.GetYaxis()->SetTitle("\\chi \\,\\,[deg]");
  frame2.GetYaxis()->SetNdivisions(4,kFALSE);
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();

  frame2.GetXaxis()->SetLabelSize(0.);
  frame2.GetYaxis()->SetLabelSize(0.);
  frame2.SetStats(kFALSE);
  frame2.GetXaxis()->SetTicks("+-");
  frame2.GetYaxis()->SetTicks("+-");
  frame2.Draw();

  eff->Smooth();
  eff->SetMinimum(0.);
  eff->Draw(" same zcol");

  ll.DrawLine(-1,90,-.92,90);
  ll.DrawLine(.92,90,1.,90);
  ll.DrawLine(-1,180,-.92,180);
  ll.DrawLine(.92,180,1.,180);
  ll.DrawLine(-1,270,-.92,270);
  ll.DrawLine(.92,270,1.,270);
  ll.DrawLine(-.5,0,-.5,10);
  ll.DrawLine(-.5,350,-.5,360);
  ll.DrawLine(0,0,0,10);
  ll.DrawLine(0,350,0,360);
  ll.DrawLine(.5,0,.5,10);
  ll.DrawLine(.5,350,.5,360);
  ll.DrawLine(-1,0,1,0);
  ll.DrawLine(1.,0.,1.,360.);
  ll.DrawLine(1.,360.,-1,360.);
  ll.DrawLine(-1,360.,-1,0.);
  const int N = 16;
  sle SLE(N);
  inverse Inverse(N);
  SLE.clear();
  
  //fit

  double ycom[40];
  for (int i=1;i<=50;i++)
    {
      for (int j=1;j<=90;j++)
	{
	  double cosx = histX->GetXaxis()->GetBinCenter(i);


          double y33 = 225.*pow(1.-pow(cosx,2),3)*7./720.;
          double y32 = 225.*pow(cosx,2)*pow(pow(cosx,2)-1.,2)*7./120.;
          double y31 = (1.-pow(cosx,2))*pow(5.*pow(cosx,2)-1.,2)*9./4.*7./12.;
          double y30 = pow(5.*pow(cosx,3)-3.*cosx,2)/4.*7.;

	  double p33 = -15.*pow(1.-pow(cosx,2),3./2.)*sqrt(7./720.);
          double p32 = -15.*cosx*(pow(cosx,2)-1.)*sqrt(7./120.);
          double p31 = -3./2.*sqrt(1.-pow(cosx,2))*(5.*pow(cosx,2)-1.)*sqrt(7./12.);
	  double p30 = 1./2.*(5.*pow(cosx,3)-3.*cosx)*sqrt(7.);

	  double chi = histX->GetYaxis()->GetBinCenter(j);
          chi *= acos(-1.)/180.;  //radians



	  ycom[0] = y30;
          ycom[1] = 2.*y31;
          ycom[2] = 2.*y32;
          ycom[3] = 2.*y33;
          ycom[4] = -2.*y31*cos(2.*chi);   // m = 1 and -1
          ycom[5] = -2.*y32*cos(4.*chi);   // m = 2 and -2
          ycom[6] = -2.*y33*cos(6.*chi);   // m = 3 and -3
	  ycom[7] = -4.*p33*p32*cos(chi); // m=3 and m=2 also M=-3 and m=-2
	  ycom[8] = -4.*p32*p31*cos(chi);// m=2 and m=1 also M=-2 and m=-1
	  ycom[9] = -4.*p31*p30*cos(chi);// m=1 and m=0 also M=-1 and m=0
	  ycom[10] = 4.*p33*p31*cos(2.*chi); // m=3 and m=1 also M=-3 and m=-1
          ycom[11] = 4.*p32*p30*cos(2.*chi); //m=2 and m=0 also m=-2 and m=0
	  ycom[12] = -4.*p33*p30*cos(3.*chi); //m=3 and m=0 also m=-3 and m=0
          ycom[13] = 4.*p32*p31*cos(3.*chi); //m=2 and m=-1 also m-2 and m=1
          ycom[14] = -4.*p33*p31*cos(4.*chi); //m=3 and m=-1 also m=-3 and m=1
          ycom[15] = 4.*p33*p32*cos(5.*chi); //m=3 and m=-2 also m=-3 and m=2


	  double yexp = corr->GetBinContent(i,j);
	  double weigh = weight->GetBinContent(i,j);
	  if (yexp == 0.) continue;
	  for (int n=0;n<N;n++)
	    {
              SLE.Y[n] += ycom[n]*yexp*weigh;
	      for (int m =0;m<N;m++)SLE.M[n][m] += ycom[n]*ycom[m]*weigh;
	    }

	}
    }

  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) Inverse.A[i][j] = SLE.M[i][j];

  SLE.solve();
  Inverse.solve();
  //check A*A-1 is the identity
  Inverse.S3000();

  for (int i=0;i<N;i++) cout << i << " " << SLE.Y[i]<< " " << sqrt(Inverse.AI[i][i]) << endl;

  double chisq = 0.;
  for (int i=1;i<=50;i++)
    {
      for (int j=1;j<=90;j++)
	{
	  double cosx = histX->GetXaxis()->GetBinCenter(i);


          double y33 = 225.*pow(1.-pow(cosx,2),3)*7./720.;
          double y32 = 225.*pow(cosx,2)*pow(pow(cosx,2)-1.,2)*7./120.;
          double y31 = (1.-pow(cosx,2))*pow(5.*pow(cosx,2)-1.,2)*9./4.*7./12.;
          double y30 = pow(5.*pow(cosx,3)-3.*cosx,2)/4.*7.;


	  double p33 = -15.*pow(1.-pow(cosx,2),3./2.)*sqrt(7./720.);
          double p32 = -15.*cosx*(pow(cosx,2)-1.)*sqrt(7./120.);
          double p31 = -3./2.*sqrt(1.-pow(cosx,2))*(5.*pow(cosx,2)-1.)*sqrt(7./12.);
	  double p30 = 1./2.*(5.*pow(cosx,3)-3.*cosx)*sqrt(7.);


	  double chi = histX->GetYaxis()->GetBinCenter(j);
          chi *= acos(-1.)/180.;  //radians

          
	  ycom[0] = y30;
          ycom[1] = 2.*y31;
          ycom[2] = 2.*y32;
          ycom[3] = 2.*y33;
          ycom[4] = -2.*y31*cos(2.*chi);   // m = 1 and -1
          ycom[5] = -2.*y32*cos(4.*chi);   // m = 2 and -2
          ycom[6] = -2.*y33*cos(6.*chi);   // m = 3 and -3
	  ycom[7] = -4.*p33*p32*cos(chi); // m=3 and m=2 also M=-3 and m=-2
	  ycom[8] = -4.*p32*p31*cos(chi);// m=2 and m=1 also M=-2 and m=-1
	  ycom[9] = -4.*p31*p30*cos(chi);// m=1 and m=0 also M=-1 and m=0
	  ycom[10] = 4.*p33*p31*cos(2.*chi); // m=3 and m=1 also M=-3 and m=-1
          ycom[11] = 4.*p32*p30*cos(2.*chi); //m=2 and m=0 also m=-2 and m=0
	  ycom[12] = -4.*p33*p30*cos(3.*chi); //m=3 and m=0 also m=-3 and m=0
          ycom[13] = 4.*p32*p31*cos(3.*chi); //m=2 and m=-1 also m-2 and m=1
          ycom[14] = -4.*p33*p31*cos(4.*chi); //m=3 and m=-1 also m=-3 and m=1
          ycom[15] = 4.*p33*p32*cos(5.*chi); //m=3 and m=-2 also m=-3 and m=2


	  double yfit = 0.;
	  for (int n=0;n<N;n++) yfit += ycom[n]*SLE.Y[n];
	  fit->SetBinContent(i,j,yfit);

	  double yexp = corr->GetBinContent(i,j);
	  double weigh = weight->GetBinContent(i,j);
	  chisq += pow(yexp-yfit,2)*weigh;


	}
    }
  double chisqnu = chisq/(50.*90.-(double)N);
  cout << "chisq = " << chisqnu << endl;
 
  pad12->cd();
  TH2S frame4("frame4","",10,-1,1,10,0,360);
  frame4.GetXaxis()->SetTitle("\\cos(\\psi)");
  //frame4.GetYaxis()->SetTitle("\\chi \\,\\,[deg]");
  frame4.GetXaxis()->CenterTitle();
  frame4.GetYaxis()->CenterTitle();
  frame4.GetYaxis()->SetLabelSize(0.);
  frame4.SetStats(kFALSE);
  frame4.GetYaxis()->SetNdivisions(4,kFALSE);
  frame4.GetXaxis()->SetTicks("+-");
  frame4.GetYaxis()->SetTicks("+-");
  frame4.Draw();

  fit->SetMinimum(0.);
  fit->Draw("ZCOL SAME");

  ll.DrawLine(-1,90,-.92,90);
  ll.DrawLine(.92,90,1.,90);
  ll.DrawLine(-1,180,-.92,180);
  ll.DrawLine(.92,180,1.,180);
  ll.DrawLine(-1,270,-.92,270);
  ll.DrawLine(.92,270,1.,270);
  ll.DrawLine(-.5,0,-.5,10);
  ll.DrawLine(-.5,350,-.5,360);
  ll.DrawLine(0,0,0,10);
  ll.DrawLine(0,350,0,360);
  ll.DrawLine(.5,0,.5,10);
  ll.DrawLine(.5,350,.5,360);
    ll.DrawLine(-1,0,1,0);
  ll.DrawLine(1.,0.,1.,360.);
  ll.DrawLine(1.,360.,-1,360.);
  ll.DrawLine(-1,360.,-1,0.);

  pad11->cd();
  TH2S frame3("frame3","",10,-1,1,10,0,360);
  frame3.GetXaxis()->SetTitle("\\cos(\\psi)");
  frame3.GetYaxis()->SetTitle("\\chi\\,\\, [deg]");
  frame3.GetXaxis()->CenterTitle();
  frame3.GetYaxis()->CenterTitle();
  frame3.SetStats(kFALSE);
  frame3.GetYaxis()->SetNdivisions(4,kFALSE);
  frame3.GetXaxis()->SetTicks("+-");
  frame3.GetYaxis()->SetTicks("+-");
  frame3.Draw();


  corr->Smooth();
  corr->SetMinimum(0.);
  corr->SetMaximum(110.);
  corr->Draw("same zcol");
 
  ll.DrawLine(-1,90,-.92,90);
  ll.DrawLine(.92,90,1.,90);
  ll.DrawLine(-1,180,-.92,180);
  ll.DrawLine(.92,180,1.,180);
  ll.DrawLine(-1,270,-.92,270);
  ll.DrawLine(.92,270,1.,270);
  ll.DrawLine(-.5,0,-.5,10);
  ll.DrawLine(-.5,350,-.5,360);
  ll.DrawLine(0,0,0,10);
  ll.DrawLine(0,350,0,360);
  ll.DrawLine(.5,0,.5,10);
  ll.DrawLine(.5,350,.5,360);
  ll.DrawLine(-1,0,1,0);
  ll.DrawLine(1.,0.,1.,360.);
  ll.DrawLine(1.,360.,-1,360.);
  ll.DrawLine(-1,360.,-1,0.);

  pad21->cd();
  frame1.GetYaxis()->SetNdivisions(4,kFALSE);
  frame1.Draw();
    histX->Smooth();
    histX->SetMinimum(0.);
    histX->SetMaximum(100.);
  histX->Draw("ZCOL same");

  ll.DrawLine(-1,90,-.92,90);
  ll.DrawLine(.92,90,1.,90);
  ll.DrawLine(-1,180,-.92,180);
  ll.DrawLine(.92,180,1.,180);
  ll.DrawLine(-1,270,-.92,270);
  ll.DrawLine(.92,270,1.,270);
  ll.DrawLine(-.5,0,-.5,10);
  ll.DrawLine(-.5,350,-.5,360);
  ll.DrawLine(0,0,0,10);
  ll.DrawLine(0,350,0,360);
  ll.DrawLine(.5,0,.5,10);
  ll.DrawLine(.5,350,.5,360);
  ll.DrawLine(-1,0,1,0);
  ll.DrawLine(1.,0.,1.,360.);
  ll.DrawLine(1.,360.,-1,360.);
  ll.DrawLine(-1,360.,-1,0.);

  float tot = SLE.Y[0] + 2.*SLE.Y[1] + 2.*SLE.Y[2] + 2.*SLE.Y[3];

  string tit[16];
  tit[0] = "0";
  tit[1] = "1";
  tit[2] = "2";
  tit[3] = "3";
  tit[4] = "1 & -1";
  tit[5] = "2 & -2";
  tit[6] = "3 & -3";
  tit[7] = "3 & 2 and -3 & -2";
  tit[8] = "2 & 1 and -2 & -1";
  tit[9] = "1 & 0 and -1 & 0";
  tit[10] = "3 & 1 and -3 & -1";
  tit[11] = "2 & 0 and -2 & 0";
  tit[12] = "3 & 0 and -3 & 0";
  tit[13] = "2 & -1 and -2 & 1";
  tit[14] = "3 & -1 and -3 & 1";
  tit[15] = "3 & -2 and -3 & 2";
  for (int i=0;i<16;i++)
    {

      cout << i << " " << SLE.Y[i]/tot << " " << sqrt(Inverse.AI[i][i]*chisqnu)/tot << " " << tit[i] << endl;
    }

  canvas.cd();
  tt.SetNDC();
  tt.SetTextSize(.04);
  tt.DrawLatex(.2,.95,"(a) measured");
  tt.DrawLatex(.65,.95,"(b) efficiency");
  tt.DrawLatex(0.2,0.5,"(c) corrected");
  tt.DrawLatex(.72,0.5,"(d) fit");




}
