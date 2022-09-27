{

  gROOT->Reset();
  TCanvas canvas("2d","",600,600);
  canvas.Divide(2,2);


  TH2S frame("frame","",10,-1,1,10,0,360);
  frame.GetXaxis()->SetTitle("cos(#psi)");
  frame.GetYaxis()->SetTitle("#chi [deg]");
  frame.SetStats(kFALSE);


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
        corr->SetBinContent(i,j,yexp/effic);
	}
    }


  canvas.cd(2);
  frame.Draw();
  eff->Smooth();
  eff->SetMinimum(0.);
  eff->Draw(" same zcol");


  const int N = 28;
  sle SLE(N);
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
          ycom[1] = y31;
          ycom[2] = y32;
          ycom[3] = y33;
          ycom[4] = y31*cos(2.*chi);   // m = 1 and -1
          ycom[5] = y31*sin(2.*chi);   // m = 1 and -1
          ycom[6] = y32*cos(4.*chi);   // m = 2 and -2
          ycom[7] = y32*sin(4.*chi);   // m = 2 and -2
          ycom[8] = y33*cos(6.*chi);   // m = 3 and -3
          ycom[9] = y33*sin(6.*chi);   // m = 3 and -3
	  ycom[10] = -p33*p32*cos(chi); // m=3 and m=2 also M=-3 and m=-2
          ycom[11] = -p33*p32*sin(chi);
	  ycom[12] = -p32*p31*cos(chi);// m=2 and m=1 also M=-2 and m=-1
          ycom[13] = -p32*p31*sin(chi);
	  ycom[14] = -p31*p30*cos(chi);// m=1 and m=0 also M=-1 and m=0
          ycom[15] = -p31*p30*sin(chi);
	  ycom[16] = p33*p31*cos(2.*chi); // m=3 and m=1 also M=-3 and m=-1
          ycom[17] = p33*p31*sin(2.*chi);
          ycom[18] = p32*p30*cos(2.*chi); //m=2 and m=0 also m=-2 and m=0
          ycom[19] = p32*p30*sin(2.*chi); 
	  ycom[20] = -p33*p30*cos(3.*chi); //m=3 and m=0 also m=-3 and m=0
	  ycom[21] = -p33*p30*sin(3.*chi); 
          ycom[22] = -p32*p31*cos(3.*chi); //m=2 and m=-1 also m-2 and m=1
          ycom[23] = -p32*p31*sin(3.*chi); 
          ycom[24] = p33*p31*cos(4.*chi); //m=3 and m=-1 also m=-3 and m=1
          ycom[25] = p33*p31*sin(4.*chi); 
          ycom[26] = -p33*p32*cos(5.*chi); //m=3 and m=-2 also m=-3 and m=2
          ycom[27] = -p33*p32*sin(5.*chi); 


	  double yexp = corr->GetBinContent(i,j);

	  for (int n=0;n<N;n++)
	    {
              SLE.Y[n] += ycom[n]*yexp;
	      for (int m =0;m<N;m++)SLE.M[n][m] += ycom[n]*ycom[m];
	    }

	}
    }

  SLE.solve();


  for (int i=0;i<N;i++) cout << i << " " << SLE.Y[i]<< endl;

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
          ycom[1] = y31;
          ycom[2] = y32;
          ycom[3] = y33;
          ycom[4] = y31*cos(2.*chi);   // m = 1 and -1
          ycom[5] = y31*sin(2.*chi);   // m = 1 and -1
          ycom[6] = y32*cos(4.*chi);   // m = 2 and -2
          ycom[7] = y32*sin(4.*chi);   // m = 2 and -2
          ycom[8] = y33*cos(6.*chi);   // m = 3 and -3
          ycom[9] = y33*sin(6.*chi);   // m = 3 and -3
	  ycom[10] = -p33*p32*cos(chi); // m=3 and m=2 also M=-3 and m=-2
          ycom[11] = -p33*p32*sin(chi);
	  ycom[12] = -p32*p31*cos(chi);// m=2 and m=1 also M=-2 and m=-1
          ycom[13] = -p32*p31*sin(chi);
	  ycom[14] = -p31*p30*cos(chi);// m=1 and m=0 also M=-1 and m=0
          ycom[15] = -p31*p30*sin(chi);
	  ycom[16] = p33*p31*cos(2.*chi); // m=3 and m=1 also M=-3 and m=-1
          ycom[17] = p33*p31*sin(2.*chi);
          ycom[18] = p32*p30*cos(2.*chi); //m=2 and m=0 also m=-2 and m=0
          ycom[19] = p32*p30*sin(2.*chi); 
	  ycom[20] = -p33*p30*cos(3.*chi); //m=3 and m=0 also m=-3 and m=0
	  ycom[21] = -p33*p30*sin(3.*chi); 
          ycom[22] = -p32*p31*cos(3.*chi); //m=2 and m=-1 also m-2 and m=1
          ycom[23] = -p32*p31*sin(3.*chi); 
          ycom[24] = p33*p31*cos(4.*chi); //m=3 and m=-1 also m=-3 and m=1
          ycom[25] = p33*p31*sin(4.*chi); 
          ycom[26] = -p33*p32*cos(5.*chi); //m=3 and m=-2 also m=-3 and m=2
          ycom[27] = -p33*p32*sin(5.*chi); 

	  double yfit = 0.;
	  for (int n=0;n<N;n++) yfit += ycom[n]*SLE.Y[n];
	  fit->SetBinContent(i,j,yfit);

	  double yexp = corr->GetBinContent(i,j);
	  chisq += pow(yexp-yfit,2);

	}
    }

  cout << "chisq = " << chisq << endl;
  canvas.cd(4);
  frame.Draw();
  fit->SetMinimum(0.);
  fit->Draw("ZCOL SAME");

  canvas.cd(3);
  frame.Draw();
  corr->Smooth();
  corr->SetMinimum(0.);
  corr->SetMaximum(110.);
  corr->Draw("same zcol");


  canvas.cd(1);
  frame.Draw();
    histX->Smooth();
    histX->SetMinimum(0.);
    histX->SetMaximum(100.);
  histX->Draw("ZCOL same");

  string tit[12];
  tit[0] = "1 & -1";
  tit[1] = "2 & -2";
  tit[2] = "3 & -3";
  tit[3] = "3 & 2 and -3 & -2";
  tit[4] = "2 & 1 and -2 & -1";
  tit[5] = "1 & 0 and -1 & 0";
  tit[6] = "3 & 1 and -3 & -1";
  tit[7] = "2 & 0 and -2 & 0";
  tit[8] = "3 & 0 and -3 & 0";
  tit[9] = "2 & -1 and -2 & 1";
  tit[10] = "3 & -1 and -3 & 1";
  tit[11] = "3 & -2 and -3 & 2";
  for (int i=0;i<12;i++)
    {
      double yi = sqrt(pow(SLE.Y[4+2*i],2)+pow(SLE.Y[4+2*i+1],2));
      double phase = atan2(SLE.Y[4+2*i+1],SLE.Y[4+2*i])*180./acos(-1.);
      cout << i << " " << yi << " " << phase << " " << tit[i] << endl;
    }
}
