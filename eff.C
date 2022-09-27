{
  gROOT->Reset();
  TH2S frame("frame","",10,0,20,10,0,.25);
  frame.GetXaxis()->SetTitle("E_{T} [MeV]");
  frame.GetYaxis()->SetTitle("efficiency");
  frame.SetStats(kFALSE);
  frame.SetTitle("^{7}Be#rightarrow p+#alpha");
  frame.Draw();
  const int N = 28;
  double x[N]={.01,.1,.2,.3,.4,.5,.6,.7,.8,1.,1.25,1.5,1.75,2.,2.5,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.};
  double y[N]={.079,.1985,.2065,.213,.225,.234,.237,.234,.228,.217,.207,.201,.1975,.193,.177,.160,.1298,.105,.0858,.0704,.0603,.0526,.0463,.0378,.0318,.028,.0251,.0225};

  TGraph g(N,x,y);
  g.Draw(" *L");
  //g.Fit("pol2");

}
