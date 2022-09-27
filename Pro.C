{
  gROOT->Reset();
  double x[8] = {-7./2.,-5./2.,-3./2.,-1./2.,1./2.,3./2.,5./2.,7./2.};
  double y[8] = {2./3.,8./21.,2./21.,0.,0.,2./21,8./21.,2./3.};


  TH2S frame ("frame","",10,-4,4,10,0,1.);
  frame.Draw();

  TGraph g(8,x,y);
  g.SetMarkerStyle(21);
  g.Draw("PL");
}
