{
  double x[4]={5.,9.23,10.,15.};
  double y[4]={2.44,1.778,1.815,2.99};

  TGraph g(4,x,y);
  g.Draw("A*");
  TF1 * funct = new TF1("funct","[0] + pow(x-[1],2)*[2]",0,20);
  funct->SetParameter(0,1.8);
  funct->SetParameter(1,10.);
  funct->SetParameter(2,.01);

  g.Fit(funct);

}
