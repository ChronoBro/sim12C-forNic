{
  gROOT->Reset();
  TCanvas canvas("jDist","",500,800);
  canvas.Divide(1,2);
  canvas.cd(1);
  double xl[4]={0,1,2,3};
  double yl[4]={1.628,5.395,11.54,28.45};
  double xj[4]={.5,1.5,2.5,3.5};
  double yj[4]={2.849,5.843,11.515,26.80};


  TH2S frame("frame","",10,-.5,5.5,10,0,40);
  frame.GetXaxis()->SetTitle("m_{l}");
  frame.GetYaxis()->SetTitle("Relative Yield");
  frame.GetXaxis()->CenterTitle();
  frame.GetYaxis()->CenterTitle();
  frame.SetStats(kFALSE);
  frame.Draw();

  TGraph gl(4,xl,yl);
  gl.Draw("B");
  canvas.cd(2);

  TH2S frame2("frame2","",10,-.5,5.5,10,0,40);
  frame2.SetStats(kFALSE);
  frame2.GetXaxis()->SetTitle("m_{j}");
  frame2.GetYaxis()->SetTitle("Relative Yield");
  frame2.GetXaxis()->CenterTitle();
  frame2.GetYaxis()->CenterTitle();
  frame2.Draw();

  TGraph gj(4,xj,yj);
  gj.Draw("B");


}
