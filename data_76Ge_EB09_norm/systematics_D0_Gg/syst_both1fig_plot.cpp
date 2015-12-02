// Ge neutron-resonance parameter systematics, ACL
{	
	// Cleaning up (necessary when not running script in root directory)	
	gROOT->Reset();
	gROOT->SetStyle("Plain");	// no fill in pads
	gStyle->SetOptTitle(0);		// no histogram title
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);		// no statistics
	
	// Deleting histograms from last run to avoid memory leak	
	m = (TH2F*)gROOT->FindObject("h1"); 
	if (m){ m->Delete();}
    m = (TH2F*)gROOT->FindObject("h2");
    if (m){ m->Delete();}
	

    // Correction factor; 
    double factor = 1.2;

	// Values for D0 taken from RIPL-3
    double Bnexp[5]     = {6.072, 6.505, 6.783, 7.416, 10.196};
    double massexp[5]   = {77.00, 75.00, 73.00, 71.00, 74.00};
    double rhoexp[5]    = {5.926E+03, 8.951E+03, 1.794E+04, 2.317E+04, 1.115E+05};
    double rhoerrexp[5] = {1.745E+03, 3.446E+03, 4.979E+03, 6.363E+03,  2.710E+04};
    double gammaexp[5]   = {115, 195, 162, 185, 195};
    double gammaexpup[5]   = {115+35, 195+60, 162+50, 185+50, 195+45};
    double gammaexplow[5]   = {115-35, 195-60, 162-50, 185-50, 195-45};
    double gammaerrexp[5]= {35, 60, 50, 50, 45};
    double Bn_err[5] = {0.};

    double Bnrobinodd[6]     = {5.698, 6.072, 6.505, 6.783, 7.416, 8.192};
    double massrobinodd[6]   = {79.00, 77.00, 75.00, 73.00, 71.00, 69.00};
    double robinodd[6]  = {factor*0.26123E+04, factor*0.54871E+04, factor*0.10215E+05, factor*0.14343E+05, factor*0.16949E+05, factor*0.17585E+05};
    
    double Bnrobineven[5]     = {8.719, 9.428, 10.196, 10.749, 11.534};
    double massrobineven[5]   = {78.00, 76.00, 74.000, 72.000, 70.000};
    double robineven[5]  = {factor*0.17117E+05, factor*0.39193E+05, factor*0.74251E+05, factor*0.94626E+05, factor*0.92454E+05};
    
    // Robin values for unknown 76Ge    
    double Bn_76Ge[1] = {9.428};
    double robin76Ge[1] = {factor*0.39193E+05};
    double robin76Geerr[1] = {0.*robin76Ge[0]};
    double robin76Geerrlow[1] = {0.*robin76Ge[0]};
   
    cout << " Estimated lower level density, 76Ge: " 
        << robin76Ge[0] << " 1/MeV. " << endl;

 	// Making graphs			
	TGraphErrors *d2rho = new TGraphErrors(5,Bnexp,rhoexp,Bn_err,rhoerrexp);
	
    TGraph *robin_graph_odd = new TGraph(6,Bnrobinodd,robinodd);
    TGraph *robin_graph_even = new TGraph(5,Bnrobineven,robineven);
	TGraphAsymmErrors *robin76Ge_graph = new TGraphAsymmErrors(1,Bn_76Ge,robin76Ge,Bn_err,Bn_err,robin76Geerrlow,robin76Geerr);
	
	// Making a canvas
    TCanvas *c1 = new TCanvas("c1","Systematics NLD and <Gg>, 76Ge",600,750);
    c1.Divide(1,2,0,0);
    
    // Making histograms for the axis
	TH2F *h1 = new TH2F("h1","Systematics, Cd, D0",100,4.710,11.94,50,1.7E03,3.7E05); //0,2.7E06
    
	
    c1.cd();
    c1.cd(1);
    c1_1.SetLogy();
	c1_1.SetLeftMargin(0.14);
	c1_1.SetBottomMargin(0.05);
	c1_1.SetRightMargin(0.05);
    
	h1->GetXaxis()->SetTitle(" Neutron separation energy S_{n} (MeV)");
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->SetTitleSize(0.055);
	h1->GetXaxis()->CenterTitle();
	h1->GetXaxis()->SetLabelSize(0.06);
	h1->GetXaxis()->SetLabelOffset(10.05);
	h1->GetXaxis()->SetTitleOffset(1.2);
    h1->GetXaxis()->SetLabelFont(42);
	h1->GetYaxis()->SetTitle(" #rho(S_{n}) (MeV^{-1})");
    h1->GetYaxis()->SetTitleFont(42);
	h1->GetYaxis()->CenterTitle();
	h1->GetYaxis()->SetLabelSize(0.06);
	h1->GetYaxis()->SetTitleSize(0.07);
    h1->GetYaxis()->SetLabelFont(42);
	h1->GetYaxis()->SetTitleOffset(0.9);
	h1->Draw();

    d2rho->SetMarkerStyle(21);
    d2rho->SetMarkerSize(0.7);
    d2rho->Draw("p");
    
    robin_graph_odd->SetMarkerStyle(25);
    robin_graph_odd->SetMarkerColor(kRed);
    robin_graph_odd->SetLineColor(kRed);
    robin_graph_odd->Draw("lp");
    
    robin_graph_even->SetMarkerStyle(25);
    robin_graph_even->SetMarkerColor(kBlue);
    robin_graph_even->SetLineColor(kBlue);
    robin_graph_even->Draw("lp");

    robin76Ge_graph->SetMarkerColor(kBlue);
    robin76Ge_graph->SetLineColor(kBlue);
    robin76Ge_graph->SetMarkerStyle(33);
    robin76Ge_graph->SetMarkerSize(1.8);
    robin76Ge_graph->Draw("p");
    

	TLatex t;
	t.SetTextSize(0.065);
    t.SetTextFont(42);
	t.DrawLatex(5.33,260500,"(a)");
	
    t.SetTextSize(0.06);
	t.DrawLatex(5.17,3605,"^{79}Ge");
	t.DrawLatex(5.6,7975,"^{77}Ge");
	t.DrawLatex(6.,13975,"^{75}Ge");
	t.DrawLatex(6.3,23975,"^{73}Ge");
 	t.DrawLatex(6.95,32000,"^{71}Ge");
 	t.DrawLatex(7.75,24000,"^{69}Ge");
  	t.DrawLatex(8.3,26000,"^{78}Ge");
    t.SetTextColor(kBlue);
  	t.DrawLatex(8.95,76000,"^{76}Ge");
    t.SetTextColor(kBlack);
  	t.DrawLatex(9.7,150000,"^{74}Ge");
  	t.DrawLatex(10.4,140000,"^{72}Ge");
  	t.DrawLatex(11.2,130000,"^{70}Ge");


    TLegend *leg = new TLegend(0.18,0.66,0.55,0.92);
    //leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.052);
    leg->AddEntry(d2rho," From n.res.data ","P");
    leg->AddEntry(robin_graph_odd," Systematics odd ","LP");	//x 0.6
    leg->AddEntry(robin_graph_even," Systematics even  ","LP");	//x 0.6
    leg->AddEntry(robin76Ge_graph," ^{76}Ge, estimated ","P");
    leg->Draw();

    
    /*** Next is for <Gg> ***/
    TH2F *h2 = new TH2F("h2","Systematics, Cd, <Gamma0>",100,4.710,11.94,100,60,410);

    c1.cd(2);
    //c1_2.SetLogy();
	c1_2.SetLeftMargin(0.14);
	c1_2.SetBottomMargin(0.16);
	c1_2.SetTopMargin(0.01);
	c1_2.SetRightMargin(0.05);

	h2->GetXaxis()->SetTitle(" Neutron separation energy S_{n} (MeV)");
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetLabelFont(42);
	h2->GetXaxis()->SetTitleOffset(1.2);
    h2->GetXaxis()->CenterTitle();
	h2->GetXaxis()->SetLabelSize(0.055);
	h2->GetXaxis()->SetTitleSize(0.055);
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetYaxis()->SetLabelFont(42);
	h2->GetYaxis()->SetTitle(" #LT#Gamma_{#gamma0}#GT (meV)");
	h2->GetYaxis()->CenterTitle();
	h2->GetYaxis()->SetLabelSize(0.055);
	h2->GetYaxis()->SetTitleOffset(1.05);
	h2->GetYaxis()->SetTitleSize(0.06);
    //	c1->SetGridx();
    //	c1->SetGridy();
	h2->Draw();
    
 	TGraphErrors *gammavalues = new TGraphErrors(5,Bnexp,gammaexp,Bn_err,gammaerrexp);
  	TGraph *gammaupper = new TGraph(5,Bnexp,gammaexpup);
   	TGraph *gammalower = new TGraph(5,Bnexp,gammaexplow);
 
	gammavalues->SetMarkerStyle(21);
    gammavalues->Draw("P");
    
    gammaupper->SetFillStyle(3002);
    gammaupper->SetLineWidth(-201);
    //gammaupper->Draw("L");
    
    gammalower->SetFillStyle(3002);
    //gammalower->SetFillColor(kBlue);
    //gammalower->SetLineColor(kBlue);
    gammalower->SetLineWidth(201);
    //gammalower->Draw("L");
    
    TF1 *fit_gamma = new TF1("fit_gamma","4.64320e+01 + x*1.55772e+01 ",6.,10.3);
    //gammavalues->Fit("fit_gamma","R");
    fit_gamma->SetLineWidth(2);
    fit_gamma->SetLineColor(kBlue);
    fit_gamma->Draw("L same");
    
    TF1 *fit_gamma2 = new TF1("fit_gamma2","1.60449e+02",6.,10.3);
    fit_gamma2->SetLineWidth(1);
    fit_gamma2->SetLineStyle(2);
    //fit_gamma2->Draw("L same");
    //gammavalues->Fit("fit_gamma2","R");

    //TF1 *fit_gammalow = new TF1("fit_gammalow","[0] + x*[1] + x*x*[2]",5,11);
    //gammaupper->Fit("fit_gammalow","R");
    //TF1 *fit_gammalow = new TF1("fit_gammalow","-6.85623e+02 + x*2.1e+02  -1.1e+01*x*x",5,11);
    //fit_gammalow->SetLineWidth(1);
    //fit_gammalow->SetLineColor(kMagenta+1);
    //fit_gammalow->Draw("same L ");
    

    double gamma76Ge[1] = {193.294}; // linear fit
    double gamma76Geerr[1] = {0.24*gamma76Ge[0]};

    double gamma76Ge_talys[1] = {295}; // talys spline fit
    double gamma76Geerr_talys[1] = {gamma76Ge_talys[0] - gamma76Ge[0]};
    

    cout << " 76Ge, estimated <Gg> = " << (4.64320e+01 + Bn_76Ge[0]*1.55772e+01) << " + " << (gamma76Ge_talys[0]-(4.64320e+01 + Bn_76Ge[0]*1.55772e+01)) << " - " << gamma76Geerr[0] << " meV. " << endl;

    TGraphErrors *gamma76Ge_graph = new TGraphErrors(1,Bn_76Ge,gamma76Ge,Bn_err,gamma76Geerr);
    TGraph *gamma76Ge_graph_talys = new TGraph(1,Bn_76Ge,gamma76Ge_talys);

    TGraphAsymmErrors *gamma76Ge_asym_graph = new TGraphAsymmErrors(1,Bn_76Ge,gamma76Ge,Bn_err,Bn_err,gamma76Geerr,gamma76Geerr_talys);

    gamma76Ge_graph->SetMarkerStyle(33);
    gamma76Ge_graph->SetMarkerColor(kBlue);
    gamma76Ge_graph->SetLineColor(kBlue);
    gamma76Ge_graph->SetMarkerSize(1.6);
    //gamma76Ge_graph->Draw("p");
    
    gamma76Ge_graph_talys->SetMarkerStyle(26);
    gamma76Ge_graph_talys->SetMarkerColor(kBlue);
    gamma76Ge_graph_talys->SetLineColor(kBlue);
    gamma76Ge_graph_talys->SetMarkerSize(1.2);
    gamma76Ge_graph_talys->Draw("p");

    gamma76Ge_asym_graph->SetMarkerStyle(33);
    gamma76Ge_asym_graph->SetMarkerColor(kBlue);
    gamma76Ge_asym_graph->SetLineColor(kBlue);
    gamma76Ge_asym_graph->SetMarkerSize(1.6);
    gamma76Ge_asym_graph->Draw("p");

	t.SetTextColor(kBlack);
    t.SetTextFont(42);
    t.SetTextSize(0.057);
    t.DrawLatex(5.33,387,"(b)");
    
    t.SetTextSize(0.052);
    t.DrawLatex(5.45,150,"^{77}Ge");
    t.DrawLatex(5.85,240,"^{75}Ge");
    t.DrawLatex(6.55,240,"^{73}Ge");
    t.DrawLatex(7.15,250,"^{71}Ge");
    t.SetTextColor(kBlue);
    t.DrawLatex(8.75,245,"^{76}Ge");
    t.SetTextColor(kBlack);
    t.DrawLatex(9.8,250,"^{74}Ge");
   //double massexp[5]   = {77.00, 75.00, 73.00, 71.00, 74.00};

    TLegend *leg2 = new TLegend(0.18,0.7,0.6,0.915);
    //leg2.SetBorderSize(0);
    leg2.SetFillColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.046);
    leg2->AddEntry(gammavalues," From n.res. data ","P");
    leg2->AddEntry(gamma76Ge_asym_graph," ^{76}Ge, estimated ","P");	
    leg2->AddEntry(gamma76Ge_graph_talys," ^{76}Ge, TALYS tab.","P");	
    leg2->AddEntry(fit_gamma," linear fit of n.res.data","L");	
    leg2->Draw();
    


	c1.Update();
	c1.Print("syst_nld_Gg_Ge.pdf");

} // END of script		

