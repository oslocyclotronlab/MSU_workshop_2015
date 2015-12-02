// Root script for plotting cross sections  and reaction rates of 75Ge(n,g)76Ge
// ACL
{
	// Cleaning up (necessary when not running script in root directory)
	gROOT->Reset();
	gROOT->SetStyle("Plain");	// no fill in pads
	gStyle->SetOptTitle(0);		// no histogram title
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);		// no statistics
	
	// Deleting histograms from last run to avoid memory leak	
	m = (TH1F*)gROOT->FindObject("h2"); 
	if (m){
        m->Delete();
    }
	m = (TH1F*)gROOT->FindObject("h3"); 
	if (m){
        m->Delete();
    }
    
    
	// Opening files, cross sections
    ifstream talys_low_xfile("rp032076.tot");

    // Opening reaction rate files
    ifstream talys_low_ratefile("astrorate.g");


	// create vectors and variables    
    // cross sections
    double e_cross_talyslow[100], cross_talyslow[100];
    
    // reaction rates
    double e_rate_talyslow[100], rate_talyslow[100];

	int i,j;
	double x,y,z;
    double d1, d2, d3, d4, d5, d6, d7, d8;
    
    string line;
    line.resize(200);
    
    // cross sections
    // talys lower
    getline(talys_low_xfile,line);
    getline(talys_low_xfile,line);
    getline(talys_low_xfile,line);
    getline(talys_low_xfile,line);
    getline(talys_low_xfile,line);
	i=0;
	x=y=z=0.;
	while(talys_low_xfile){
		talys_low_xfile >> d1 >> d2;
		e_cross_talyslow[i]   = d1;
		cross_talyslow[i] = d2;
		i++;
	}
    talys_low_xfile.close();

    // reaction rates
    // talys lower
    getline(talys_low_ratefile,line);
    getline(talys_low_ratefile,line);
	i=0;
	x=y=z=0.;
	while(talys_low_ratefile){
		talys_low_ratefile >> d1 >> d2 >> d3;
		e_rate_talyslow[i]   = d1;
		rate_talyslow[i] = d2;
		i++;
	}
    talys_low_ratefile.close();

    // Making a canvas
	TCanvas *c = new TCanvas("c","Cross sections and rates, (n,g)76Ge",900,450);
    c->Divide(2,1,0,0);

	// making graphs
    TGraph *xs_talyslow_graph = new TGraphErrors(59,e_cross_talyslow,cross_talyslow);
    TGraph *rate_talyslow_graph = new TGraphErrors(30,e_rate_talyslow,rate_talyslow);
    
    // Histograms, just for making nice x&y axis
    TH2F *h2 = new TH2F("h2"," ",100,5.50E-03,1.300E+00,100,1.25832E+00,3.31595E+03); // cross section
    TH2F *h3 = new TH2F("h3"," ",100,0.08,9.9,100,2.25E+04,2.31E+09); // reaction rate
	
	c->cd(1);
	c_1->SetLeftMargin(0.17);
	c_1->SetRightMargin(0.08);
	c_1->SetBottomMargin(0.15);
	c_1->SetTopMargin(0.05);
    c_1->SetLogx();
    c_1->SetLogy();
    
    h2->GetXaxis()->SetLabelSize(0.053);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetYaxis()->SetLabelSize(0.053);
    h2->GetYaxis()->SetLabelFont(42);
    //h2->GetYaxis()->SetLabelOffset(1.0);
    h2->GetYaxis()->SetTitleSize(0.057);
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetTitleSize(0.057);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->CenterTitle();
    h2->GetYaxis()->CenterTitle();
    h2->GetYaxis()->SetTitleOffset(1.1);
    h2->GetXaxis()->SetTitle("E_{n} (MeV)");
    h2->GetXaxis()->SetLabelOffset(0.0015);
    h2->GetXaxis()->SetTitleOffset(1.15);
    h2->GetYaxis()->SetTitle("#sigma(E_{n}) (mb)");
    h2->Draw();

	
    // TALYS PREVIOUS UPPER/LOWER
    //xs_talyslow_graph->SetFillStyle(3003);
    //xs_talyslow_graph->SetLineWidth(1002);
    xs_talyslow_graph->SetLineWidth(2);
    xs_talyslow_graph->Draw("C");
    
    TLegend *leg2 = new TLegend(0.44,0.72,0.9,0.93);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);
    leg2->AddEntry(xs_talyslow_graph," EB2009 ","L");
    leg2->SetTextSize(0.042);
    leg2->SetTextFont(42);
    //leg2->Draw();
    

    TLatex t;
	t->SetTextSize(0.055);
    t->SetTextFont(42);
    t->DrawLatex(    9.00E-03,1.51595E+03," ^{75}Ge(n,#gamma)^{76}Ge");
    
    
    
    c->cd(2);
	c_2->SetLeftMargin(0.14);
	c_2->SetRightMargin(0.03);
	c_2->SetBottomMargin(0.15);
    c_2->SetTopMargin(0.05);
	c_2->SetLogx();
	c_2->SetLogy();
        
    h3->GetXaxis()->SetLabelSize(0.053);
    h3->GetYaxis()->SetLabelSize(0.053);
    h3->GetXaxis()->SetTitleSize(0.057);
    h3->GetYaxis()->SetTitleSize(0.057);
    h3->GetXaxis()->SetLabelFont(42);
    h3->GetYaxis()->SetLabelFont(42);
    h3->GetXaxis()->SetTitleFont(42);
    h3->GetYaxis()->SetTitleFont(42);
    h3->GetXaxis()->CenterTitle(42);
    h3->GetYaxis()->CenterTitle(42);
    h3->GetXaxis()->SetTitleOffset(1.15);
    h3->GetYaxis()->SetTitleOffset(1.2);
    h3->GetXaxis()->SetLabelOffset(0.0001);

    h3->GetXaxis()->SetTitle("T (10^{9} K)");
    h3->GetYaxis()->SetTitle("N_{A}#LT#sigmav#GT (cm^{3} s^{-1} mol ^{-1} )");
    //h3->GetYaxis()->SetLabelOffset(1.07);
    h3->Draw();
    
    // PREVIOUS TALYS LIMITS
    rate_talyslow_graph->SetLineWidth(2);
    rate_talyslow_graph->Draw("C");

    TLegend *leg3 = new TLegend(0.44,0.72,0.9,0.93);
	leg3->SetBorderSize(0);
	leg3->SetFillColor(0);
	leg3->AddEntry(rate_talyslow_graph," TALYS low  ","L");
	leg3->SetTextSize(0.042);
	leg3->SetTextFont(42);
	//leg3->Draw();

    t->DrawLatex(    0.12,0.8E+09," ^{75}Ge(n,#gamma)^{76}Ge rate");
    
    c->Print("cross_section_rate_76Ge.pdf");
    c->Print("cross_section_rate_76Ge.eps");
    
    
	
} // END of script
