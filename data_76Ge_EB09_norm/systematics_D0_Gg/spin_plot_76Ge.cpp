// Spin distribution 76Ge, EB09, and fraction propulated in beta decay
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
	

    // Read robin spin distribution (BSFG, EB2009)
    ifstream spinfile("spindis_76Ge.rbn");
    
	// Vectors
    double ex[50], spin[50];
    int i,j;
    double d1,d2;
    double sum_intensity;
    double sum_1_3;
    double sum_0_4;
    
    i=0;
    while(spinfile){
        spinfile >> d1 >> d2;
        ex[i]   = d1;
        spin[i] = d2;
        sum_intensity += d2;
        if(ex[i]>=1 && ex[i]<=3){
            sum_1_3 += d2;
        }
        if(ex[i]>=0 && ex[i]<=4){
            sum_0_4 += d2;
        }
        i++;
    }
    
    cout << " Sum of all intensities = " << sum_intensity << endl;
    cout << " Sum of J=1-3 intensities = " << sum_1_3 << endl;
    cout << " Ratio, J=1-3/all = " << sum_1_3/sum_intensity << endl;
    cout << " Sum of J=0-4 intensities = " << sum_0_4 << endl;
    cout << " Ratio, J=0-4/all = " << sum_0_4/sum_intensity << endl;
    
    //Assumed experimental distribution
    double spin_ex[3] = {1.,2,3};
    double spin_I[3]  = {0.333,0.333,0.333};
    
 	// Making graphs			
	TGraph *spindistr = new TGraph(41,ex,spin);
    TGraph *exp_spindistr = new TGraph(3,spin_ex,spin_I);

	// Making a canvas
	TCanvas *c1 = new TCanvas("c1","Spin distribution, EB09, 76Ge",600,500);
    	
    c1.cd();
	//c1->SetLogy();
    c1.SetLeftMargin(0.14);
    c1.SetBottomMargin(0.12);
	spindistr->GetXaxis()->SetTitle(" spin J");
    spindistr->GetXaxis()->SetTitleFont(42);
	spindistr->GetXaxis()->CenterTitle();
	spindistr->GetXaxis()->SetLabelSize(0.04);
	spindistr->GetXaxis()->SetTitleOffset(1.2);
    spindistr->GetXaxis()->SetLabelFont(42);
	spindistr->GetYaxis()->SetTitle(" relative intensity per spin");
    spindistr->GetYaxis()->SetTitleFont(42);
	spindistr->GetYaxis()->CenterTitle();
	spindistr->GetYaxis()->SetLabelSize(0.04);
	spindistr->GetYaxis()->SetLabelSize(0.04);
    spindistr->GetYaxis()->SetLabelFont(42);
	spindistr->GetYaxis()->SetTitleOffset(1.5);
    spindistr->GetXaxis()->SetRangeUser(0.0,17.);
    spindistr->GetYaxis()->SetRangeUser(0.0,0.25);
    spindistr->SetMarkerStyle(21);
	spindistr->Draw("ACP");
        
    //exp_spindistr->Draw("L same");

 
	TLatex t;
	t.SetTextSize(0.04);
    t.SetTextFont(42);
	t.DrawLatex(6.17,0.23,"^{76}Ge spin distr., EB09");


    TLegend *leg = new TLegend(0.16,0.69,0.57,0.88);
    //leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.038);
    //leg->AddEntry(d2rho," From neutron res. data ","P");
    //leg->Draw();

	c1.Update();
	c1.Print("spindis_76Ge.pdf");

} // END of script

