void plotLDF(int antenna_height, int theta_d){
	std::cout << "starting" << std::endl;
	string prefix = string("/Users/wissels/Dropbox/MountainTop/harms_sims/%dkm/paper/DeltaPsi0.04/DecayID10/DecayHeight0/AntennaHeight%d/Zenith%d/Energy1e17/");
	std::cout << prefix.c_str() << std::endl;
	//int zens[7] = {60,70,75,80,85,87,89};

	TCanvas *can = new TCanvas("can", "", 800, 1200);
	can->Divide(1, 2);
	can->cd(1);
	TGraph *gMyLDF = new TGraph(80);
	//for( int theta_d_i = 0; theta_d_i < 7; theta_d_i++){
	//	int theta_d = zens[theta_d_i];
	TFile f(Form(Form("%s/out.root", prefix.c_str()), antenna_height,antenna_height*1000, theta_d));
	std::cout << f.GetName() << std::endl;
	for(int antenna=1; antenna<81; antenna++){
		std::cout << theta_d << " " << antenna << std::endl;
		TGraph *g1 = (TGraph*)f.Get(Form("Antenna%d_time_x", antenna));
		TGraph *g2 = (TGraph*)f.Get(Form("Antenna%d_time_y", antenna));
		TGraph *g3 = (TGraph*)f.Get(Form("Antenna%d_time_z", antenna));
		TGraph *g4 = (TGraph*)f.Get(Form("Antenna%d_freq_x", antenna));
		TGraph *g5 = (TGraph*)f.Get(Form("Antenna%d_freq_y", antenna));
		TGraph *g6 = (TGraph*)f.Get(Form("Antenna%d_freq_z", antenna));

		g2->Draw(Form("L%s", antenna==1 ? "A":"SAME" ));
		double max = -9999;
		for( int iPoint=0; iPoint < g2->GetN(); iPoint++){
			double x, y;
			g2->GetPoint(iPoint, x, y);
			if( antenna == 2)
				std::cout << iPoint << " " << x << " " << y << " " << max<< std::endl;
			if( TMath::Abs(y) > max){
				max = TMath::Abs(y);
			}
		}
	 	std::cout << antenna-1 << " "<< ( (float) antenna)*0.04 << " " << max << " " << g2->GetN() << std::endl;	
		gMyLDF->SetPoint(antenna-1, ( (float) antenna-1)*0.04, max);
	}
	can->cd(2);
	gMyLDF->SetLineColor(kRed);
	gMyLDF->Draw("AL");
	TGraph *gHarmsLDF = (TGraph*)f.Get("grLDF");
	gHarmsLDF->Draw("PSAME");
}

void compareLDFs(int antenna_height){
	string prefix = string("/Users/wissels/Dropbox/MountainTop/harms_sims/%dkm/paper/DeltaPsi0.04/DecayID10/DecayHeight0/AntennaHeight%d/Zenith%d/Energy1e17/");
	
	int theta_d[9] = {55,60,65,70,75,80,85,87,89};
	for( int it=0; it<9; it++){
		std::cout << theta_d[it] << std::endl;
		TFile f(Form(Form("%s/out.root", prefix.c_str()), antenna_height,antenna_height*1000, theta_d[it]));
		TGraph *gHarmLDF = (TGraph*)f.Get("grLDF");
		gHarmLDF->SetLineColor(it+1);
		gHarmLDF->SetFillColor(kWhite);
		gHarmLDF->SetTitle(Form("%d^{#circ}", theta_d[it]));
		gHarmLDF->Draw(Form("L%s", it==0 ? "A" : "SAME"));
	}
}

void writeTGraph(TGraph *g, char *finame){
	//std::ofstream file(finame, ios::out);
    FILE* pFile;
    pFile = fopen(finame, "w");
    std::cout << "writing graph " << finame << std::endl;
	for(int ip=0; ip<g->GetN(); ip++){
		double x, y;
		g->GetPoint(ip, x, y);
        	std::fprintf(pFile, "%1.10e,%1.10e\n", x,y);
		//file << x << ", "<< y << std::endl;
	}
	fclose(pFile);
}


int dumpData(int antenna_height){
	std::cout << "starting" << std::endl;
	string prefix = string("/Users/wissels/Dropbox/MountainTop/harms_sims/%dkm/paper/DeltaPsi0.04/DecayID10/DecayHeight0/AntennaHeight%d/Zenith%d/Energy1e17/");
	std::cout << prefix.c_str() << std::endl;
	int zens[9] = {55,60,65,70,75,80,85,87,89};
	for( int theta_d_i = 0; theta_d_i < 9; theta_d_i++){
		int theta_d = zens[theta_d_i];
		TFile f(Form(Form("%s/out.root", prefix.c_str()), antenna_height,antenna_height*1000, theta_d));
		std::cout << f.GetName() << std::endl;
		for(int antenna=1; antenna<81; antenna++){
			std::cout << theta_d << " " << antenna << std::endl;
			TGraph *g1 = (TGraph*)f.Get(Form("Antenna%d_time_x", antenna));
			TGraph *g2 = (TGraph*)f.Get(Form("Antenna%d_time_y", antenna));
			TGraph *g3 = (TGraph*)f.Get(Form("Antenna%d_time_z", antenna));
			TGraph *g4 = (TGraph*)f.Get(Form("Antenna%d_freq_x", antenna));
			TGraph *g5 = (TGraph*)f.Get(Form("Antenna%d_freq_y", antenna));
			TGraph *g6 = (TGraph*)f.Get(Form("Antenna%d_freq_z", antenna));
			std::cout << 1 << std::endl;	
			writeTGraph(g1, Form(Form("%s/csv/antenna_%d_time_x.csv", prefix.c_str(), antenna), antenna_height,antenna_height*1000, theta_d));
			writeTGraph(g2, Form(Form("%s/csv/antenna_%d_time_y.csv", prefix.c_str(), antenna), antenna_height,antenna_height*1000, theta_d));
 			writeTGraph(g3, Form(Form("%s/csv/antenna_%d_time_z.csv", prefix.c_str(), antenna), antenna_height,antenna_height*1000, theta_d));
			writeTGraph(g4, Form(Form("%s/csv/antenna_%d_freq_x.csv", prefix.c_str(), antenna), antenna_height,antenna_height*1000, theta_d));
			writeTGraph(g5, Form(Form("%s/csv/antenna_%d_freq_y.csv", prefix.c_str(), antenna), antenna_height,antenna_height*1000, theta_d));
 			writeTGraph(g6, Form(Form("%s/csv/antenna_%d_freq_z.csv", prefix.c_str(), antenna), antenna_height,antenna_height*1000, theta_d));
			std::cout << 1 << std::endl;
		}
		
	}
   	return 0;
}
