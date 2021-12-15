// A C++ program for interfacing with MicroOmegas.
// Returns phenomenologically interesting results of the spectrum, including:
// B_s -> mu+ mu-
// Omega h^2
// m_h
// g-2
// B_mu -> tau nu
// SUSY Scale masses
// SUSY Scale couplings
// CDM-Nucleon scattering rates
// Rl23
// LSP Mass

// We write a class that contains a function to calculate these values for an input point in the (p)MSSM.

#include "dm_compute.h" 

using namespace std;

DM::DM(string micro_omega_dir, string softSUSY_dir, string temporary_dir)
{
    DM::ss_dir = softSUSY_dir;
    DM::mo_dir = micro_omega_dir;
    DM::temp_dir = temporary_dir;
};

string DM::getError(string tempDir, string threadID)
{
    ifstream slhaInput3;
    string line;
    slhaInput3.open(tempDir + "/slha_output" + threadID+".txt");
    vector<string> split_line;
    while(getline(slhaInput3, line))
        if(line.find("SPINFO") != string::npos)
        {
            while(getline(slhaInput3, line))
            {
                if(line.find("Block") != string::npos)
                    break;
                boost::trim(line);
                if(line[0] == '4')
                {
                    vector<string> split1, split2;
                    boost::split(split1, line, boost::is_any_of("["));
                    boost::split(split2, split1[1], boost::is_any_of("]"));
                    return split2[0];		
                }				
                if(line[0] == '3')
                {
                    vector<string> split1;
                    boost::split(split1 , line, boost::is_any_of("#"));
                    return split1[0].substr(2);
                }
            }
            break;
	}
    return "Could not find error";
};

string DM::format_string(string x)
{
	boost::erase_all(x, "MSSM");
	boost::erase_all(x, "(Q)");
	boost::erase_all(x, "^");
	boost::erase_all(x, "_");
	boost::replace_all(x, "'", "prime");
	return x;
};

// What a mess
// This parses the slha formatted output file to grab all parameters at all scales
//each element in returned vector is a different scale
// each element in the map corresponds to a parameter and its value (see phypar in .h for keys)
vector<map<string, double>> DM::get_running_params(string x)
{
	bool entered = false;
	vector<map<string,double>> all_scales;
	vector<string> output;
	boost::split(output, x, boost::is_any_of("\n"));
	for(int i=0; i<output.size(); i++)
		if(output[i].find("Block gauge") != string::npos)
		{
			map<string,double> running_params;
			vector<string> split_line;
			boost::split(split_line, output[i], boost::is_any_of(" "));
			int q_index;
			for(int j=0; j<split_line.size(); j++)
				if(split_line[j] == "Q=")
					q_index = j;
			double q_value = atof(split_line[q_index + 1].c_str());
			while(i != output.size()-2 && output[i+1].find("Block gauge") == string::npos)
			{
				string line = format_string(output[i+1]);
                                vector<string> new_split_line;
				boost::split(new_split_line, line, boost::is_any_of(" "));
				for(int j=0; j<phypar.size(); j++)
				{
					int hashtag_index = -1;
					bool param_bool = false;
                                        int split_line_length = new_split_line.size();
					for(int k=0; k<split_line_length;k++)
					{
                                            string element = new_split_line[k];
                                            boost::trim(element);
					    if(element == "#") hashtag_index = k;
                                            // Erase empty elements from the vector
                                            if(element == "")
                                            {
                                                new_split_line.erase(new_split_line.begin() + k);
                                                k -= 1; split_line_length -=1;
                                            }
					    if(element == phypar[j]) param_bool = true;
					    if(param_bool && hashtag_index != -1)
                                            {
                                                /*if(phypar[j]=="Ac" && line.find("Ac") != string::npos)
                                                { 
                                                    for(int n=0; n <new_split_line.size();n++) cout << new_split_line[n] << ",";
                                                    cout << endl;
                                                    cout << hashtag_index << " " << new_split_line[hashtag_index-1] << endl;
                                                }*/                                               
                                                element = new_split_line[hashtag_index-1];
                                                boost::trim(element);
						running_params.insert(pair<string, double>(phypar[j], atof(element.c_str())));
                                            }
					}    
				}
				running_params.insert(pair<string,double>("Q", q_value));
				entered = true;
				i++;		
			}
			all_scales.push_back(running_params);
		}		
	return all_scales;
};

map<string, double> DM::get_weak_params(string x)
{
    map<string, double> weak_params;
    vector<string> output;
    boost::split(output, x, boost::is_any_of("\n"));
    for(int i=0; i<output.size();i++)
    {
        if(output[i].find("Block MASS") != string::npos)
            while(i != output.size()-2 && output[i+1].find("Block alpha") == string::npos)
            {
                string line = format_string(output[i+1]);
                vector<string> split_line;
                boost::split(split_line, line, boost::is_any_of(" "));
                //for(int k=0; k < split_line.size();k++) cout << split_line[k] << endl;
                for(int j=0; j<weakpar.size(); j++)
                {
                    int hashtag_index = -1;
                    bool param_bool = false;
                    for(int k=0; k < split_line.size();k++)
                    {
                        string element = split_line[k];
                        boost::trim(element);
                        if(element=="#") hashtag_index = k;
                        //if(element=="") split_line.erase(split_line.begin()+k);
                        if(element==weakpar[j]) param_bool = true;
                        //cout << param_bool << " " << hashtag_index << endl;
                        if(param_bool && hashtag_index != -1)
                        {
                            element = split_line[hashtag_index-3];
                            boost::trim(element);
                            weak_params.insert(pair<string, double>(weakpar[j], atof(element.c_str())));
                        }
                    }
                }
                i++;
            }
        if(output[i].find("Block nmix") != string::npos)
        {
            // Add neutralino mixing matrix scale
            vector<string> split_header;
            boost::split(split_header, output[i], boost::is_any_of(" "));
            int scale_index;
            for(int j=0; j<split_header.size(); j++)
                if(split_header[j] == "matrix")
                    scale_index = j+1;
            vector<string> split_scale;
            boost::split(split_scale, split_header[scale_index], boost::is_any_of("="));
            double Qnmix = atof(split_scale[1].c_str());
            weak_params.insert(pair<string,double>("Qnmix", Qnmix));
            // Add neutralino mixing matrix elements for LSP components
            while(i != output.size()-2 && output[i+1].find("Block Umix") == string::npos)
            {
                string line = format_string(output[i+1]);
                vector<string> split_line;
                boost::split(split_line, line, boost::is_any_of(" "));
                for(int j=0; j<weakpar.size(); j++)
                {
                    int hashtag_index = -1;
                    bool param_bool = false;
                    for(int k=0; k < split_line.size();k++)
                    {
                        string element = split_line[k];
                        boost::trim(element);
                        if(element=="#") hashtag_index = k;
                        if(element=="") split_line.erase(split_line.begin()+k);
                        if(element==weakpar[j]) param_bool = true;
                        if(param_bool && hashtag_index != -1)
                        {
                            element = split_line[hashtag_index-1];
                            boost::trim(element);
                            weak_params.insert(pair<string, double>(weakpar[j], atof(element.c_str())));
                        }
                    }
                }
                i++;
            }
        }
    }
    return weak_params;
}

vector<double> DM::get_gut_scale_gauge_couplings(string x, vector<map<string,double>> running_params)
{
    double mx_scale;
    vector<string> output;
    vector<double> gauge_couplings = {-1,-1,-1,-1};
    boost::split(output, x, boost::is_any_of("\n"));
    for(int i=0; i<output.size(); i++)
    {
        if(output[i].find("Block EXTPAR") != string::npos)
        {
            vector<string> split_line;
            // Line after Block EXTPAR is always MX scale
            boost::split(split_line, output[i+1], boost::is_any_of(" "));
            for(int k=0; k<split_line.size();k++)
                if(split_line[k] == "#")
                    mx_scale = atof(split_line[k-2].c_str());
            break;
        }
    }
    // Loop through running params to get gauge couplings at scales
    // just below and just above mx scale, then linearly interpolate.
    for(int i=0; i<running_params.size(); i++)
    {
        if(running_params[i]["Q"] < mx_scale && running_params[i+1]["Q"] > mx_scale)
        {
            double lower_Q, upper_Q, lower_g, lower_gprime, lower_g3, upper_g, upper_gprime, upper_g3, epsilon3;
            lower_Q = running_params[i]["Q"];
            upper_Q = running_params[i+1]["Q"];
            lower_g = running_params[i]["g"];
            lower_gprime = running_params[i]["gprime"];
            lower_g3 = running_params[i]["g3"];
            upper_g = running_params[i+1]["g"];
            upper_gprime = running_params[i+1]["gprime"];
            upper_g3 = running_params[i+1]["g3"];
            // linear interpolation
            double g_mx = (upper_g - lower_g) / log10(upper_Q/lower_Q) * log10(mx_scale / lower_Q) + lower_g;
            double gprime_mx = (upper_gprime - lower_gprime) / log10(upper_Q/lower_Q) * log10(mx_scale / lower_Q) + lower_gprime;
            double g3_mx = (upper_g3 - lower_g3) / log10(upper_Q/lower_Q) * log10(mx_scale / lower_Q) + lower_g3;
            // Grand unified normalization, constant g2 = sqrt(5/3) g' 
            gauge_couplings[0] = mx_scale;
            gauge_couplings[1] =  sqrt(5./3.) * gprime_mx;
            gauge_couplings[2] = g_mx;
            gauge_couplings[3] = g3_mx;
        }
    }
    return gauge_couplings;
}

string DM::run_calc(vector<double> input_point, string thread_id)
{
	double m0, m12, a0, tanb, inm1, inm2, inm3, ina0, inat, inab, inata, inmel, inmmul, inmtal, inmer, inmmur, inmtar, inmq1l, inmq2l, inmq3l, inmur, inmcr, inmtr, inmdr, inmsr, inmbr, inmmu, inma, inmhu2, inmhd2;
	int sgnmu;
	int num_inputs = input_point.size();
	ostringstream content;
	// 2 Dimensional theory with X_t = 0, can compare to published results.
	//INPUTS: {m0, m12}
	if(num_inputs == 3)
	{
		m0 = input_point[0];
		m12 = input_point[1];
        	sgnmu = round(input_point[2]);
		tanb = 9.0;
		// Necessary for X_t = 0
		a0 = 0.0;
		inmhu2 = pow(m0,2);
		inmhd2 = pow(m0,2);
		inm1 = m12;
		inm2 = m12;
		inm3 = m12;
		ina0 = a0;
		inat = a0;
		inab = a0;
		inata = a0;
		inmel = m0;
		inmmul = m0;
		inmtal = m0;
		inmer = m0;
		inmmur = m0;
		inmtar = m0;
		inmq1l = m0;
		inmq2l = m0;
		inmq3l = m0;
		inmur = m0;
		inmcr = m0;
		inmtr = m0;
		inmdr = m0;
		inmsr = m0;
		inmbr = m0;
	}
	// 4 Dimensional, mSUGRA.
	// INPUTS: {m0, m12, a0, tanb, sign mu}
	if(num_inputs == 5)
	{
		m0 = input_point[0];
		m12 = input_point[1];
		a0 = input_point[2];
		tanb = input_point[3];
        	sgnmu = round(input_point[4]);
		// See Martin SUSY Primer EQ 8.1.8, 8.1.9
		inmhu2 = pow(m0,2);
		inmhd2 = pow(m0,2);
		inm1 = m12;
		inm2 = m12;
		inm3 = m12;
		ina0 = a0;
		inat = a0;
		inab = a0;
		inata = a0;
		inmel = m0;
		inmmul = m0;
		inmtal = m0;
		inmer = m0;
		inmmur = m0;
		inmtar = m0;
		inmq1l = m0;
		inmq2l = m0;
		inmq3l = m0;
		inmur = m0;
		inmcr = m0;
		inmtr = m0;
		inmdr = m0;
		inmsr = m0;
		inmbr = m0;
	}
	// INPUTS: {m1, m2, m3, mmu, mA, at, ab, ata, mel, mtal, mer, mtar, mq1l, mq3l, mur, mtr, mdr, mbr, tanb}
	if(num_inputs == 19)
	{
		ina0 = 0;
		tanb = input_point[18];
		inm1 = input_point[0];
		inm2 = input_point[1];
		inm3 = input_point[2];
		inmmu = input_point[3];
		inma = input_point[4];
		inat = input_point[5];
		inab = input_point[6];
		inata = input_point[7];
		inmel = input_point[8];
		inmmul = inmel;
		inmtal = input_point[9];
		inmer = input_point[10];
		inmmur = inmer;
		inmtar = input_point[11];
		inmq1l = input_point[12];
		inmq2l = inmq1l;
		inmq3l = input_point[13];
		inmur = input_point[14];
		inmcr = inmur;
		inmtr = input_point[15];
		inmdr = input_point[16];
		inmsr = inmdr;
		inmbr = input_point[17];
	}
	ostringstream ss;
	if(num_inputs ==3 || num_inputs==5)
		ss << "BLOCK MODSEL      # Select Model\n"
		   << " 1 1              # mSUGRA\n"
		   << " 11 32            # number of points to evaluate in the spectrum\n"
		   << " 12 10.000e+17    # Max scale to evaluate couplings\n" 
		   << "BLOCK SMINPUTS    # Standard Model Inputs\n" 
		   << " 1 1.27934e+02    # alpha^(-1) SM MSbar (MZ)\n" 
		   << " 2 1.16637e-05    # G_Fermi\n" 
		   << " 3 1.17200e-01    # alpha_s (MZ) SM MSbar\n"
		   << " 4 9.11876e+01    # MZ(Pole)\n" 
		   << " 5 4.25000e+00    # mb(mb) SM MSbar\n"
		   << " 6 1.74200e+02    # mtop(pole)\n" 
		   << " 7 1.77700e+00    # mtau (pole)\n" 
		   << "BLOCK MINPAR      # Input Parameters\n" 
		   << " 3 " << scientific << tanb << "\n" 
		   << " 4 " << sgnmu << "\n" 
		   << " 5 " << ina0 << "\n" 
		   << "BLOCK EXTPAR      # Optional Input Parameters\n" 
		   << " 1 " << inm1 << "\n" 
		   << " 2 " << inm2 << "\n" 
		   << " 3 " << inm3 << "\n" 
		   << " 11 " << inat << "\n" 
		   << " 12 " << inab << "\n" 
		   << " 13 " << inata << "\n" 
		   << " 21 " << inmhu2 << "\n" 
		   << " 22 " << inmhd2 << "\n" 
		   << " 31 " << inmel << "\n" 
		   << " 32 " << inmmul << "\n" 
		   << " 33 " << inmtal << "\n" 
		   << " 34 " << inmer << "\n" 
		   << " 35 " << inmmur << "\n" 
		   << " 36 " << inmtar << "\n" 
		   << " 41 " << inmq1l << "\n" 
		   << " 42 " << inmq2l << "\n" 
		   << " 43 " << inmq3l << "\n" 
		   << " 44 " << inmur << "\n" 
		   << " 45 " << inmcr << "\n" 
		   << " 46 " << inmtr << "\n" 
		   << " 47 " << inmdr << "\n" 
		   << " 48 " << inmsr << "\n" 
		   << " 49 " << inmbr << "\n";
	else
		ss << "BLOCK MODSEL      # Select Model\n"
		   << " 1 0              # Universal\n"
		   << " 11 32            # number of points to evaluate in the spectrum\n" 
		   << " 12 1.000e+17    # Max scale to evaluate couplings\n"
		   << "BLOCK SMINPUTS    # Standard Model Inputs\n" 
		   << " 1 1.27934e+02    # alpha^(-1) SM MSbar (MZ)\n" 
		   << " 2 1.16637e-05    # G_Fermi\n" 
		   << " 3 1.17200e-01    # alpha_s (MZ) SM MSbar\n"
		   << " 4 9.11876e+01    # MZ(Pole)\n" 
		   << " 5 4.25000e+00    # mb(mb) SM MSbar\n"
		   << " 6 1.74200e+02    # mtop(pole)\n" 
		   << " 7 1.77700e+00    # mtau (pole)\n" 
                   << "BLOCK MINPAR      # Input Parameters\n" 
		   << " 3 " << scientific << tanb << "\n" 
		   << "BLOCK EXTPAR      # Optional Input Parameters\n"
                   << " 0 -1 \n # SUSY scale parameterization \n" 
		   << " 1 " << inm1 << "\n"
		   << " 2 " << inm2 << "\n" 
		   << " 3 " << inm3 << "\n" 
		   << " 11 " << inat << "\n" 
		   << " 12 " << inab << "\n" 
		   << " 13 " << inata << "\n" 
		   << " 23 " << inmmu << "\n" 
		   << " 26 " << inma << "\n" 
		   << " 31 " << inmel << "\n" 
		   << " 32 " << inmmul << "\n" 
		   << " 33 " << inmtal << "\n" 
		   << " 34 " << inmer << "\n" 
		   << " 35 " << inmmur << "\n" 
		   << " 36 " << inmtar << "\n" 
		   << " 41 " << inmq1l << "\n" 
		   << " 42 " << inmq2l << "\n" 
		   << " 43 " << inmq3l << "\n" 
		   << " 44 " << inmur << "\n" 
		   << " 45 " << inmcr << "\n" 
		   << " 46 " << inmtr << "\n" 
		   << " 47 " << inmdr << "\n" 
		   << " 48 " << inmsr << "\n" 
		   << " 49 " << inmbr << "\n";
	// Write input to a file.
	ofstream writefile;
	writefile.open(temp_dir + "/slha_input" + thread_id + ".txt");
	writefile << ss.str();
	writefile.close();
	// run the softsusy command
	string cmd = ss_dir + string("/softpoint.x leshouches <") + temp_dir + string("/slha_input") + thread_id + string(".txt >") + temp_dir + string("/slha_output") + thread_id +string(".txt");
	int error = system(cmd.c_str());
	if(error != 0)
	{
		cout << getError(temp_dir, thread_id) << endl;
	        cmd = string("rm ") + temp_dir + string("/slha_input") + thread_id + string(".txt");
	        system(cmd.c_str());
	        cmd = string("rm ") + temp_dir + string("/slha_output") + thread_id + string(".txt");
	        system(cmd.c_str());
		return "Error";
	}	
	// delete the input file
	cmd = string("rm ") + temp_dir + string("/slha_input") + thread_id + string(".txt");
	//system(cmd.c_str());
	// read the output file
	ifstream ss_readfile; ostringstream ss_content; string line;
	ss_readfile.open(temp_dir + "/slha_output" +thread_id + ".txt");
	while(!ss_readfile.eof())
	{	
		getline(ss_readfile, line);
		ss_content << line << "\n";
	}
        ss_readfile.close();
	// Get the running parameters from softsusy output.
	vector<map<string,double>> running_params = get_running_params(ss_content.str());
        // Get the weak scale masses, couplings, and mixings from softsusy output
        map<string, double> weak_params = get_weak_params(ss_content.str());
	/// Interpolate running parameters to extract GUT scale gauge couplings.
        vector<double> gut_scale_gauge_couplings = get_gut_scale_gauge_couplings(ss_content.str(), running_params);
	// Parameters at a single scale for FlexibleEFTHiggs input
	map<string,double> ran_params = running_params[0];
        //for (auto const& pair: ran_params) std::cout << "{" << pair.first << ": " << pair.second << "}\n";
        //for (auto const& pair: weak_params) std::cout << "{" << pair.first << ": " << pair.second << "}\n";
        
	// run the MicrOmegas command
	cmd = mo_dir + string("/main ") + temp_dir + string("/slha_output") + thread_id + string(".txt ") + temp_dir + " " + thread_id;
	error = system(cmd.c_str());
	if(error != 0)
	{
		cout << "There was an error with micrOmegas." << endl;
		return "Error";
	}
        // Add to content (Weak scale soft masses, weak scale masses, micrOmegas result)
        for(int i=0; i <phypar.size(); i++) content << ran_params[phypar[i]] << ",";
        for(int i=0; i <weakpar.size(); i++) content << weak_params[weakpar[i]] << ",";
	// read the micrOmegas output file
	ifstream readfile; vector<string> micrOmegas_content;
	readfile.open(temp_dir + "/micrOmegas_output" +thread_id + ".txt");
	while(!readfile.eof())
	{	
		getline(readfile, line);
		content << line;
	}
	//delete the micrOmegas output file
	cmd = string("rm ") + temp_dir + string("/micrOmegas_output") + thread_id + string(".txt");
	system(cmd.c_str());
	// delete the SoftSUSY output file
	cmd = string("rm ") + temp_dir + string("/slha_output") + thread_id + string(".txt");
	system(cmd.c_str());
//  Omega, g-2, b->sgamma, b->sgamma(SM), B->tau+nu, B_s -> mu+mu-, D->taunu, D_s+->mu+nu, dRho, RL23, m_dm, SIPsig,
//SDPsig, SINsig, SDNsig, {Gluino, mA, M_H, M_h, M_H+-, eL, muL tauL, eR, muR, tauR, nu_e, nu_mu, nu_tau, uL, cL, tL,
// uR, tR, dL, sL, tL, dR, bR, tanb, neutralino 1, neutralino 2, neutralino 3, neutralino 4, chargino 1, chargino 2}(MZ)
        content << "," << gut_scale_gauge_couplings[0] << "," << gut_scale_gauge_couplings[1] << "," << gut_scale_gauge_couplings[2] << "," << gut_scale_gauge_couplings[3];

	return content.str();
};

