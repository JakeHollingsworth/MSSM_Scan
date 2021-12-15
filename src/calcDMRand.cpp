#include "dm_compute.h"
#include <chrono>
#include <thread>
#include <unistd.h>
#include <limits.h>

/*
Randomly sample mssm parameter space and calculate Higgs masses.
*/

vector<string> readConfig()
{
	vector<string> parameters = {};
	ifstream configFile;
	// This .cpp file uses relative paths!
	configFile.open("../config.txt");
	if(!configFile)
	{
		cout << "The parent directory must contain a file config.txt." << endl;
		exit(1);
	}
	string line;
	while(getline(configFile, line))
		if(line.length()!= 0 && line[0] != '#')
		{
			string indicator = line.substr(0,2);
			vector<string> split_content;
			string content;
			boost::split(split_content, line, boost::is_any_of(" "))[1];
			content = split_content[1];
			boost::trim(content);
			if(indicator != "##")
                                parameters.push_back(content);
		}
	configFile.close();
	return parameters;
}

int main(int argc, char* argv[])
{
	// Read in parameters from config file
	vector<string> parameters = readConfig();
	string micrOmegasDir = parameters[0];
        string softSUSYDir = parameters[1];
	string param_space = parameters[2];
	string finalOutputDir = parameters[3];
	string temporaryDir_root = parameters[4];
	int number_threads = stoi(parameters[5]);
	int num_to_generate = stoi(parameters[6]);
        // Get the batch ID number to ensure unique seeds
        int batchNum;
        if(argc < 1)
	{
	    cout << "You need to pass a seed with this program" << endl;
            return 0;
	}
        else
            batchNum = atoi(argv[1]);
        string temporaryDir = temporaryDir_root + "/" + "working_directory_"+to_string(batchNum);
        string create_temporaryDir = string("mkdir " + temporaryDir);
        system(create_temporaryDir.c_str());
	DM dm(micrOmegasDir, softSUSYDir, temporaryDir);
        // For pMSSM bounds, we copy those used here https://arxiv.org/pdf/1407.4130.pdf
        // For cMSSM bounds, we copy those used here: https://arxiv.org/pdf/1305.2914.pdf
	// Only used for 2D and 4D theories
	vector<double> m0_bounds = {0.0, 10000.0};
	vector<double> m12_bounds = {0.0, 10000.0};
        // sign Mu randomly chosen below.
	// Used in all theories
	vector<double> tanb_bounds_pmssm = {1., 60.0};
        vector<double> tanb_bounds_cmssm = {1.5, 50.0};
	// Only used for 19D theory (pMSSM).
        // sign determined below
	vector<double> mmu_bounds = {100.0, 4000.0};
	vector<double> ma_bounds = {100.0, 4000.0};
        // sign determined below
	vector<double> m1_bounds = {50.0, 4000.0};
        // sign determined below
	vector<double> m2_bounds = {100.0, 4000.0};
	vector<double> m3_bounds = {400.0, 4000.0};
	vector<double> at_bounds = {-4000.0, 4000.0};
	vector<double> ab_bounds = {-4000.0, 4000.0};
	vector<double> ata_bounds = {-4000.0, 4000.0};
	vector<double> mel_bounds = {100.0, 4000.0};
	vector<double> mtal_bounds = {100.0, 4000.0};
	vector<double> mer_bounds = {100.0, 4000.0};
	vector<double> mtar_bounds = {100.0, 4000.0};
	vector<double> mq1l_bounds = {400.0, 4000.0};
	vector<double> mq3l_bounds = {200.0, 4000.0};
	vector<double> mur_bounds = {400.0, 4000.0};
	vector<double> mtr_bounds = {200.0, 4000.0};
	vector<double> mdr_bounds = {400.0, 4000.0};
	vector<double> mbr_bounds = {200.0, 4000.0};
	// Randomly generated points within the bounds defined above.
	stringstream ss;
   	uniform_real_distribution<double> unif(0,1);
	default_random_engine eng;
        // This ensures a unique seed.
	eng.seed(time(0) + batchNum);
	int start_s=time(NULL);
	cout << "Beginning loop" << endl;
        stringstream output;
        string fail_key = "-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1";
	#pragma omp parallel for num_threads(number_threads)
	for(int i=0; i < num_to_generate; i++)
	{
		string thread_id = to_string(omp_get_thread_num());
		double sgnmu=1.0;
		double m0, m12, a0, tanb_pmssm, tanb_cmssm, m1, m2, m3, at, ab, ata, mel, mtal, mer, mtar, mq1l, mq3l, mur, mtr, mdr, mbr, mmu, ma;
		ofstream thread_file;	
 		thread_file.open(temporaryDir + "/outputjh" + thread_id + ".txt", ios_base::app);
		if(unif(eng) < .5)
			sgnmu=-1.0;
		m0 = unif(eng) * (m0_bounds[1] - m0_bounds[0]) + m0_bounds[0];
		m12 = unif(eng) * (m12_bounds[1] - m12_bounds[0]) + m12_bounds[0];
                vector<double> a0_bounds = {-6 * m0, 6 * m0};
		a0 = unif(eng) * (a0_bounds[1] - a0_bounds[0]) + a0_bounds[0];
		tanb_pmssm = unif(eng) * (tanb_bounds_pmssm[1] - tanb_bounds_pmssm[0]) + tanb_bounds_pmssm[0];
		tanb_cmssm = unif(eng) * (tanb_bounds_cmssm[1] - tanb_bounds_cmssm[0]) + tanb_bounds_cmssm[0];
		mmu = unif(eng) * (mmu_bounds[1] - mmu_bounds[0]) + mmu_bounds[0];
                if(unif(eng) < .5)
                    mmu = -mmu;
	        ma = unif(eng) * (ma_bounds[1] - ma_bounds[0]) + ma_bounds[0];
		m1 = unif(eng) * (m1_bounds[1] - m1_bounds[0]) + m1_bounds[0];
                if(unif(eng) < .5)
                    m1 = -m1;
		m2 = unif(eng) * (m2_bounds[1] - m2_bounds[0]) + m2_bounds[0];
                if(unif(eng) < .5)
                    m2 = -m2;
		m3 = unif(eng) * (m3_bounds[1] - m3_bounds[0]) + m3_bounds[0];
		at = unif(eng) * (at_bounds[1] - at_bounds[0]) + at_bounds[0];
		ab = unif(eng) * (ab_bounds[1] - ab_bounds[0]) + ab_bounds[0];
		ata = unif(eng) * (ata_bounds[1] - ata_bounds[0]) + ata_bounds[0];
		mel = unif(eng) * (mel_bounds[1] - mel_bounds[0]) + mel_bounds[0];
		mtal = unif(eng) * (mtal_bounds[1] - mtal_bounds[0]) + mtal_bounds[0];
		mer = unif(eng) * (mer_bounds[1] - mer_bounds[0]) + mer_bounds[0];
		mtar = unif(eng) * (mtar_bounds[1] - mtar_bounds[0]) + mtar_bounds[0];
		mq1l = unif(eng) * (mq1l_bounds[1] - mq1l_bounds[0]) + mq1l_bounds[0];
		mq3l = unif(eng) * (mq3l_bounds[1] - mq3l_bounds[0]) + mq3l_bounds[0];
		mur = unif(eng) * (mur_bounds[1] - mur_bounds[0]) + mur_bounds[0];
		mtr = unif(eng) * (mtr_bounds[1] - mtr_bounds[0]) + mtr_bounds[0];
		mdr = unif(eng) * (mdr_bounds[1] - mdr_bounds[0]) + mdr_bounds[0];
		mbr = unif(eng) * (mbr_bounds[1] - mbr_bounds[0]) + mbr_bounds[0]; 
		vector<double> mssm_point;
		string result;
		if(param_space == "2D")
		{
			mssm_point = {m0, m12, sgnmu};
			result = dm.run_calc(mssm_point, thread_id);
			if(result != "Error")
				thread_file << m0 << "," << m12 << "," << result << endl;
			else
				thread_file << m0 << "," << m12 << "," << fail_key << endl;			
		}			
		if(param_space == "cMSSM")
		{
			mssm_point = {m0, m12, a0, tanb_cmssm, sgnmu};
			result = dm.run_calc(mssm_point, thread_id);
			if(result != "Error")			
				thread_file << m0 << "," << m12 << "," << a0 << "," << tanb_cmssm << "," << sgnmu << "," << result << endl;
			else
				thread_file << m0 << "," << m12 << "," << a0 << "," << tanb_cmssm << "," << sgnmu << "," << fail_key << endl;
		}
		if(param_space == "pMSSM")
		{	
			mssm_point = {m1, m2, m3, mmu, ma, at, ab, ata, mel, mtal, mer, mtar, mq1l, mq3l, mur, mtr, mdr, mbr, tanb_pmssm};
			result = dm.run_calc(mssm_point, thread_id);
			if(result != "Error")
				thread_file << m1 << "," << m2 << "," << m3 << "," << mmu << "," << ma << "," << at << "," << ab << "," << ata << "," << 
					mel << "," << mtal << "," << mer << "," << mtar << "," << mq1l << "," << mq3l << "," << mur << "," << mtr << "," <<
					mdr << "," << mbr << "," << tanb_pmssm << "," << result << endl;
			else
				thread_file << m1 << "," << m2 << "," << m3 << "," << mmu << "," << ma << "," << at << "," << ab << "," << ata << "," << 
					mel << "," << mtal << "," << mer << "," << mtar << "," << mq1l << "," << mq3l << "," << mur << "," << mtr << "," <<
					mdr << "," << mbr << "," << tanb_pmssm << "," << fail_key << endl;
		}		
	}
	// Combine all thread files into one node file.
	ofstream writefile;
        char hostname[HOST_NAME_MAX];
        // Node and time uniquely identify an output file, can combine all files later using a python script so that names do not matter.
        gethostname(hostname, HOST_NAME_MAX);
        int clock = time(0);
	writefile.open(finalOutputDir + "/" +  hostname + "-" + to_string(clock) + ".txt", ios_base::app);
        cout << "Saving results to" <<finalOutputDir<< endl;
	for(int i=0; i<number_threads; i++)
	{
		string line;
		ifstream readfile(temporaryDir + "/outputjh" + to_string(i) + ".txt");		
		while(getline(readfile, line))
			writefile << line << endl; 
		readfile.close();
		string cmd = string("rm " + temporaryDir + "/outputjh" + to_string(i) + ".txt"); 
		system(cmd.c_str());
	}	
	writefile.close();
	int stop_s=time(NULL);
	cout << "time per point: " << double(stop_s-start_s) / num_to_generate << " seconds" << endl;
}
