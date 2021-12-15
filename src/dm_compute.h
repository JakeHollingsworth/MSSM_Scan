#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <ctime>
#include <random>
#include <algorithm>
#include <omp.h>
#include <memory>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

class DM
{
public:
    string mo_dir, ss_dir, temp_dir; 
    DM(string micro_omega_dir,string softSUSY_dir, string temporary_dir);
    vector<map<string, double>> get_running_params(string x);
    map<string, double> get_weak_params(string x);
    string getError(string tempDir, string threadID);
    string run_calc(vector<double> input_point, string thread_id);
    vector<double> get_gut_scale_gauge_couplings(string x, vector<map<string,double>> running_params);
    string format_string(string x);
    vector<string> phypar = {"Au", "Ac", "At", "Ad", "As", "Ab", "Ae", "Amu", "Atau", "beta", "g", "gprime", "g3", "M1", "M2", "M3", "mA2", "mu", "mH12", "mH22", "muR", "mcR", "mtR", "mdR", "msR", "mbR", "meR", "mmuR", "mtauR", "mqL1", "mqL2", "mqL3", "meL", "mmuL", "mtauL", "vev", "Yt", "Yb", "Ytau"};
    // A0 is pseudoscalar mass, not univ. trilinear coupling.
    vector<string> weakpar = {"MW", "h0", "H0", "A0", "H+", "~g", "~neutralino(1)", "~neutralino(2)", "~neutralino(3)", "~neutralino(4)", "~chargino(1)", "~chargino(2)", "~dL", "~uL", "~sL", "~cL", "~b1", "~t1", "~eL", "~nueL", "~muL", "~numuL", "~stau1", "~nutauL", "~dR", "~uR", "~sR", "~cR", "~b2", "~t2", "~eR", "~muR", "~stau2", "N{1,1}", "N{1,2}", "N{1,3}", "N{1,4}"};
};
