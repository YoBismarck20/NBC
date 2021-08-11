//
//  main.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//
#include "Read.hpp"
#include "NB.hpp"
#include "Diskutil.hpp"
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <unordered_map>
#include <boost/program_options.hpp>
using namespace std;
namespace p_opt = boost::program_options;

string usageMsg("./NB.run [mode: train/classify/benchmark] [source dir] [options]\n");
const int CONF_THRESHOLD = 100000;
const int DEF_KMER_SIZE = 6;
const string GEN_SEQ_EXTEN(".fna");
const string KMER_EXTEN(".kmr");
const string PROG_VER("NB v. 0.1.5a-dev.");
const string DEF_SAVEDIR("./NB_save");
unordered_map<string, vector<double>* > confidence_list;
vector<path> read_filenames;
string resultfile("");


