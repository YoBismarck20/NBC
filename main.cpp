//
//  main.cpp
//  NB_C++
//
//  Created by Alexandru Cristian on 06/04/2017.
//
//
#include "Read.hpp"
#include <chrono>
#include "NB.hpp"
#include "Diskutil.hpp"
#include <cstring>
#include <cstdlib>
#include <cstdint>
#include <unordered_map>
#include <boost/program_options.hpp>
#include <thread>
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

int getMapKey(string kmer_s) {
    int kmer;
    kmer = 0;
    kmer_s.erase(std::remove(kmer_s.begin(), kmer_s.end(), '\n'), kmer_s.end());
    for (string::iterator it = kmer_s.begin(); it != kmer_s.end(); it++) {
        kmer<<=2;
        switch(*it) {
        case 'A':
        case 'a':
            kmer+=0;
            break;
        case 'C':
        case 'c':
            kmer+=1;
            break;
        case 'G':
        case 'g':
            kmer+=2;
            break;
        case 'T':
        case 't':
            kmer+=3;
            break;
        }
    }
    return kmer;
}

void trainNB(NB &nb, path srcdir, string extension, unsigned int nbatch,
             uint64_t memoryLimit) {
    unsigned int count = 1, counter = 0;
    uint64_t usedMemory = 0;
    vector<tuple<string, path, path> > result =
        Diskutil::getTrainingGenomePaths(srcdir, extension);
    string cls_s="-1";
    Class<int> *current = NULL;
    for(vector<tuple<string, path, path> >::iterator iter = result.begin();
            iter != result.end(); iter++, counter++) {
        unsigned int savefileSize = 0;
        bool loadedNewClass = false;

        // FASTA files not needed for training, so just add up the kmr file size
        int genomeSize = Diskutil::getFileSize(get<1>(*iter));

        if(cls_s.compare(get<0>(*iter)) != 0) {

            cls_s = get<0>(*iter);

            current = nb.getClass(cls_s);

            if(current == NULL) {
                path save_file = path(nb.getSavedir().native()
                                      + path::preferred_separator
                                      + cls_s + "-save.dat");
                current = new Class<int>(cls_s,
                                         nb.getKmerSize(),
                                         save_file);
                nb.addClass(current);
            }

            savefileSize = Diskutil::getFileSize(current->getSavefilePath());
            loadedNewClass = true;
        }

        if(memoryLimit != 0 && usedMemory + genomeSize + savefileSize> memoryLimit) {
            nb.processClassUpdates();
            usedMemory = 0;
            cls_s = "-1"; // This will force the next iteration to add this savefile's size again
        }

        Genome *genome = new Genome(get<1>(*iter), get<2>(*iter));
        current->queueGenome(genome);
        if(loadedNewClass) {
            nb.addClassToUpdateQueue(current);
            usedMemory += savefileSize;
        }
        usedMemory += genomeSize;

        if(nbatch != 0 && counter % nbatch == 0) {
            nb.processClassUpdates();
        }
    }

    nb.processClassUpdates();
}


void printConfidences() {
    if(resultfile.empty()) {
        return;
    }

    std::ofstream out(resultfile);
    int nclasses;
    nclasses = confidence_list.begin()->second->size();
    out<<"Filename,";
    for(unordered_map<string, vector<double>* >::iterator it = confidence_list.begin(); it != confidence_list.end(); it++) {
        out<<it->first<<",";
    }
    out<<"\n";
    for(int count = 0; count < nclasses; count++) {
        out<<read_filenames[count]<<",";
        for(unordered_map<string, vector<double>* >::iterator it = confidence_list.begin(); it != confidence_list.end(); it++) {
            out<<(*(it->second))[count]<<",";
        }
        out<<"\n";
    }
    out.close();
}

void addToPosteriorList(Genome::pqueue &queue) {
    double confidence;
    string pred_class;

    while(!queue.empty()) {
        pred_class = queue.top().second->getId();
        confidence = queue.top().first;
        queue.pop();

        if (confidence_list.find(pred_class) == confidence_list.end()) {
            confidence_list[pred_class] = new vector<double>;
        }

        confidence_list[pred_class]->push_back(confidence);
    }
}

int printClassifierResults(vector<Read*> reads,
                           vector<tuple<string, path, path> > result) {
//Need to ask how this algorithm works
    unsigned int correct=0, total = reads.size();

    for(int i=0; i < total; i++) {
        string pred_class;
        double posterior, prior;
        Read::pqueue queue;
        if (Read::STORE_ALL_NUMERATORS) {
            queue = reads[i]->getConfidences();

            pred_class = queue.top().second->getId();
            posterior = queue.top().first;
        } else {
            Read::score score = reads[i]->getMaximum();

            pred_class = score.second->getId();
            prior = score.first;
        }

        cout<<"Genome with class "<<get<0>(result[i]);
        cout<<", predicted "<<pred_class;
        if (Read::STORE_ALL_NUMERATORS) {
            addToPosteriorList(queue);
        }
        cout<<'\n';
        cout.flush();

        if(pred_class.compare(get<0>(result[i])) == 0) {
            correct++;
        }
        cout << "Correct is properly outputted" <<"\n";
        /*
        int position=1;
        while(!queue.empty() && pred_class.compare(get<0>(result[i])) != 0){
          queue.pop(); position++;
          pred_class = queue.top().second->getId();
        }

        if(pred_class.compare(get<0>(result[i])) != 0){
          cout<<"[ERROR] Actual class not in queue.\n";
        }else{
          cout<<"Actual class was on position "<<position<<", score "<<queue.top().first<<"\n";
        }
        */
    }
    //printConfidences();

    for(vector<Read*>::iterator iter = reads.begin();
            iter != reads.end(); iter++) {
        delete *iter;
    }
    cout <<"Deleted everything" <<"\n";

    return correct;
}

void processReadData(Read* read, int i) {
    read->processData(i);
}
void classifyNB(NB &nb, path srcdir, string extension, unsigned int nbatch,
                uint64_t memoryLimit, unsigned int nthread) {
    vector<Read*> reads;
    // Variables for computing kmer counts
    unordered_map<int, int> *kmer_counts;
    kmer_counts = new unordered_map<int, int>;
    char buf[1234];
    char count[30];
    int key=0;
    int val=0;
    int sizeOfTemp=0;
    FILE *fp;

    // end of varaibles for kmer counts
    unsigned int correct=0, counter=0, total=0;
    uint64_t usedMemory = 0;
    string filename;
    int sizeOfStandby;
    vector<tuple<string, path, path> > result =
        Diskutil::getTrainingGenomePaths(srcdir, extension);
    vector<Read*> standby;
    vector<thread*> processor(nthread);

    std::ofstream out("outputv2");
    for(vector<tuple<string, path, path> >::iterator iter =
                result.begin(); iter != result.end(); iter++, counter++) {
        auto start2 = std::chrono::high_resolution_clock::now();

        filename=get<2>(*iter).string();
        cout << filename << "\n";
        std::fstream newFile;
        newFile.open(filename, ios::in);
        //delete this ASAP
        if(newFile.is_open()) {
            printf("File is Open\n");
            string s;
            string header;
            string line;
            string finalString;
            int iterator;
            iterator=0;

            while(getline(newFile,line)) {
                if (line.at(0) =='>' && s.empty()) {

                    header=line;
                    header.append("\\n");
                }
                else {
                    //Due to testing for Alex's Code; need to modify this clause
                    //Once done; replace evertyhing with s.append(line);
                    s.append(line);
                    finalString=header+s;
                    while(finalString.find("\'") != std::string::npos) {
                        finalString.replace(finalString.find("\'"),1,"\"");
                    }
                    //Read *read= new Read(finalString,s);
                    cout << "This is the final string" << "\n";
                    cout << finalString << "\n";
                    //Test Code; remove ASAP

                    if(memoryLimit != 0 && usedMemory > memoryLimit) {

                        auto start = std::chrono::high_resolution_clock::now();
                        auto stop2 = std::chrono::high_resolution_clock::now();
                        auto duration2 = std::chrono::duration_cast<std::chrono::minutes>(stop2 - start2);
                        cout << "This is the time it took to load in this batch: " << duration2.count() << "\n";
                        nb.classify(reads);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
                        cout << "This is the time it took for this batch: " << duration.count() << "\n";
                        correct += printClassifierResults(reads, result);
                        usedMemory = 0;

                        iter = result.erase(result.begin(), result.begin() + reads.size());
                        reads.clear();
                        auto start2 = std::chrono::high_resolution_clock::now();

                    }
                    //Going to try to batch the read and then push back onlt after if it reaches the max amunt of threads
                    //
                    if(sizeOfTemp < 10000) {
                        // No matter how much things there is, it will always be in the file...
                        //standby.push_back(read);
                        sizeOfTemp += 1;
                        cout << "Size of tempFile" << sizeOfTemp << "\n";
                        string target = "echo -e \'" + finalString+"\' | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin";
                        out<< target << "\n";
                        out << "echo \"#\"\n";
                        out << "echo @" + filename + "\n";
                        //out << "echo $" <<  "\n";
                    }
                    else {
                        //remove this afer finishing getting the strings
                        //This is the part where it doesnt do it for all....
                        out.close();
                        fp = popen("parallel -k -j 16 < outputv2","r");
                        sizeOfTemp = 0;
                        while (fgets(buf,sizeof(buf),fp)) {
                            if(*buf =='>') {
                                string str(buf);
                                string sval = str.substr(1);
                                cout <<  "Integer: "<<sval<<"\n";
                                val=stoi(sval);
                            } else if (*buf == '#') {
                                Read *read= new Read(finalString,s);
                                cout << "Before getting Kmer counts from current map" << "\n";
                                read->obtainKmerCountsByInput(*kmer_counts);
                                cout << "After getting counts"  << "\n";
                                kmer_counts = new unordered_map<int,int>;
                                cout << "end of a Read" << "\n" ;
                                reads.push_back(read);
                                usedMemory += read->countMapSize();
                            } else {
                                cout << buf << "\n";
                                key=getMapKey(buf);
                                (*kmer_counts)[key]=val;
                            }
                        }
                        //exit(0);
                        cout << "Starting batch processing" << "\n";
                        //standby.clear();
                        //standby.push_back(read);
                        sizeOfTemp = 0;
                    }
//  	  reads.push_back(read);
//	  usedMemory += stringMemory2;
                    header = line;
                    header.append("\n");
                    s = "";
                    iterator++;
                    total++;
                }
            }
            newFile.close();
        }
    }
    string finalString = "standby";
    string s = "standby2";
    int debugOnlyReadCount = 0;
    fp = popen("parallel -k -j 16 < outputv2","r");
    while (fgets(buf,sizeof(buf),fp) != NULL ) {
        if(*buf =='>') {
            string str(buf);
            string sval = str.substr(1);
            cout <<  "Integer: "<<sval<<"\n";
            val=stoi(sval);
        } else if (*buf == '#') {
            Read *read= new Read(finalString,s);
            cout << "Before getting Kmer counts from current map" << "\n";
            read->obtainKmerCountsByInput(*kmer_counts);
            cout << "After getting counts"  << "\n";
            kmer_counts = new unordered_map<int,int>;
            cout << "end of a Read" << "\n" ;
            reads.push_back(read);

            debugOnlyReadCount += 1;
            cout <<"Number of reads read in: "<< debugOnlyReadCount << "\n";
            usedMemory += read->countMapSize();

        } else if (*buf == '@') {
            //Debug ONLY
            cout << "This is the file read in" << "\n";
            cout << buf << "\n";
        } else if (*buf == '$') {
            cout << "This is debugging" <<"\n";


        } else {
            cout << buf << "\n";
            key=getMapKey(buf);
            (*kmer_counts)[key]=val;
        }
    }
    int error;
    error = pclose(fp);
    cout << "This is the potential error "<< error << "\n";
    cout <<"Size of reads: " << reads.size() << "\n";
    auto start = std::chrono::high_resolution_clock::now();
    nb.classify(reads);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
    cout << "This is the time it took for this batch: " << duration.count() << "\n";

    correct += printClassifierResults(reads, result);
    //iter = result.erase(result.begin(), result.begin() + reads.size());
    reads.clear();

    cout<<"Accuracy: "<<correct*1.0/total<<"\n";
    cout<<"Total correct: "<<correct<<"\n";
    cout<<"Total classified: "<<total<<"\n";
}


int main(int argc, char* argv[]) {
    unsigned int nbatch, nthreads, kmersize;
    uint64_t memLimit;
    string kmer_ext, srcdir, mode, savedir;
    bool print_posterior;

    p_opt::options_description generic("Generic options");
    generic.add_options()
    ("help,h", "Print help message")
    ("version,v", "Print version information");

    p_opt::options_description hidden("Hidden options");
    hidden.add_options()
    ("mode", "Sets mode of program, train or classify.")
    ("srcdir", "Path to source folder");

    p_opt::options_description visible("Allowed options");
    visible.add_options()
    ("savedir,s", p_opt::value<string>(&savedir)->default_value(DEF_SAVEDIR),
     "Path to save folder")
    ("kmersize,k", p_opt::value<unsigned int>(&kmersize)->default_value(DEF_KMER_SIZE),
     "Kmer size used in count files")
    ("memlimit,m", p_opt::value<uint64_t>(&memLimit)->default_value(0),
     "Cap memory use to a predefined value (KBs).")
    ("nthreads,t", p_opt::value<unsigned int>(&nthreads)->default_value(1),
     "Number of threads to spawn, 1 by default")
    ("ext,e", p_opt::value<string>(&kmer_ext)->default_value(KMER_EXTEN),
     "Extension of kmer count files, \".kmr\" by default")
    ("nbatch,n", p_opt::value<unsigned int>(&nbatch)->default_value(0),
     "Number of genomes to load at one time in memory, \
all at once by default")
    ("p_posterior,p", p_opt::value<bool>(&print_posterior)->default_value(Genome::STORE_ALL_NUMERATORS),
     "Print posteriors for every classified read. This \
flag increases the classifier's memory usage and is not compatible with the memory cap flag.")
    ("resf,r", p_opt::value<string>(&resultfile)->default_value(""),
     "Output path for posterior logs. Log files will not be created by default.");

    p_opt::positional_options_description pos_args;
    pos_args.add("mode", 1);
    pos_args.add("srcdir", 1);

    p_opt::options_description cmdline_options;
    cmdline_options.add(generic).add(visible).add(hidden);

    p_opt::variables_map opt_map;
    p_opt::store(
        p_opt::command_line_parser(argc, argv).options(cmdline_options)
        .positional(pos_args).run(),
        opt_map);
    p_opt::notify(opt_map);

    if(opt_map.count("version")) {
        cout<<PROG_VER<<"\n";
        return 1;
    }

    if(opt_map.count("help") || opt_map.count("mode") == 0
            || opt_map.count("srcdir") == 0) {
        cout<<usageMsg<<"\n"<<generic<<"\n"<<visible<<"\n";
        return 1;
    }

    srcdir = opt_map["srcdir"].as<string>();
    mode = opt_map["mode"].as<string>();

    create_directories(savedir);
    NB nb(kmersize, path(savedir), nthreads);
    nb.debug_flag = NB::Debug::LOG_SOME;

    if (!resultfile.empty()) {
        print_posterior = true;
    }

    Genome::STORE_ALL_NUMERATORS = print_posterior;

    if(mode.compare("train") == 0) {
        cout<<"Train mode.\n";

        trainNB(nb, path(srcdir), kmer_ext, nbatch, memLimit);

        cout<<"Training complete.\n";

    } else if(mode.compare("classify") == 0) {
        cout<<"Classify mode.\n";

        classifyNB(nb, path(srcdir), kmer_ext, nbatch, memLimit, nthreads);

    } else if(mode.compare("benchmark") == 0) {
        trainNB(nb, path(srcdir), kmer_ext, nbatch, memLimit);

        classifyNB(nb, path(srcdir+"_test"), kmer_ext, nbatch, memLimit,nthreads);
    } else {
        cout<<usageMsg<<"\n"<<generic<<"\n"<<visible<<"\n";
        return 1;
    }
    return 0;
}

