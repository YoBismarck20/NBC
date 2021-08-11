//
///
//
// Reads.cpp
// NB_C++
//
//
//
//
//

#include "Read.hpp"
bool Read::STORE_ALL_NUMERATORS = false;

Read::Read(string _finalString, string _sequences) {
    finalString = _finalString;
    sequences = _sequences;
}

Read::~Read() {
    unload();
}

void Read::print_map(std::unordered_map<int, int> const &m)
{
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

size_t Read::countMapSize() {
    size_t count = 0;
    for (unsigned i = 0; i < getKmerCounts().bucket_count(); ++i) {
        size_t bucket_size = getKmerCounts().bucket_size(i);
        if (bucket_size == 0) {
            count++;
        }
        else {
            count += bucket_size;
        }
    }
    return count;
}
int Read::getMapKey(string kmer_s) {
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

void Read::read_data(FILE *stream) {
    char buf[1234];
    char count[30];
    int key=0;
    int val=0;
    kmer_counts = new unordered_map<int, int>;
    while (fgets(buf,sizeof(buf),stream)) {
        if(*buf =='>') {
            string str(buf);
            string sval = str.substr(1);
            val=stoi(sval);
        }
        else {
            key=getMapKey(buf);
            (*kmer_counts)[key]=val;
        }
    }
}

void Read::processData(int i) {
    FILE *output;
    string theI = to_string(i);
    //output=popen(("echo -e \'" + finalString+"\' | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin").c_str(),"r");
    string cmd1 = "echo -e \'" + finalString+"\' > testfile" + theI;
    string cmd2 = "./jellyfish-linux count -m 15 -s 10M testfile"+ theI + " -o testfiledump"+theI;
    string cmd3 = "./jellyfish-linux dump testfiledump"+theI + "> testfileout" + theI;
    std::system(cmd1.c_str());
    std::system(cmd2.c_str());
    std::system(cmd3.c_str());
    //cout << "Pipe " << theI << "has finished jellfish" << "\n";
    output=popen(("cat testfileout"+theI).c_str(),"r");
    //output=popen(("echo -e \'" + finalString+"\' > file" + theI + " | ./jellyfish-linux count -m 15 -s 10M file"+ theI + " -o fileout"+theI + " | ./jellyfish-linux dump fileout"+theI + "> filedump"+ theI+ "| cat filedump"+theI).c_str(),"r");
    read_data(output);
    pclose(output);
    deleteFinalString();
    kmersLoaded=true;
}

void Read::obtainKmerCountsByInput(std::unordered_map<int, int> &m) {
    print_map(m);
    kmer_counts = &m;

}

unordered_map<int,int>& Read::getKmerCounts() {
    /*
    if (!kmersLoaded){
    	processData();
    }*/
    //For testing
    return *kmer_counts;
}
string Read::getFinalString() {
    return finalString;
}
string Read::getSequences() {
    return sequences;
}



void Read::deleteFinalString() {
    finalString.clear();
}
int Read::getSize() {
    return getSequences().size() + getFinalString().size(); //+ countMapSize();
}
unordered_map<int,int>* Read::getKmerCount() {
    return kmer_counts;
}
void Read::unload() {
    if(kmersLoaded) {
        delete kmer_counts;
        kmersLoaded=false;
    }
}

void Read::computeClassificationNumerator(Class<int>* cl) {
    accumulator_set<double, features<tag::sum_kahan> > current, sumfrq;

    current(cl->getNGenomes_lg());
    ostringstream strs;
    if(NB::debug_flag == NB::Debug::LOG_ALL) {
        strs<<"("<<cl->getId()<<"): "<<sum_kahan(current);
    }

    for(unordered_map<int, int>::iterator freq = getKmerCounts().begin();
            freq != getKmerCounts().end(); freq++) {
        sumfrq(freq->second);
        current(freq->second * cl->getFreqCount_lg(freq->first));
        if(NB::debug_flag == NB::Debug::LOG_ALL) {
            strs<<" + "<<freq->second<<" * "<<cl->getFreqCount_lg(freq->first);
        }
    }
    current(-sum_kahan(sumfrq) * cl->getSumFreq_lg());
    if(NB::debug_flag == NB::Debug::LOG_ALL) {
        strs<<" - "<<sum_kahan(sumfrq)<<" * "<<cl->getSumFreq_lg()<<" = ";
        strs<<sum_kahan(current)<<"\n";
    }

    numeratorAccess.lock();

    if(!Read::STORE_ALL_NUMERATORS) {
        double candidateNumerator = sum_kahan(current);
        if(maximumNumeratorClass == NULL
                || maximumNumerator < candidateNumerator) {
            maximumNumerator = candidateNumerator;
            maximumNumeratorClass = cl;
        }
    } else {
        numerator.push(make_pair(sum_kahan(current), cl));
        if(NB::debug_flag == NB::Debug::LOG_ALL) {
            cout<<strs.str();
        }
    }

    numeratorAccess.unlock();
}
Read::pqueue Read::getConfidences() {
    accumulator_set<double, features<tag::sum_kahan> >  p_denominator(1);
    double max, denominator;

    pqueue num_cpy = numerator;

    max = num_cpy.top().first;
    vector<score> cache;
    cache.push_back(num_cpy.top());
    num_cpy.pop();

    while(!num_cpy.empty()) {
        p_denominator(exp(num_cpy.top().first - max));
        cache.push_back(num_cpy.top());
        num_cpy.pop();
    }
    denominator = max + log(sum_kahan(p_denominator));

    pqueue confidence;
    for(vector<score>::iterator num=cache.begin();
            num != cache.end(); num++) {
        double confidence_lg = num->first - denominator;
        confidence.push(make_pair(confidence_lg, num->second));
        //processFreqs(getKmerCounts(), num->second, confidence_lg);
    }

    return confidence;
}

Read::score Read::getMaximum() {
    numeratorAccess.lock();
    cout << maximumNumerator << "\n";
    cout <<maximumNumeratorClass <<"\n";
    Read::score sc = make_pair(maximumNumerator, maximumNumeratorClass);
    numeratorAccess.unlock();
    return sc;
}

string Read::charAt(int pos) {
    return string(1, getSequences()[pos]);
}

long long int Read::max(long long int a,
                        long long int b,
                        long long int c,
                        long long int d) {
    if(a>=b && a>=c && a>=d) {
        return a;
    } else if(b>=a && b>=c && b>=d) {
        return b;
    } else if(c>=b && c>=a && c>=d) {
        return c;
    } else {
        return d;
    }
}
