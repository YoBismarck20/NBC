#include <iostream>
#include <typeinfo>
#include <fstream>
#include <string>
#include <unordered_map>
using namespace std;

void print_map(std::unordered_map<int, int> const &m)
{
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}
int getMapKey(string kmer_s) {
    int kmer;
    kmer = 0;
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
unordered_map<int, int>& read_data(FILE *stream) {
    unordered_map<int,int> *kmer_counts;
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
            cout << val << "\n";
        }
        else {
            key=getMapKey(buf);
            (*kmer_counts)[key]=val;
            cout << buf;
        }
    }
    return *kmer_counts;
}
int main(int argc, char *argv[]) {
    unordered_map<int,int> *kmer_counts;
    if(argc==1) {
        printf("\nNeed only one argument: the file path for now");
        return -1;
    }
    fstream newFile;
    newFile.open(argv[1], ios::in);
    if(newFile.is_open()) {
        printf("File is Open\n");
        string s;
        string header;
        string line;
        string finalString;
        int iterator;
        iterator=0;
        while(getline(newFile,line)) {
            if(iterator > 2) {
                break;
            }
            if(line.at(0) == '>' && !s.empty()) {
                finalString=header+s;
                kmer_counts = new unordered_map<int, int>;

                cout << header;
//				system(("echo -e \"" + finalString+"\" | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin").c_str());
                FILE *output;
                output=popen(("echo -e \"" + finalString+"\" | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin").c_str(),"r");
                *kmer_counts=read_data(output);
                print_map(*kmer_counts);
                pclose(output);
                header = line;
                header.append("\n");
                s = "";
                iterator++;
            }
            else if (line.at(0) =='>' && s.empty()) {
                header=line;
                header.append("\n");
            }
            else {
                s.append(line);
            }
        }
        newFile.close();

    }
    return 0;
}
