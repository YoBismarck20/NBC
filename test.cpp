#include <iostream>
#include <typeinfo>
#include <fstream>
#include <string>
#include <unordered_map>
#include "Read.hpp"
#include "Genome.hpp"
#include <boost/program_options.hpp>
#include "NB.hpp"
#include "Diskutil.hpp"
#include <boost/algorithm/string/replace.hpp>
void replaceAll(std::string& source, const std::string& from, const std::string& to)
{
    std::string newString;
    newString.reserve(source.length());  // avoids a few memory allocations

    std::string::size_type lastPos = 0;
    std::string::size_type findPos;

    while(std::string::npos != (findPos = source.find(from, lastPos)))
    {
        newString.append(source, lastPos, findPos - lastPos);
        newString += to;
        lastPos = findPos + from.length();
    }

    // Care for the rest after last occurrence
    newString += source.substr(lastPos);

    source.swap(newString);
}


int main(int argc, char *argv[]) {
    if(argc==1) {
        printf("\nNeed only one argument: the file path for now");
        return -1;
    }
    else if (argc  == 4) {
        printf("\nNeed only 2 argumenet: the file path AND the kmer count  path");
    }
    std::fstream newFile;
    newFile.open(argv[1], ios::in);
    vector<Read*> reads;
    uint64_t usedMemory = 0;
    uint64_t memoryLimit = 1000;
    if(newFile.is_open()) {
        printf("File is Open\n");
        string s;
        string header;
        string line;
        string finalString;
        int iterator;
        iterator=0;
        while(getline(newFile,line)) {
            if(line.at(0) == '>' && !s.empty()) {
                finalString=header+s;
                Read *read= new Read(finalString,s);

                cout << "This is what was read:" << finalString <<"\n";

                read->processData();
//				system(("echo -e \"" + finalString+"\" | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin").c_str());
                int resultantMemory = read->getSize();
                usedMemory += (uint64_t) resultantMemory;
                cout << usedMemory << "\n";
                read->processData();
                reads.push_back(read);
                if(memoryLimit != 0 && usedMemory > memoryLimit) {
                    cout << "LIMIT REACHED" <<"\n";
                    usedMemory=0;
                }
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
                /* Remove everything below it */
                finalString=header+s;

                while(finalString.find("\'") != std::string::npos) {
                    finalString.replace(finalString.find("\'"),1,"\"");
                }

                Read *read= new Read(finalString,s);

                cout << "This is what was read:" << finalString <<"\n";

                read->processData();
//				system(("echo -e \"" + finalString+"\" | ./jellyfish-linux count -m 15 -s 10M /dev/stdin -o /dev/stdout | ./jellyfish-linux dump /dev/stdin").c_str());
                int resultantMemory = read->getSize();
                usedMemory += (uint64_t) resultantMemory;
                read->processData();
                read->print_map(read->getKmerCounts());
                cout << usedMemory << "\n";
                reads.push_back(read);
                if(memoryLimit != 0 && usedMemory > memoryLimit) {
                    cout << "LIMIT REACHED" <<"\n";
                    usedMemory=0;
                }
                header = line;
                header.append("\n");
                s = "";

            }
        }

        cout <<"This is the amount of memory needed: " << usedMemory << "\n";
        newFile.close();

    }
    else {


        cout <<"Not open " << "\n";
    }
    cout << "This is the Genome's Version" << "\n";
    Genome *genome = new Genome (argv[2],argv[1]);
    genome->loadKmerCounts();
    genome->print_map(genome->getKmerCounts());
    return 0;
}
