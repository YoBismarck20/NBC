#include <mutex>
#include<iostream>
#include<string>
#include<stdio.h>
#include<fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>



using namespace std;



void print_map(std::unordered_map<int, int> const &m)
{
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

int getMapKey(string kmer_s){
	int kmer;
	kmer = 0;
	kmer_s.erase(std::remove(kmer_s.begin(), kmer_s.end(), '\n'), kmer_s.end());
	for (string::iterator it = kmer_s.begin(); it != kmer_s.end(); it++){
		kmer<<=2;
		switch(*it){
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

unordered_map<int, int>& getComputedKmers(FILE* fp){
	char buf[1234];
	unordered_map<int, int> *kmer_counts;
	kmer_counts = new unordered_map<int, int>;
	char count[30];
	int key=0;
	int val=0;
	int num = 0;
	while (fgets(buf,sizeof(buf),fp)){
		if(*buf =='>'){
			cout << "Integer"<<"\n";
			string str(buf);
			string sval = str.substr(1);
			val=stoi(sval);
		}else if (*buf == '#') {
			kmer_counts = new unordered_map<int,int>;
			num += 1;
			cout << "Number ID: " << num << "\n";
			cout << "end of a Read" << "\n" ;
		} else {
			cout << buf << "\n";
			key=getMapKey(buf);
			(*kmer_counts)[key]=val;
		}
	}
	return *kmer_counts;

}

int main(int argc, char* argv[]){
	unordered_map<int, int> *kmer_counts;
	FILE *fp;
	fp = popen("parallel -k -j 16 < outputv2","r");
	*kmer_counts = getComputedKmers(fp);
	print_map(*kmer_counts);
	/*
	while (fgets(buf,sizeof(buf),fp)){
		if(*buf =='>'){
			cout << "Integer"<<"\n";
			string str(buf);
			string sval = str.substr(1);
			val=stoi(sval);
		}else if (*buf == '#') {
			print_map(*kmer_counts);
			kmer_counts = new unordered_map<int,int>;
			cout << "end of a Read" << "\n" ;
		} else {
			cout << buf << "\n";
			key=getMapKey(buf);
			(*kmer_counts)[key]=val;
		}
	}
	*/
}

