#ifndef Read_hpp
#define Read_hpp

#include <vector>
#include <utility>
#include <unordered_map>
#include <queue>
#include <mutex>
#include<iostream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/accumulators/statistics.hpp>
using namespace boost::accumulators;
using namespace std;


template <class T>
class Class;

class Read {
private:
public:
  typedef pair<double, Class<int>* > score;
  typedef priority_queue<score, vector<score>, less<score> > pqueue;
  /**
   * Instantiates a Read instance.
   * @param  finalString denotes the string to run the jellyfish command on.
   * @param  the reads's own sequence.
   */
  
  Read(string finalString,string sequences);
  ~Read();
  void print_map(std::unordered_map<int, int> const &m);
  int getMapKey(string kmer_s);
  void read_data(FILE *stream);
  void processData(int i);
  unordered_map<int,int>& getKmerCounts();
  void loadSequences();
  void deleteFinalString();
  void computeClassificationNumerator(Class<int>* cls);
  size_t countMapSize();
  string charAt(int pos);
  string getSequences();
  pqueue getConfidences();
  string getFinalString();
  void obtainKmerCountsByInput(std::unordered_map<int, int> &m);
  void unload();
  int getSize();
  score getMaximum();
  static bool STORE_ALL_NUMERATORS;
  unordered_map<int,int>* getKmerCount();
protected:
  bool kmersLoaded = false, sequenceLoaded = false;
  pqueue numerator;
  string finalString,sequences;
  string *sequence;
  mutex numeratorAccess, loadKmersLock, loadSequenceLock;
  unordered_map<int, int> *kmer_counts;
  long long int max(long long int a,
                    long long int b,
                    long long int c,
                    long long int d);
  double maximumNumerator;
  Class<int>* maximumNumeratorClass = NULL;
};

#include "Class.hpp"
#include "NB.hpp"
#endif /* Read_hpp */
