
#ifndef bpe_h
#define bpe_h

#include "basic.h"
#include <libcdsBasics.h>
#include <BitSequence.h>
#include "factorization.h"


#define REPAIR_HDR 7

using namespace std;
using namespace cds_utils;
using namespace cds_static;
//class BPE 
class BPE: public BitSequence {
private: 
  //void save(ofstream * o);
  //gonzalo's
  int free_building_sructs();
  void prepare();
  void createRule(ulong x, ulong y,int FREQ);
  int ini_repair(); 
  int ini_repair(ulong limFreq); 
  void ini_compact();
  ulong huffmanSizeOfC(bool full);
  ulong DeltaSizeOfC();
  
  ulong *s_lengths; //length of each rule, comprimida
  ulong *s_ones;    //number of 1s, idem.
  
  ulong *gn_symb;
  ulong *gn_ones;
  ulong *gn_off;
  ulong factor;
  ulong integers;
  
  ulong b;
  ulong s; //sample period in C.
  
  ulong *symbols; //here is the compressed data
  ulong m; //length of symbols
  BitSequenceRG *BRR; //bitmap to indicate the structure of symbols_new
  ulong *symbols_new; // table to decompress de new codes
  ulong  symbols_new_len; // len of symbols_new
  //ulong n; //We use length, inherited. 
  ulong nbits;
  
  //WHEN SOME_RULES==TRUE
  ulong bits_r;//bits of rules lenghts
  ulong bits_o;//bits of ones 
  ulong shift;
  ulong max_value;
/////////////////////////////////////////////////////////////////////////////////////////////////7
  ulong shift_it();
  void new_value(ulong *symbols_pair,ulong *symbols_new_value,ulong *k1,ulong *j,ulong pos, int *max_prof, ulong *symbols_new_prof);
  int compress_pair_table();
  void compress_final_seq();
  void save_lengths();
  void ini_bes();
  ulong max_prof;
  ulong largo_max_real;
  ulong max_ones_real;
  ulong largo_max; //param
  ulong max_ones; //param
  int * b_1;
  int * b_2;
  ulong max_bt;
  int cod_rules; 
  void BuildRank();
  virtual bool C_ARB_GNS_access(uint i) const; //max length rules
  virtual bool C_ARB_GNS_SOME_access(uint i) const;//some rules
  virtual bool C_DAC_access(uint i) const;//some rules
  
  //virtual uint C_ARBrank(uint i);
  virtual uint C_ARB_GNS_rank(uint i) const;
  virtual uint C_ARB_GNS_SOME_rank(uint i) const;
  virtual uint C_DAC_rank(uint i) const;
  ulong * C_dispairall(); 
  int T;
  bool SOME_rules;
  bool MAXLEN_rules;
  bool DAC_rules;
  BitSequenceRG *B_rules; 
  ulong *v_lengths;
  ulong *v_ones;
  ulong n_rules;
  factorization * dac_lengths;
  factorization * dac_ones;
public: 
  ulong max_assigned;
  void save(ofstream & of) const;
  static BPE * load(ifstream & f) ;
  uint getSize() const;
  ulong maxbits_rules();
  void re_set_max_rule(ulong tbits);
  void re_set_max_rule_len(ulong newmaxlen);
  void re_set_SOME_rules2(ulong param1,ulong param2); //nuevo metodo
  void re_set_DAC_for_rules(ulong param1); 
  void re_build_samples(ulong tbits);
  void build_GN_samples(ulong period);
  BPE * getNext();
  void complete();
  void MaxBT(); 
  ulong getMaxProf();
  uint getLargoMaxReal();
  BPE(ulong *a, ulong n, bool _verbose, int rules, int T, int limfreq);//DV
  BPE(ulong *a, ulong n, bool _verbose, int rules, int T);//DV
  BPE();
  ~BPE(); //destructor
  ulong * dispairall(); 
  void cotas();
  
  uint size() const;
  uint sizePlain() const;
  uint sizeRulesLength() const;
  uint sizeSamples() const;
  uint DictSize();
  uint minSize();
  uint huffSize();
  uint DeltaSize();
  uint DACSize();
  uint DACVARSize();
  double entropyOfC();
  /*load-save functions*/
  int save(FILE *f);
  int load(FILE *f);
	static BPE * my_load(FILE * fp);
  
  BPE(FILE *f, int *error); 
  
  ulong GetSymbol(int i);
  uint getM();
  ulong getNbits();
  virtual bool access(const size_t i) const ;
  virtual size_t rank1(const size_t  i) const;
};

#endif
