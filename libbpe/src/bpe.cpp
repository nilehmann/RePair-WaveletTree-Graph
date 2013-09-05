/* 
   Bitsequences supporting rank and selected compressed with RePair.
   Daniel Valenzuela.
 */
#include "records.h"
#include "hash.h"
#include "heap.h"

#include "bpe.h"

#include "math.h"
#include "assert.h"
#include "maloc.h"
#include "mylist.h"
#include "huff.h"

#include <sys/times.h>
//rcanovas
#define FACT_RANK 20

#define LENGTH(X) ((X) < 2 ) ? 1   : rule_length[(X)-2] 
#define ONES(X)   ((X) < 2 ) ? (X) : rule_ones[(X)-2] 

/* These variables are used only in construction time */
static ulong *symbols_pair; 
static ulong *data;
static ulong *symbols_new_bit; 
static ulong symbols_new_bit_len;
bool verbose=true;
static bool testing=true;

int minsize = 256; // to avoid many Reallocs at small sizes, should be ok as is

float RPfactor = 0.75; // 1/extra space overhead; set closer to 1 for smaller and slower execution.
int u; // |text| and later current |C| with gaps

int *C; // compressed text

int c;  // real |C|

int alph; // max used terminal symbol

//int n; // |R|
int new_symbol; // |R|

int PRNC = 0;
int PRNR = 0;
int PRNP = 0;
int PRNL = 0;

Tlist *myL; // |myL| = c;

Thash Hash; // hash table of pairs

Theap Heap; // special heap of pairs

Trarray Rec; // records

bool verbose_rules = false;

int max_freq;

mylist *reemplazos;
ulong * rule_length;
ulong * rule_ones;

#define NULL_POS  ((ulong)-1)
#define NULL_DATA ((ulong)-1)


ulong BPE::shift_it(){
  max_freq=0;
  new_symbol=2;
  ulong  min=C[0],max=C[0];
  for(ulong i=0; i<length ; i++) {
    if(min > C[i]) min = C[i];
    if(max < C[i]) max = C[i];
  }
  if(min != 0 ) {
    for(ulong i=0; i<length ; i++) 
      C[i] = C[i]-min;
    max=max-min;
  }
  this->shift=min;
  this->max_value=max;
  this->max_assigned=max;
  return 0;
}
//GN
void BPE::prepare()
{
  { int i,id;
    Tpair pair;
    c = u = length;
    alph = 0;
    for (i=0;i<u;i++) 
	{ if (C[i] > alph) alph = C[i];
	}
    new_symbol = ++alph;
    Rec = createRecords(RPfactor,minsize);
    Heap = createHeap(u,&Rec,RPfactor,minsize);
    Hash = createHash(256*256,&Rec);
    myL = (Tlist*)Malloc(u*sizeof(Tlist));
    assocRecords (&Rec,&Hash,&Heap,myL);
    for (i=0;i<c-1;i++) 
	{ pair.left = C[i]; pair.right = C[i+1];
	  id = searchHash (Hash,pair);
    	  if (id == -1) // new pair, insert
	     { id = insertRecord (&Rec,pair);
	       myL[i].next = -1;
	     }
	  else 
	     { myL[i].next = Rec.records[id].cpos;
	       myL[myL[i].next].prev = i;
	       incFreq (&Heap,id);
	     }
	  myL[i].prev = -id-1;
	  Rec.records[id].cpos = i;
if (PRNL && (i%10000 == 0)) printf ("Processed %i chars\n",i);
	}
    purgeHeap (&Heap);
  }

}

void prnSym(int c)
   { printf("%i",c); 
   }

void prnC (void)

  { int i = 0;
    printf ("C[1..%i] = ",c);
    while (i<u)
      { prnSym(C[i]);
  	printf (" ");
	i++; 
	if ((i<u) && (C[i] < 0)) i = -C[i]-1; 
      }
    printf ("\n\n");
// printf ("Physical C[1..%i] = ",u);
// for (i=0;i<u;i++)
//       { prnSym(C[i]);
//  	printf (" ");
//       }
//     printf ("\n\n");
  }

void prnRec (void)

  { int i;
    printf ("Active pairs:\n");
    for (i=0;i<Rec.size;i++)
        { printf ("\t(");
	  prnSym(Rec.records[i].pair.left);
          printf (",");
	  prnSym(Rec.records[i].pair.right);
	  printf ("), %i occs\n", Rec.records[i].freq);
	}
    printf ("\n");
  }

void BPE::createRule(ulong x, ulong y, int FREQ) {
  reemplazos->push(x,y);  //best solution in terms of memory.
  if(verbose_rules)
    printf("\ncreated rule : %i --> %li %li \n\n", 
        new_symbol ,x,y );
}

int BPE::ini_repair(){
  ini_repair(0);
}

//wrapper for Gonazlo Navarro's RePair
int BPE::ini_repair(ulong limFreq) 
{
  reemplazos = new mylist();
  //createRule(data[pos], data[successor(pos)],FREQ);
  { int oid,id,cpos;
    Trecord *rec,*orec;
    Tpair pair;
if (PRNC) prnC();
    while (new_symbol+1 > 0)
      { 
if (PRNR) prnRec();
	oid = extractMax(&Heap);
	if (oid == -1) break; // the end!!
 	orec = &Rec.records[oid];
  //artificial end!!!
  if(orec->freq<= limFreq) break;//the other end!:
 	cpos = orec->cpos;
  createRule(orec->pair.left,orec->pair.right,orec->freq);
if (PRNP) 
    { printf("Chosen pair %i = (",new_symbol);
      prnSym(orec->pair.left);
      printf(",");
      prnSym(orec->pair.right);
      printf(") (%i occs)\n",orec->freq);
    }
	while (cpos != -1)
	   { int ant,sgte,ssgte; 
		// replacing bc->e in abcd, b = cpos, c = sgte, d = ssgte
	     if (C[cpos+1] < 0) sgte = -C[cpos+1]-1; 
	     else sgte = cpos+1; 
	     if ((sgte+1 < u) && (C[sgte+1] < 0)) ssgte = -C[sgte+1]-1;
	     else ssgte = sgte+1; 
		// remove bc from myL
	     if (myL[cpos].next != -1) myL[myL[cpos].next].prev = -oid-1;
	     orec->cpos = myL[cpos].next;
	     if (ssgte != u) // there is ssgte
		{ 	// remove occ of cd
		  pair.left = C[sgte]; pair.right = C[ssgte];
		  id = searchHash(Hash,pair);
		  if (id != -1) // may not exist if purgeHeap'd
	             { if (id != oid) decFreq (&Heap,id); // not to my pair!
		       if (myL[sgte].prev != NullFreq) //still exists(not removed)
 		          { rec = &Rec.records[id];
		            if (myL[sgte].prev < 0) // this cd is head of its list
		               rec->cpos = myL[sgte].next;
		            else myL[myL[sgte].prev].next = myL[sgte].next;
		            if (myL[sgte].next != -1) // not tail of its list
		               myL[myL[sgte].next].prev = myL[sgte].prev;
			  }
		     }
			// create occ of ed
	          pair.left = new_symbol;
		  id = searchHash(Hash,pair);
	          if (id == -1) // new pair, insert
		     { id = insertRecord (&Rec,pair);
 		       rec = &Rec.records[id];
	               myL[cpos].next = -1;
		     }
	          else 
 		     { incFreq (&Heap,id);
 		       rec = &Rec.records[id]; 
	               myL[cpos].next = rec->cpos;
	               myL[myL[cpos].next].prev = cpos;
	             }
	          myL[cpos].prev = -id-1;
	          rec->cpos = cpos;
		}
	     if (cpos != 0) // there is ant
		{ 	// remove occ of ab
	          if (C[cpos-1] < 0) 
		     { ant = -C[cpos-1]-1; 
		       if (ant == cpos) // sgte and ant clashed -> 1 hole
			  ant = cpos-2;
		     }
	          else ant = cpos-1; 
		  pair.left = C[ant]; pair.right = C[cpos];
		  id = searchHash(Hash,pair);
		  if (id != -1) // may not exist if purgeHeap'd
	             { if (id != oid) decFreq (&Heap,id); // not to my pair!
		       if (myL[ant].prev != NullFreq) //still exists (not removed)
 		          { rec = &Rec.records[id];
		            if (myL[ant].prev < 0) // this ab is head of its list
		                 rec->cpos = myL[ant].next;
		            else myL[myL[ant].prev].next = myL[ant].next;
		            if (myL[ant].next != -1) // it is not tail of its list
		               myL[myL[ant].next].prev = myL[ant].prev;
			  }
		     }
			// create occ of ae
	          pair.right = new_symbol;
		  id = searchHash(Hash,pair);
	          if (id == -1) // new pair, insert
		     { id = insertRecord(&Rec,pair);
 		       rec = &Rec.records[id];
	               myL[ant].next = -1;
	             }
	          else 
	             { incFreq (&Heap,id);
 		       rec = &Rec.records[id];
	               myL[ant].next = rec->cpos;
	               myL[myL[ant].next].prev = ant;
	             }
	          myL[ant].prev = -id-1;
	          rec->cpos = ant;
		}
	     C[cpos] = new_symbol;
	     if (ssgte != u) C[ssgte-1] = -cpos-1;
	     C[cpos+1] = -ssgte-1;
	     c--;
 	     orec = &Rec.records[oid]; // just in case of Rec.records Realloc'd
	     cpos = orec->cpos;
	   }
if (PRNC) prnC();
	 removeRecord (&Rec,oid);
	 new_symbol++;
	 purgeHeap(&Heap); // remove freq 1 from heap
	 if (c < RPfactor * u) // compact C
	    { int i,ni;
	      i = 0;
	      for (ni=0;ni<c-1;ni++) 
		{ C[ni] = C[i];
		  myL[ni] = myL[i];
		  if (myL[ni].prev < 0) 
		     { if (myL[ni].prev != NullFreq) // real ptr
			  Rec.records[-myL[ni].prev-1].cpos = ni; 
		     }
		  else myL[myL[ni].prev].next = ni;
		  if (myL[ni].next != -1) myL[myL[ni].next].prev = ni;
		  i++; if (C[i] < 0) i = -C[i]-1;
		}
	      C[ni] = C[i];
	      u = c;
	      C = (int*)Realloc (C, c * sizeof(int));
	      myL = (Tlist*)Realloc (myL, c * sizeof(Tlist));
              assocRecords (&Rec,&Hash,&Heap,myL);
	    }
       }
     return 0;
   }
  return 0;
}


// Algorithm to compress the dictionary, by Rodrigo Gonzalez.
// Modified to store lengths and ones.
int BPE::compress_pair_table() 
{
  //DV interfacing
  max_value = 1;
  max_assigned = new_symbol-1;

  //symbols = data;


  ulong aux,j,i;
  aux=0;

  j=max_assigned;
  mynodo *a;
  symbols_pair = new ulong[2*(max_assigned-max_value)];
  if(!symbols_pair) return -1;
  while (!reemplazos->empty()) {
    a = reemplazos->pop();
    symbols_pair[2*(j-max_value-1)]=a->v1;
    symbols_pair[2*(j-max_value-1)+1]=a->v2;
    j--;
    free(a);
  }
  delete reemplazos;
  
  rule_length = new ulong[max_assigned-max_value];
  rule_ones = new ulong[max_assigned-max_value];
  
  if(!rule_length) return -1;
  if(!rule_ones) return -1;
  largo_max_real = 0;
  max_ones_real = 0;
  for(j=2;j<=max_assigned;j++){
    ulong x = symbols_pair[2*(j-max_value-1)];
    ulong y = symbols_pair[2*(j-max_value-1)+1];

    rule_length [j-2 ] = 0; 
    rule_length [j-2 ] += LENGTH(x); 
    rule_length [j-2 ] += LENGTH(y);
    rule_ones[j-2 ] = 0; 
    rule_ones[j-2 ] += ONES(x); 
    rule_ones[j-2 ] += ONES(y);
    if( rule_length [j-2 ] > largo_max_real)  
      largo_max_real = rule_length[j-2];  
    if( rule_ones [j-2 ] > max_ones_real)  
      max_ones_real = rule_length[j-2];  
  }
  //

  /////
  /* Compress table of pairs */
  //ulong *symbols_pair_tmp = (ulong*) malloc(sizeof(ulong)*(max_assigned-max_value)); // U array.
  ulong * symbols_pair_tmp = new ulong[max_assigned-max_value];
  if(!symbols_pair_tmp) return -1;
  for (i =0 ; i< (max_assigned-max_value) ; i++) 
    symbols_pair_tmp[i]=0;

  for (i =0 ; i< (max_assigned-max_value) ; i++) {
    aux=symbols_pair[2*i];
    if (aux > max_value) {  
      symbols_pair_tmp[aux-max_value-1]++;
    }
    aux=symbols_pair[2*i+1];
    if (aux > max_value) { //unin v
      symbols_pair_tmp[aux-max_value-1]++; 
    }
  }
  j=0; //number of refrenced rules.
  for (i =0 ; i< (max_assigned-max_value); i++) {
    if (symbols_pair_tmp[i] != 0) 
      j++;
  }
  symbols_new_len = 2*(max_assigned-max_value)-j;
  
  symbols_new = new ulong[symbols_new_len];
  symbols_new_bit = new ulong[symbols_new_len+(max_assigned-max_value)];
  ulong * symbols_new_value = new ulong[(max_assigned-max_value)];
  ulong * symbols_new_prof = new ulong[(max_assigned-max_value)];

  if(!symbols_new) return -1;
  if(!symbols_new_bit) return -1;
  if(!symbols_new_value) return -1;
  if(!symbols_new_prof) return -1;


  for (i =0 ; i<((symbols_new_len+(max_assigned-max_value))/W+1);i++)
    symbols_new_bit[i]=0;
  for (i =0 ; i<symbols_new_len;i++)
    symbols_new[i]=0;
  for (i =0 ; i<(max_assigned-max_value);i++)
    symbols_new_value[i]=0;
  for (i =0 ; i<(max_assigned-max_value);i++)
    symbols_new_prof[i]=0;
  j=1;
  ulong k1=0;
  int max_prof;
  int max_max=0;
  for (i =0 ; i< (max_assigned-max_value) ; i++) {
    if (symbols_pair_tmp[i] == 0) {
      // expand_rule:::::
      max_prof=1;
      symbols_new_value[i]=j; bitset(symbols_new_bit,j-1); j++;
      new_value(symbols_pair,symbols_new_value,&k1,&j,i,&max_prof, symbols_new_prof);
      symbols_new_prof[i] = max_prof;
      if(max_prof > max_max ) max_max = max_prof;
    }
  }
  this->max_prof = 1+max_max; 
  delete[](symbols_new_prof);
  symbols_new_bit_len = j;
  ulong *symbols_new_bit_aux = new ulong [(symbols_new_bit_len/W+1)];

  for (i =0 ; i<symbols_new_bit_len/W+1;i++) {
    symbols_new_bit_aux[i]= symbols_new_bit[i];
  }
  delete[](symbols_new_bit);
  symbols_new_bit=symbols_new_bit_aux;

  delete[](symbols_pair_tmp);
  ulong cuentame=0,cuentame2=0;

  BRR=new BitSequenceRG((uint*)symbols_new_bit, symbols_new_bit_len , 20);
  delete[]symbols_new_bit;

  ulong ones = BRR->rank1(symbols_new_bit_len);

  ulong nb1 = bits(BRR->rank1(symbols_new_bit_len)+max_value+1);
  this->nbits=nb1;

  ulong maped, unm;
  ulong * symbols_aux;
  symbols_aux = new ulong[((m*nbits/W)+1)];

  for(i=0;i<((m*nbits/W)+1);i++){
    symbols_aux[i]=0;
  }
  if(!symbols_aux) return -1;
  for (i =0 ; i< m ; i++) {
    aux = symbols[i];
    if (aux > max_value) {

      unm = symbols_new_value[aux-max_value-1]+max_value;
      maped = (BRR->rank1(unm-2))-1;
      SetField(symbols_aux,nbits,i,maped+2);
      cuentame2++;
    } else {
      SetField(symbols_aux,nbits,i,aux);
      cuentame++;
    }
  }
  free(symbols);
  symbols=symbols_aux;
  symbols_aux = new ulong[((symbols_new_len*nbits)/W+1)];
  if(!symbols_aux) return -1;
  for(i=0;i<((symbols_new_len*nbits)/W+1);i++){
    symbols_aux[i]=0;
  }
  for (i =0 ; i< symbols_new_len ; i++) {
    if(symbols_new[i] > max_value){ 
      unm = symbols_new[i];
      maped = (BRR->rank1(unm-2))-1;
      SetField(symbols_aux,nbits,i,maped+2);
    }
    else{
      SetField(symbols_aux,nbits,i,symbols_new[i]);
    } 
  }
  delete[](symbols_new);
  symbols_new=symbols_aux;


  delete[](symbols_pair); 

  //DV: largos y y unos de las reglas:
  s_lengths = new ulong[max_assigned-max_value];
  s_ones = new ulong[max_assigned-max_value];

  for (i =0 ; i<(max_assigned-max_value);i++){
    s_lengths[i]=-1;
    s_ones[i]=-1;
  }
  if(!s_lengths) return -1;
  if(!s_ones) return -1;

  for (i =0 ; i<(max_assigned-max_value);i++){
    s_lengths[BRR->rank1(symbols_new_value[i]-1)-1] = rule_length[i]; 
    s_ones   [BRR->rank1(symbols_new_value[i]-1)-1] = rule_ones[i]; 
  }


  delete[](rule_length);
  delete[](rule_ones);
  delete[](symbols_new_value);  
  return 0;
}

void BPE::compress_final_seq() 
{
}
void BPE::save_lengths() 
{
}

ulong My_bitget(ulong *e, ulong p){
  return   ((((e)[(p)/W] >> ((p)%W))) & 1);
}

//expand rule.
void  BPE::new_value(ulong *symbols_pair,ulong *symbols_new_value,ulong *k1,ulong *j,ulong pos, int *m_p, ulong * symbols_new_prof) {
  int mpi, mpd;
  mpi=0;mpd=0;
  ulong izq,der;
  izq=symbols_pair[2*pos];
  der=symbols_pair[2*pos+1];

  if (izq>max_value) {
    izq=izq-max_value-1;
    if (symbols_new_value[izq] == 0) {
      symbols_new_value[izq]=*j; bitset(symbols_new_bit,*j-1); (*j)++;
      new_value(symbols_pair,symbols_new_value,k1,j,izq,&mpi, symbols_new_prof);
      symbols_new_prof[izq] = mpi;
    } else {
      symbols_new[*k1]=symbols_new_value[izq]+max_value; (*j)++; (*k1)++;
      mpi = 1+symbols_new_prof[izq];
    }
  } else {
    symbols_new[*k1]=izq;(*j)++;(*k1)++;
    mpi = 1;
  }

  if (der>max_value) {
    der=der-max_value-1;
    if (symbols_new_value[der] == 0) {
      symbols_new_value[der]=*j; bitset(symbols_new_bit,*j-1); (*j)++;
      new_value(symbols_pair,symbols_new_value,k1,j,der,&mpd,symbols_new_prof);
      symbols_new_prof[der] = mpd;
    } else {
      symbols_new[*k1]=symbols_new_value[der]+max_value; (*j)++; (*k1)++;
      mpd = 1+symbols_new_prof[der];
    }
  } else {
    symbols_new[*k1]=der;(*j)++;(*k1)++;
    mpd = 1;
  }
  *m_p = *m_p+max(mpi,mpd);

}

BPE::BPE(){
}
int BPE::free_building_sructs(){  
  destroyHash(&Hash);
  destroyHeap(&Heap);
  destroyRecords(&Rec);
  free(myL);
  return 0;
}
BPE::BPE(ulong *a, ulong n, bool _verbose, int rules , int eT) {
  BPE(a, n, _verbose, rules , eT, 1);
}

//rules 
//      = -1  : every length and ones are stored
//      = -2  : just nbits  (or nbits/2 ? ) for each.
//      >=0    : max numbers of bits for each length  (and for each ones)
BPE::BPE(ulong *__a, ulong _n, bool _verbose, int rules , int eT, int limfreq) {
  this->T = eT;
  this->cod_rules = rules;
  this->gn_symb=NULL;
  this->gn_ones=NULL;
  this->gn_off=NULL;

  this->bits_r=-1;
  this->bits_o=-1;
  this->SOME_rules=false;//default!
  this->B_rules=NULL;
verbose=_verbose;
  /*
  //data=a;
  data = (ulong*) malloc((n+1)*sizeof(ulong));
  if(!data) return;
  for(ulong i=0;i<n;i++){
    data[i] = a[i];
    //printf("%lu",data[i]);
  }
  data[n] = NULL_DATA;
  this->length=n;
  //naux=n-1;
*  this->length=_n;
  if (verbose) { printf("Shift:\n"); fflush(stdout); }
*/
  this->length=_n;
  C = (int*)Malloc(length*sizeof(int));
  for(ulong i=0;i<length;i++)
    C[i] = __a[i];
  //free(a);
  delete[](__a);
  if (verbose) { printf("Shift:\n"); fflush(stdout); }
  shift_it();
  if (verbose) { printf("Prepare.\n"); fflush(stdout); }
  prepare();
   if (ini_repair(limfreq) != 0)
	{ fprintf (stderr,"Error in repair itself\n");
	  exit(1);
	}
   if(free_building_sructs() !=0)
   { fprintf (stderr,"Error freeing structures\n");
     exit(1);
   }

     //fprintf (stderr,"RePair succeeded\n\n");
     //fprintf (stderr,"   Original ints: %i\n",length);
     //fprintf (stderr,"   Number of rules: %i\n",n-alph);
     //fprintf (stderr,"   Final sequence length: %i\n",c);
     //fprintf (stderr,"   Compression ratio: %0.2f%%\n",
	   // (2.0*(n-alph)+c)*(float)bits(n-1)/(float)(length*bits(alph-1))*100.0);
     //exit(0);
  //
  //if (verbose) { printf("RePair.\n"); fflush(stdout); }
  //ini_repair();


  if (verbose) { printf("Compact Data\n"); fflush(stdout); }
  ini_compact();

  if (verbose) { printf("Compress pair table:\n"); fflush(stdout); }
  if(compress_pair_table()) return;
  compress_final_seq();//por ahora plano...


  build_GN_samples(T);

  this->largo_max = largo_max_real;

  this->ini_bes();//required for rank
  this->ones = this->rank1(length-1);
  //this->ones = length;//cuando rank no este bien implementad
  SOME_rules=false;//default!
  B_rules=NULL;

}


//get max space posible to invert in samples
ulong BPE::maxbits_rules(){
  return bits(largo_max_real)*(max_assigned-max_value);
}

void BPE::re_set_max_rule_len(ulong bits_max){
  MAXLEN_rules = true;
  
  SOME_rules=false;
  DAC_rules = false;

  if(bits(largo_max_real)<=bits_max){
    this->largo_max= largo_max_real; 
    return;
  } 
  this->largo_max = (1<<bits_max)-1;
  if((this->largo_max)<2) 
    this->largo_max=2;
}
void BPE::re_set_max_rule(ulong targ_bits){
  
  SOME_rules=false;
  
  ulong bits_max = targ_bits /(2*(max_assigned-max_value));
  if(bits_max<2){ 
    this->largo_max=2;
    return;
  }
  if(bits(largo_max_real)<=bits_max){
    this->largo_max= largo_max_real; 
    return;
  } 
  this->largo_max = (1<<bits_max)-1;

  ulong real =(((2*bits(largo_max))*(max_assigned-max_value))); 

  assert(real<=targ_bits);
}


//PARAM2 UNUSED
// Actual implementation of Max-Depth (MD)
// technique to store ones and lengths
void BPE::re_set_SOME_rules2(ulong param1, ulong param2){
  
  SOME_rules = true;
  
  DAC_rules = false;
  MAXLEN_rules = false;
  
  ulong nmaxdepth=0;
  //ulong HMAX = ((max_prof-2))/param1;
  ulong HMAX = param1;
  //printf("no expandir mas de %lu\n",HMAX);
  ulong level=0;
  ulong j,i,aux,k,aux2,aux1,rank_aux;
  
  largo_max = largo_max_real;
  
  ulong *b_rules_bit = new ulong [((max_assigned+1)/W+1)]; //+1?
  for(int i=0;i<((max_assigned+1)/W+1);i++)
    b_rules_bit[i]=0;
 
 /*
 Whole final sequence requieres too much space. 
 Just begining from param!
  */
  if(param1==0){
    for(int i=0;i<m;i++){
      aux = GetField(symbols,nbits,i); //lread from C 
      if(aux>=2)
        bitset(b_rules_bit,aux);//marco la regla que quiero guardar
    }
  }
  //

///////////////
//dispairall:
  ulong *info = (ulong*) malloc(sizeof(ulong)*this->length);
  ulong *levels = (ulong*) malloc(sizeof(ulong)*this->length);
  ulong rank_rule;
  i=length-1;//pos en T
  j=m-1;//pos en C

  ulong k_ch=0;
  while (j+1>0) {
    aux=GetField(symbols,nbits,j); //symbols[j];
    k=0;
    j--;
    if (aux <=max_value) {
      info[i]=aux+shift;
      i--; 
    } else {
      level=1;
      k=0;//usa el ppio de info como buffer para expandir niveles...
      rank_rule=aux-max_value;
      aux1=BRR->select1(rank_rule);
      assert(BRR->access(aux1));//root of the corrsp. tree.
      aux2=1; //contamos exceso de unos. cuando llega a cero terminamos el arbol.
      ulong nr=aux1-rank_rule+1; //posicion asociada en RS == symbols_new.
      rank_aux=nr;
      while (aux2!=0) {
        if (BRR->access(aux1))
          aux2++;
        else {
          aux2--;
          info[k]=GetField(symbols_new,nbits,rank_aux);//level 1
          levels[k]=1;
          rank_aux++;
          k++;
        }
        aux1++;
      }
    }
    while (k > 0) {
      aux=info[k-1];//POP
      level = levels[k-1];
      //if(!(level%HMAX))
      if((level==HMAX))
      {
        if(aux>=2)
          bitset(b_rules_bit,aux);//marco la regla que quiero guardar 
      } 
      //  printf("rule: %lu, level: %lu\n",aux,level);     
      if (aux <=max_value) {
        info[i]=aux+shift;
        i--;k--; 
      } else {
        level++;
        k_ch=k;
        if(level>nmaxdepth)
          nmaxdepth=level;
        rank_rule=aux-max_value-1;
        aux1=BRR->select1(rank_rule+1);

        assert(BRR->access(aux1));
        aux2=1;
        rank_aux=aux1-BRR->rank1(aux1)+1;
        while (aux2!=0) {
          if (BRR->access(aux1))
            aux2++;
          else {
            aux2--;
            //printf("level:%lu\n",level);
            info[k-1]=GetField(symbols_new,nbits,rank_aux);
            levels[k-1]=level;
            rank_aux++;
            k++; 
          }
          aux1++;
        }
        k--;
      }
    }
  }
  free(info);
  free(levels);
  //printf("nmd:%lu, omd:%lu\n",nmaxdepth,max_prof);
//////////////////////////
  if(B_rules!=NULL){
    delete[](v_lengths);
    delete[](v_ones);
    delete(B_rules);
  }
  //B_rules=new BBitRankW32Int(b_rules_bit,max_assigned, true, 20); 
  //B_rules=new BBitRankW32Int(b_rules_bit,max_assigned, 20); 
  B_rules=new BitSequenceRG((uint*)b_rules_bit,max_assigned+1, 20); 
  delete []b_rules_bit;
  
  n_rules = B_rules->rank1(max_assigned);//=='getOnes()'
  //v_lengths = (ulong*) malloc(n_rules*sizeof(ulong));
  //v_ones = (ulong*) malloc(n_rules*sizeof(ulong));
  v_lengths = new ulong[n_rules];
  v_ones = new ulong[n_rules];
  ulong my_longest_r=0;
  ulong my_longest_o=0;
  for(int i=0;i<n_rules;i++){
    aux = B_rules->select1(i+1);
    v_lengths[i] = s_lengths[aux-max_value-1];
    v_ones[i]    = s_ones[aux-max_value-1];
    if(v_lengths[i]>my_longest_r)
      my_longest_r = v_lengths[i];
    if(v_ones[i]>my_longest_o)
      my_longest_o = v_ones[i];
  }
  //printf("my l_r: %lu my l_o: %lu, largo maximo real: %lu\n",my_longest_r, my_longest_o,largo_max_real);
  if(verbose)
    printf("marcamos %lu reglas!\n",n_rules);
  //bits_v = bits(my_longest_r);
  bits_r = bits(my_longest_r);
  bits_o = bits(my_longest_o);
  //free(s_ones);//todo:descomentar.
  //free(s_lengths);
}


void BPE::re_set_DAC_for_rules(ulong param1){
  
  DAC_rules = true;
  
  SOME_rules = false;
  MAXLEN_rules = false;
  
  ulong nmaxdepth=0;
  //ulong HMAX = ((max_prof-2))/param1;
  ulong HMAX = param1;
  //printf("no expandir mas de %lu\n",HMAX);
  ulong level=0;
  ulong j,i,aux,k,aux2,aux1,rank_aux;
  
  largo_max = largo_max_real;
  
  ulong *b_rules_bit = new ulong [((max_assigned+1)/W+1)]; //+1?
  for(int i=0;i<((max_assigned+1)/W+1);i++)
    b_rules_bit[i]=0;
 
 /*
 Whole final sequence requieres too much space. 
 Just begining from param!
  */
  if(param1==0){
    for(int i=0;i<m;i++){
      aux = GetField(symbols,nbits,i); //lread from C 
      if(aux>=2)
        bitset(b_rules_bit,aux);//marco la regla que quiero guardar
    }
  }
  //

///////////////
//dispairall:
  ulong *info = (ulong*) malloc(sizeof(ulong)*this->length);
  ulong *levels = (ulong*) malloc(sizeof(ulong)*this->length);
  ulong rank_rule;
  i=length-1;//pos en T
  j=m-1;//pos en C

  ulong k_ch=0;
  while (j+1>0) {
    aux=GetField(symbols,nbits,j); //symbols[j];
    k=0;
    j--;
    if (aux <=max_value) {
      info[i]=aux+shift;
      i--; 
    } else {
      level=1;
      k=0;//usa el ppio de info como buffer para expandir niveles...
      rank_rule=aux-max_value;
      aux1=BRR->select1(rank_rule);
      assert(BRR->access(aux1));//root of the corrsp. tree.
      aux2=1; //contamos exceso de unos. cuando llega a cero terminamos el arbol.
      ulong nr=aux1-rank_rule+1; //posicion asociada en RS == symbols_new.
      rank_aux=nr;
      while (aux2!=0) {
        if (BRR->access(aux1))
          aux2++;
        else {
          aux2--;
          info[k]=GetField(symbols_new,nbits,rank_aux);//level 1
          levels[k]=1;
          rank_aux++;
          k++;
        }
        aux1++;
      }
    }
    while (k > 0) {
      aux=info[k-1];//POP
      level = levels[k-1];
      //if(!(level%HMAX))
      if((level==HMAX))
      {
        if(aux>=2)
          bitset(b_rules_bit,aux);//marco la regla que quiero guardar 
      } 
      //  printf("rule: %lu, level: %lu\n",aux,level);     
      if (aux <=max_value) {
        info[i]=aux+shift;
        i--;k--; 
      } else {
        level++;
        k_ch=k;
        if(level>nmaxdepth)
          nmaxdepth=level;
        rank_rule=aux-max_value-1;
        aux1=BRR->select1(rank_rule+1);

        assert(BRR->access(aux1));
        aux2=1;
        rank_aux=aux1-BRR->rank1(aux1)+1;
        while (aux2!=0) {
          if (BRR->access(aux1))
            aux2++;
          else {
            aux2--;
            //printf("level:%lu\n",level);
            info[k-1]=GetField(symbols_new,nbits,rank_aux);
            levels[k-1]=level;
            rank_aux++;
            k++; 
          }
          aux1++;
        }
        k--;
      }
    }
  }
  free(info);
  free(levels);
  //printf("nmd:%lu, omd:%lu\n",nmaxdepth,max_prof);
//////////////////////////
  if(B_rules!=NULL){
    delete(B_rules);
    delete[](v_lengths);
    delete[](v_ones);
  }
  //B_rules=new BBitRankW32Int(b_rules_bit,max_assigned, true, 20); 
  //B_rules=new BBitRankW32Int(b_rules_bit,max_assigned, 20); 
  B_rules=new BitSequenceRG((uint*)b_rules_bit,max_assigned+1, 20); 
  delete []b_rules_bit;
  
  n_rules = B_rules->rank1(max_assigned);//=='getOnes()'

  
  ////// new code 
  // moved to the class
  //factorization * dac_lengths;
  //factorization * dac_ones;
  
  uint * tmp_lengths = new uint[n_rules];
  uint * tmp_ones = new uint[n_rules];
  uint dac_sigma_ones, dac_sigma_lengths;
  dac_sigma_ones = 0;
  dac_sigma_lengths = 0;
  //uint ones_value, ones_offset;
  for(int i=0; i<n_rules; i++){
    aux = B_rules->select1(i+1);
    tmp_lengths[i] = s_lengths[aux-max_value-1];
    tmp_ones[i]    = s_ones[aux-max_value-1] + 1; // Plus 1 because DAC doesn't handle zeros
    //ones_value     = s_ones[aux-max_value-1];
    //ones_offset = tmp_lengths[i] - ones_value + 1; // +1 couse worst case ones==length, and dac do not handle zeros.
    //tmp_ones[i]    = ones_offset; 
    if(tmp_lengths[i] > dac_sigma_lengths) 
      dac_sigma_lengths = tmp_lengths[i];
    if(tmp_ones[i] > dac_sigma_ones) 
      dac_sigma_ones = tmp_ones[i];
    
    if(tmp_lengths[i]*tmp_ones[i] == 0 ) {
      cout << "warning, a zero in DAC in BPE" << endl;
      exit(-1);
    }
  }
  
  dac_lengths = new factorization(tmp_lengths,n_rules);
  dac_ones = new factorization(tmp_ones,n_rules);
  delete[] tmp_lengths;
  delete[] tmp_ones;
  // DAC HAS A BUG WHEN SIGMA IS SMALLER THAN 16. 
  // potencial improvement: allow different technique (SOME/DAC) for lengths/ones.
  if (dac_sigma_lengths < 16 || dac_sigma_ones < 16) {
    
    cout << endl;
    cout << "******************************************" << endl;
    cout << "******************************************" << endl;
    cout << "DAC will not work for rule sampling in BPE" << endl;
    cout << "We will resort to SOME technique." << endl;
    cout << "******************************************" << endl;
    cout << "******************************************" << endl;
    cout << endl;
    
    SOME_rules = true;
    DAC_rules = false;
    MAXLEN_rules = false;
  }
  
  if(SOME_rules) {
    v_lengths = new ulong[n_rules];
    v_ones = new ulong[n_rules];
    ulong my_longest_r=0;
    ulong my_longest_o=0;
    for(int i=0;i<n_rules;i++){
      aux = B_rules->select1(i+1);
      v_lengths[i] = s_lengths[aux-max_value-1];
      v_ones[i]    = s_ones[aux-max_value-1];
      if(v_lengths[i]>my_longest_r)
        my_longest_r = v_lengths[i];
      if(v_ones[i]>my_longest_o)
        my_longest_o = v_ones[i];
    }

    bits_r = bits(my_longest_r);
    bits_o = bits(my_longest_o);
    //free(s_ones);//todo:descomentar.
    //free(s_lengths);

    delete(dac_lengths);
    delete(dac_ones);
  } 
}


//assume contiguous = true
void BPE::build_GN_samples(ulong _factor){

  if(_factor!=0)
    this->factor=_factor;
  else
    this->factor=5;

  b=32;
  s=b*this->factor;

  if(gn_symb) delete[](gn_symb);
  if(gn_ones) delete[](gn_ones);
  if(gn_off)  delete[](gn_off);
  
  gn_symb= new ulong[1+length/s];
  gn_ones= new ulong[1+length/s];
  gn_off= new ulong[1+length/s];

  for(uint i=0;i<1+length/s;i++){
    gn_symb[i]=0; 
    gn_ones[i]=0; 
    gn_off[i]=0; 
  }

  ulong ones,pos;
  ones=0; pos=0;
  ulong aux,curr_l,curr_o,r_aux;
  ulong index=1;
  for(uint i=0; i<m; i++){
    aux = GetField(symbols,nbits,i); //data[i] 
    if(aux < 2){
      curr_l = 1;
      curr_o = aux;
    }else{
      r_aux = aux-2;
      curr_l = (s_lengths[r_aux]);
      curr_o = (s_ones[r_aux]);
    }

    while(pos+curr_l>= s*index){
      gn_symb[index] = i;
      gn_off [index] = (s*index)-pos; // belongs [1,..N]  
      gn_ones[index] = ones; //tmp, primero rank, dp access
      index++;
    } 
    pos  += curr_l;
    ones += curr_o;  
  }
  assert(pos==length);
}

ulong BPE::getNbits(){
  return this->nbits;
}

ulong BPE::DeltaSizeOfC(){
  ulong val, total,k,p,aux;
  total=0;
  for(uint i=0;i<m;i++){
    k=0;
    val=GetField(symbols,nbits,i);
    aux = val;
    while(aux) {
      aux >>= 1;
      k++;
    }
    aux = k;
    p = 0;
    while(aux) {
      aux >>= 1;
      p++;
    }

    total+= (2*p+k-2);
  }
  return total;
}

ulong BPE::huffmanSizeOfC(bool full){
  ulong maxC=0;
  for(uint i=0;i<m;i++) 
  {
    if(GetField(symbols,nbits,i)>maxC)
      maxC = GetField(symbols,nbits,i);
  }
  maxC++; 

  uint* freqs;
  freqs = (uint *)malloc(sizeof(int)*(maxC));
  for(uint i=0;i<(maxC);i++) freqs[i]=0;

  for(uint i=0;i<m;i++) 
  {
    ulong val = GetField(symbols,nbits,i);
    freqs[val]++;
  }

  THuff Hacc;
  uint i;

  int difs=0;
  for(i=0;i<maxC ; i++){
    if(freqs[i]>0) difs++;
  }


  Hacc = createHuff (freqs,(int)(maxC-1));
  ulong resp = Hacc.total;

  if(full){ 
    resp+=8*sizeHuff(Hacc);
    resp+=8*sizeof(Hacc);
  }
  free(freqs);
  freeHuff(Hacc);
  return resp;
}

void BPE::ini_compact(){
     
     int i = 0;
     m = 0 ;
     
     while (i<u)
        { 
          m++;
          i++; if ((i < u) && (C[i] < 0)) i = -C[i]-1;
        }
     assert(m==c);
     symbols = (ulong*) Malloc (m*sizeof(ulong));
     for(i=0;i<m;i++) symbols[i]=-6969;
     i=0;
     m=0;
     while (i<u)
        { 
          assert(C[i]>=0);
          symbols[m] = C[i];
          m++;
          i++; if ((i < u) && (C[i] < 0)) i = -C[i]-1;
        }
    free(C);
}


ulong *BPE::dispairall() {
  return C_dispairall();
}

ulong *BPE::C_dispairall() {
  ulong *info = (ulong*) malloc(sizeof(ulong)*this->length);
  ulong j,i,aux,k,aux2,aux1,rank_aux;
  ulong rank_rule;
  i=length-1;//pos en T
  j=m-1;//pos en C

  while (j+1>0) {
    aux=GetField(symbols,nbits,j); //symbols[j];
    k=0;
    j--;
    if (aux <=max_value) {
      info[i]=aux+shift;
      i--; 
    } else {

      k=0;//usa el ppio de info como buffer para expandir niveles...
      rank_rule=aux-max_value;
      aux1=BRR->select1(rank_rule);
      assert(BRR->access(aux1));//root of the corrsp. tree.
      aux2=1; //contamos exceso de unos. cuando llega a cero terminamos el arbol.
      ulong nr=aux1-rank_rule+1; //posicion asociada en RS == symbols_new.
      rank_aux=nr;
      while (aux2!=0) {
        if (BRR->access(aux1))
          aux2++;
        else {
          aux2--;
          info[k]=GetField(symbols_new,nbits,rank_aux);
          rank_aux++;
          k++;
        }
        aux1++;
      }
    }
    while (k > 0) {
      aux=info[k-1];
      if (aux <=max_value) {
        info[i]=aux+shift;
        i--;k--; 
      } else {
        rank_rule=aux-max_value-1;
        aux1=BRR->select1(rank_rule+1);

        assert(BRR->access(aux1));
        aux2=1;
        rank_aux=aux1-BRR->rank1(aux1)+1;
        while (aux2!=0) {
          if (BRR->access(aux1))
            aux2++;
          else {
            aux2--;
            info[k-1]=GetField(symbols_new,nbits,rank_aux);
            rank_aux++;
            k++; 
          }
          aux1++;
        }
        k--;
      }
    }
  }
  return info;
}

//adaptandome a libcds
void BPE::save(ofstream & of) const{
    //length=1;
  saveValue(of, REPAIR_HDR); //para estandarizar con libcds. Tipo del bitmap
  (saveValue (of,length));
  saveValue(of,m);
  saveValue(of,max_value);
  saveValue(of,max_assigned);
  saveValue(of,shift);
  saveValue(of,nbits);
  saveValue(of,symbols,((m/W)*nbits+((m%W)*nbits)/W+1));
  saveValue(of,symbols_new_len);
  saveValue(of,symbols_new,(symbols_new_len*nbits)/W+1);
  saveValue(of,s_lengths,max_assigned-max_value);
  saveValue(of,s_ones   ,max_assigned-max_value);
  saveValue(of,s);

  saveValue(of,gn_symb,1+length/s);
  saveValue(of,gn_ones,1+length/s);
  saveValue(of,gn_off,1+length/s);

  saveValue(of,max_prof);
  saveValue(of,largo_max_real);
  saveValue(of,max_ones_real);
  saveValue(of,largo_max);
  saveValue(of,bits_r);
  saveValue(of,bits_o);
  saveValue(of,cod_rules);

  saveValue<bool>(of,SOME_rules);
  saveValue<bool>(of,MAXLEN_rules);
  saveValue<bool>(of,DAC_rules);
  BRR->save(of);
  //B_rules->save(of);
}

//adaptandome a libcds
BPE * BPE::load(ifstream & f) {
  uint rd = loadValue<uint>(f);
  if (rd != REPAIR_HDR) return NULL;  //Para estandarizar con libcds
  BPE * ret = NULL;
  try
  {
  
    ret = new BPE();
  
  ret->length = loadValue<ulong>(f);
  ret->m = loadValue<ulong>(f);
  ret->max_value = loadValue<ulong>(f);
  ret->max_assigned = loadValue<ulong>(f);
  ret->shift = loadValue<ulong>(f);
  ret->nbits = loadValue<ulong>(f);
  //modelo: ret->O = loadValue<uint>(f,ret->O_len);
  
  ret->symbols = loadValue<ulong>(f,((ret->m/W)*ret->nbits+((ret->m%W)*ret->nbits)/W+1));

  //
  ret->symbols_new_len = loadValue<ulong>(f);
  
  ret->symbols_new = loadValue<ulong>(f,((ret->nbits)*(ret->symbols_new_len))/W+1);

  ret->s_lengths = loadValue<ulong>(f,ret->max_assigned-ret->max_value);

  ret->s_ones = loadValue<ulong>(f,ret->max_assigned-ret->max_value);

  //
  ret->s = loadValue<ulong>(f);

  ret->gn_symb = loadValue<ulong>(f,1+ret->length/ret->s);

  ret->gn_ones = loadValue<ulong>(f,1+ret->length/ret->s);

  ret->gn_off = loadValue<ulong>(f,1+ret->length/ret->s);


  ret->max_prof = loadValue<ulong>(f);
  ret->largo_max_real = loadValue<ulong>(f);
  ret->max_ones_real = loadValue<ulong>(f);
  ret->largo_max = loadValue<ulong>(f);
  ret->bits_r = loadValue<ulong>(f);
  ret->bits_o = loadValue<ulong>(f);
  ret->cod_rules = loadValue<int>(f);

  ret->SOME_rules = loadValue<bool>(f);
  ret->MAXLEN_rules = loadValue<bool>(f);
  ret->DAC_rules = loadValue<bool>(f);

  //ret->BRR = new BitSequenceRG(f);//TODO FIX LOAD WITH NEW LIBCDS
  ret->BRR = BitSequenceRG::load(f);
  //ret->B_rules = BitSequenceRG::load(f);

  //ret->ones = ret->length;
  ret->ini_bes();
  ret->ones = ret->rank1((ret->length)-1);


  //ret->SOME_rules = false; //tmp para comparar...
  verbose=false;
  return ret;
    
    //ret = new BitSequenceRRR();
    //uint type = loadValue<uint>(f);
    // TODO:throw an exception!
    //if(type!=RRR02_HDR) {
    //  abort();
    //}
   /*
   ret->length = loadValue<size_t>(f);
    ret->ones = loadValue<size_t>(f);
    ret->C_len = loadValue<uint>(f);
    ret->C_field_bits = loadValue<uint>(f);
    ret->O_len = loadValue<uint>(f);
    ret->O_bits_len = loadValue<uint>(f);
    ret->sample_rate = loadValue<uint>(f);
    ret->C = loadValue<uint>(f,uint_len(ret->C_len,ret->C_field_bits));
    ret->O = loadValue<uint>(f,ret->O_len);
    ret->create_sampling(ret->sample_rate);
   */
    return ret;
  }
  catch(exception e) {
    delete ret;
  }
  return NULL;
}


uint BPE::getSize() const{
  return this->size(); // call to the actual function
}
//soon deprecated
int BPE::save(FILE *f) {
  if (f == NULL) return 20;
  if (fwrite (&length,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&m,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&max_value,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&max_assigned,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&shift,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&nbits,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (symbols,sizeof(ulong),((m/W)*nbits+((m%W)*nbits)/W+1),f) != ((m/W)*nbits+((m%W)*nbits)/W+1)) return 21;
  if (fwrite (&symbols_new_len,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (symbols_new,sizeof(ulong),(symbols_new_len*nbits)/W+1,f) != (symbols_new_len*nbits)/W+1) return 21;
  if (fwrite (s_lengths,sizeof(ulong),max_assigned-max_value,f) != ((max_assigned-max_value))) return 21;
  if (fwrite (s_ones   ,sizeof(ulong),max_assigned-max_value,f) != ((max_assigned-max_value))) return 21;
  if (fwrite (&s,sizeof(ulong),1,f) != 1) return 21;

  if (fwrite (gn_symb,sizeof(ulong),1+length/s,f) != ((1+length/s))) return 21;
  if (fwrite (gn_ones,sizeof(ulong),1+length/s,f) != ((1+length/s))) return 21;
  if (fwrite (gn_off,sizeof(ulong),1+length/s,f) != ((1+length/s))) return 21;

  if (fwrite (&max_prof,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&largo_max_real,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&max_ones_real,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&largo_max,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&bits_r,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&bits_o,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&cod_rules,sizeof(int),1,f) != 1) return 21;

  //if (BRR->save(f) !=0) return 21; //TODO FIX save with new libcds
  printf("BPE successfully saved\n");
  return 0;
}

  int BPE::load(FILE *f) {
    if(f)
      return 0;
    else
      return 0;
  }

BPE::BPE(FILE *f, int *error) {
  *error = BPE::load(f);
}

BPE::~BPE() {
  free(b_1);
  free(b_2);

  delete[](symbols);
  delete[](symbols_new);

  delete[](gn_symb);
  delete[](gn_ones);
  delete[](gn_off);
  if(s_lengths)
    delete[](s_lengths);//
  if(s_ones)
    delete[](s_ones);//
  if(SOME_rules){
    delete(B_rules);
    delete[](v_lengths);
    delete[](v_ones);
  }
  delete BRR;
}

uint BPE::DictSize()
{
  ulong size=0; 
  size+=BRR->getSize();
  size+=((symbols_new_len*nbits+8-1)/8); 
  size+=sizeof(BPE);

  return size;

}

uint BPE::sizePlain() const
{
  ulong size=0; 
  size+=BRR->getSize();
  size+=((symbols_new_len*nbits+8-1)/8); 
  ulong C2=  (m*nbits+8-1)/8;

  size+=C2;

  size+=sizeof(BPE);


  return size;

}

uint BPE::getLargoMaxReal(){
  return largo_max_real;
}

//in bytes
uint BPE::sizeRulesLength() const{
  if(SOME_rules){
    return (((n_rules*(bits_r+bits_o))/8+1)+B_rules->getSize());
  }else if(MAXLEN_rules){
  if(largo_max>2)
    return (((2*bits(largo_max))*(max_assigned-max_value))+8-1)/8; 
  else 
    return 0;
  }else if(DAC_rules) {
    uint ans = 0;
    ans += dac_lengths->getSize(); 
    ans += dac_ones->getSize(); 
    ans += B_rules->getSize();
    return ans; 
  }
  else {
    cout << "NO METHOD FOR RULES STORAGE DEFINED! "<< endl;
    exit(-1);
  }
}

//in bytes
uint BPE::sizeSamples() const{
  return 3*(1+length/s)*sizeof(ulong);
}

// actual size function!!!!
uint BPE::size() const
{
  ulong size=0; 
  size+=sizePlain();
  size+=sizeRulesLength();
  size+=sizeSamples();
  return size;

}

double BPE::entropyOfC(){
  ulong maxC=0;
  for(uint i=0;i<m;i++) 
  {
    if(GetField(symbols,nbits,i)>maxC)
      maxC = GetField(symbols,nbits,i);
  }
  maxC++; 

  uint* freqs;
  freqs = (uint *)malloc(sizeof(int)*(maxC));
  for(uint i=0;i<(maxC);i++) freqs[i]=0;

  double total=0;
  for(uint i=0;i<m;i++) 
  {
    freqs[GetField(symbols,nbits,i)]++;
  }
  for(uint i=0;i<maxC;i++) 
  { 
    if(freqs[i]!=0)
      total += ((1.0*freqs[i]/m)* log(1.0*m/freqs[i])/log(2.0));
  }
  free(freqs);
  return total;
}

uint BPE::getM()
{
  return (this->m);
}

uint BPE::DeltaSize()
{
  ulong size=0;
  size+=BRR->getSize();//size+=BRR->SpaceRequirement();
  size+=((symbols_new_len*nbits+8-1)/8); 
  size+=(DeltaSizeOfC()+8-1)/8;
  size+=sizeof(BPE);

  return size;
}

uint BPE::huffSize()
{
  ulong size=0;
  size+=BRR->getSize();//size+=BRR->SpaceRequirement();
  size+=((symbols_new_len*nbits+8-1)/8); 
  size+=(huffmanSizeOfC(1)+8-1)/8;
  size+=sizeof(BPE);

  return size;
}

void BPE::cotas()

{
  double ent = entropyOfC();
  double e_huff    = 1.0*huffmanSizeOfC(0)/m;
  double full_huff = 1.0*huffmanSizeOfC(1)/m;
  
  printf("\n");

  printf("BPC in C:::\n");

  printf("plain    : %lu\n",nbits);
  printf("Huff     : %f\n" ,e_huff);
  printf("Huff Full: %f\n" ,full_huff);
  printf("Entro    : %f\n",ent); 
  printf("\n");


}

uint BPE::minSize()
{
  double entro;
  uint size=0;
  size+=BRR->getSize();//size+=BRR->SpaceRequirement();
  size+=((symbols_new_len*nbits+8-1)/8); 

  entro = entropyOfC();
  size+=(int)((m*entro+8-1)/8); 
  size+=sizeof(BPE);

  return size;
}


ulong BPE::getMaxProf(){
  return max_prof;
}

bool BPE::access(const size_t mi) const
{
  if(SOME_rules) return this->C_ARB_GNS_SOME_access(mi);
  else if(MAXLEN_rules) return this->C_ARB_GNS_access(mi);
  else if(DAC_rules) return this->C_DAC_access(mi);
}

size_t  BPE::rank1(const size_t mi) const
{
  if(mi+1==0) return 0; 
  if(mi > length) return ones;
  if(SOME_rules) return this->C_ARB_GNS_SOME_rank(mi);
  else if (MAXLEN_rules) return this->C_ARB_GNS_rank(mi);
  else if(DAC_rules) return this->C_DAC_rank(mi);
}

void BPE::ini_bes(){
  b_1= (int*) malloc(this->max_prof*sizeof(int));
  b_2= (int*) malloc(this->max_prof*sizeof(int));
  for(uint i=0;i<this->max_prof;i++){
    b_1[i]=0;
    b_2[i]=0;
  }
}

bool BPE::C_DAC_access(uint _mi) const{
  ulong aux,aux1,aux2,rank_aux;
  ulong curr_l;
  ulong pos;

  uint mi = _mi;
  uint mj = gn_symb[mi/s];
  pos = s*(mi/s)-gn_off[mi/s]; 

  int level=-1;
  bool baje = false;
  ulong rank_rule;
  bool conocida=true;
  ulong mpos;
  while(true){
    aux = GetField(symbols,nbits,mj); //lread from C 
    rank_rule=aux-max_value;
    conocida = B_rules->access(aux);
    if(aux <2){ 
      curr_l = 1;
      conocida=true;
    } 
    else{  
      conocida = B_rules->access(aux); 
      if(conocida){
        mpos = B_rules->rank1(aux)-1;
        //curr_l =  v_lengths[mpos];
        curr_l = dac_lengths->access(mpos+1);
      }
    }
    while(conocida && pos + (curr_l) <= mi )
    {
      pos+= curr_l; 
      mj++;
      aux = GetField(symbols,nbits,mj);
      rank_rule=aux-max_value;
      if(aux <2){ 
        curr_l = 1;
        conocida=true;
      } 
      else{  
        conocida = B_rules->access(aux); 
        if(conocida){
          mpos = B_rules->rank1(aux)-1;
          //curr_l =  v_lengths[mpos];
          curr_l = dac_lengths->access(mpos+1);
        }
      }
    }
    level++;
    while(level>=0){
      if(!baje){
        if (aux <=max_value && pos ==mi) return aux; 
        aux1=BRR->select1(rank_rule);
        assert(BRR->access(aux1));
        aux2=1; //number of 1s
        rank_aux=aux1-BRR->rank1(aux1)+1; 
      }
      else{
        //backtrack
        aux1 = b_1[level];
        aux2 = b_2[level];
        while(aux2!=0 && BRR->access(aux1) ){
          aux2++;
          aux1++;
        }
        aux1++;
        aux2--;
        baje = false; 
        rank_aux=aux1-BRR->rank1(aux1);  
        if(BRR->access(aux1)) rank_aux++; 
      }
      if(aux2!=0){
        aux = GetField(symbols_new,nbits,rank_aux);
        rank_rule=aux-max_value;
        if(aux <2){ 
          curr_l = 1;
          conocida=true;
        } 
        else{  
          conocida = B_rules->access(aux); 
          if(conocida){
            mpos = B_rules->rank1(aux)-1;
            //curr_l =  v_lengths[mpos];
            curr_l = dac_lengths->access(mpos+1);
          }
        }
      }
      while(aux2!=0 && conocida && pos+curr_l<=mi)
      {
        if (BRR->access(aux1))
          aux2++;
        else{
          aux2--;
          pos+= curr_l;
          rank_aux++;
          if(aux2==0) break;
          aux = GetField(symbols_new,nbits,rank_aux);
          rank_rule=aux-max_value;
          if(aux <2){ 
            curr_l = 1;
            conocida=true;
          } 
          else{  
            conocida = B_rules->access(aux); 
            if(conocida){
              mpos = B_rules->rank1(aux)-1;
              //curr_l =  v_lengths[mpos];
              curr_l = dac_lengths->access(mpos+1);
            }
          }
        }
        aux1++;
      }
      if(aux2==0){
        level-=1;
        baje = true; 
      }
      else{
        b_1[level] = aux1;
        b_2[level] = aux2;
        level ++;
      }
    }
    mj++;
    baje = false;
  }
  // This never happens.
  exit(-1);

  return 0;
}
bool BPE::C_ARB_GNS_access(uint _mi) const{
  ulong aux,aux1,aux2,rank_aux;
  ulong curr_l;
  ulong pos;

  uint mi = _mi;
  uint mj = gn_symb[mi/s];
  pos = s*(mi/s)-gn_off[mi/s]; 

  int level=-1;
  bool baje = false;
  ulong rank_rule;
  while(true){
    aux = GetField(symbols,nbits,mj); //lread from C 
    rank_rule=aux-max_value;
    if(aux <2) 
      curr_l = 1;
    else  
      curr_l =  s_lengths[aux-max_value-1];

    while(pos + (curr_l) <= mi && curr_l < largo_max ){
      pos+= curr_l; 
      mj++;
      aux = GetField(symbols,nbits,mj);
      rank_rule=aux-max_value;
      if(aux <2) 
        curr_l = 1;
      else
        curr_l =  s_lengths[aux-max_value-1];
    }
    level++;
    while(level>=0){
      if(!baje){
        if (aux <=max_value && pos ==mi) return aux; 
        aux1=BRR->select1(rank_rule);
        assert(BRR->access(aux1));
        aux2=1;//nro de unos. cuando llega a cero se debe expqandir solo uno mas.
        rank_aux=aux1-BRR->rank1(aux1); //attempt to solve a bug
        if(BRR->access(aux1)) rank_aux++; 
      }
      else{//backtrack
        aux1 = b_1[level];
        aux2 = b_2[level];
        while(aux2!=0 && BRR->access(aux1) ){
          aux2++;
          aux1++;
        }
        aux1++;
        aux2--;
        baje = false; 
        rank_aux=aux1-BRR->rank1(aux1); 
        if(BRR->access(aux1)) rank_aux++; 
      }
      if(aux2!=0){
        aux = GetField(symbols_new,nbits,rank_aux);
        rank_rule=aux-max_value;
        if(aux<2)
          curr_l=1;
        else{
          curr_l =  s_lengths[aux-max_value-1];
        }
      }
      while(aux2!=0 && pos+curr_l<=mi && curr_l < largo_max ){
        if (BRR->access(aux1))
          aux2++;
        else{
          aux2--;
          pos+= curr_l;
          rank_aux++;
          if(aux2==0) break;
          aux = GetField(symbols_new,nbits,rank_aux);
          rank_rule=aux-max_value;
          if(aux<2)
            curr_l=1;
          else{
            curr_l =  s_lengths[aux-max_value-1];
          }
        }
        aux1++;
      }
      if(aux2==0){
        level-=1;
        baje = true; 
      }
      else{
        b_1[level] = aux1;
        b_2[level] = aux2;
        level ++;
      }
    }
    //terminamos de expandir un symbolo de C y no retornamos.
    mj++;
    baje = false;
  }
  exit(-1);

  return 0;
}

bool BPE::C_ARB_GNS_SOME_access(uint _mi) const{
  ulong aux,aux1,aux2,rank_aux;
  ulong curr_l;
  ulong pos;

  uint mi = _mi;
  uint mj = gn_symb[mi/s];
  pos = s*(mi/s)-gn_off[mi/s]; 

  int level=-1;
  bool baje = false;
  ulong rank_rule;
  bool conocida=true;
  ulong mpos;
  while(true){
    aux = GetField(symbols,nbits,mj); //lread from C 
    rank_rule=aux-max_value;
    conocida = B_rules->access(aux);
    if(aux <2){ 
      curr_l = 1;
      conocida=true;
    } 
    else{  
      conocida = B_rules->access(aux); 
      if(conocida){
        mpos = B_rules->rank1(aux)-1;
        curr_l =  v_lengths[mpos];
      }
    }
    while(conocida && pos + (curr_l) <= mi )
    {
      pos+= curr_l; 
      mj++;
      aux = GetField(symbols,nbits,mj);
      rank_rule=aux-max_value;
      if(aux <2){ 
        curr_l = 1;
        conocida=true;
      } 
      else{  
        conocida = B_rules->access(aux); 
        if(conocida){
          mpos = B_rules->rank1(aux)-1;
          curr_l =  v_lengths[mpos];
        }
      }
    }
    level++;
    while(level>=0){
      if(!baje){
        if (aux <=max_value && pos ==mi) return aux; 
        aux1=BRR->select1(rank_rule);
        assert(BRR->access(aux1));
        aux2=1; //number of 1s
        rank_aux=aux1-BRR->rank1(aux1)+1; 
      }
      else{
        //backtrack
        aux1 = b_1[level];
        aux2 = b_2[level];
        while(aux2!=0 && BRR->access(aux1) ){
          aux2++;
          aux1++;
        }
        aux1++;
        aux2--;
        baje = false; 
        rank_aux=aux1-BRR->rank1(aux1);  
        if(BRR->access(aux1)) rank_aux++; 
      }
      if(aux2!=0){
        aux = GetField(symbols_new,nbits,rank_aux);
        rank_rule=aux-max_value;
        if(aux <2){ 
          curr_l = 1;
          conocida=true;
        } 
        else{  
          conocida = B_rules->access(aux); 
          if(conocida){
            mpos = B_rules->rank1(aux)-1;
            curr_l =  v_lengths[mpos];
          }
        }
      }
      while(aux2!=0 && conocida && pos+curr_l<=mi)
      {
        if (BRR->access(aux1))
          aux2++;
        else{
          aux2--;
          pos+= curr_l;
          rank_aux++;
          if(aux2==0) break;
          aux = GetField(symbols_new,nbits,rank_aux);
          rank_rule=aux-max_value;
          if(aux <2){ 
            curr_l = 1;
            conocida=true;
          } 
          else{  
            conocida = B_rules->access(aux); 
            if(conocida){
              mpos = B_rules->rank1(aux)-1;
              curr_l =  v_lengths[mpos];
            }
          }
        }
        aux1++;
      }
      if(aux2==0){
        level-=1;
        baje = true; 
      }
      else{
        b_1[level] = aux1;
        b_2[level] = aux2;
        level ++;
      }
    }
    mj++;
    baje = false;
  }
  // This never happens.
  exit(-1);

  return 0;
}

uint BPE::C_DAC_rank(uint _mi) const{
  ulong aux,aux1,aux2,rank_aux;
  ulong curr_l;
  ulong pos;

  ulong curr_o;
  ulong sum = 0;

  uint mi = _mi;
  uint mj = gn_symb[mi/s];
  pos = s*(mi/s)-gn_off[mi/s]; 
  sum = gn_ones[mi/s];

  int level=-1;
  bool baje = false;
  ulong rank_rule;
  bool conocida=true;
  ulong mpos;
  while(true){
    aux = GetField(symbols,nbits,mj); //lread from C 
    rank_rule=aux-max_value;
    conocida = B_rules->access(aux);//assert(!(aux<2) OR conocida=false)
    if(aux <2){ 
      curr_l = 1;
      conocida=true;
      curr_o = aux;
    } 
    else{  
      conocida = B_rules->access(aux); 
      if(conocida){
        mpos = B_rules->rank1(aux)-1;
        //curr_l =  v_lengths[mpos];
        //curr_o =  v_ones[mpos];
        curr_l = dac_lengths->access(mpos+1);
        curr_o = dac_ones->access(mpos+1)-1;
      }
    }
    while(conocida && pos + (curr_l) <= mi )
    {
      pos+= curr_l; 
      sum+=curr_o;
      mj++;
      aux = GetField(symbols,nbits,mj);
      rank_rule=aux-max_value;
      if(aux <2){ 
        curr_l = 1;
        conocida=true;
        curr_o = aux;
      } 
      else{  
        conocida = B_rules->access(aux); 
        if(conocida){
          mpos = B_rules->rank1(aux)-1;
          //curr_l =  v_lengths[mpos];
          //curr_o =  v_ones[mpos];
          curr_l = dac_lengths->access(mpos+1);
          curr_o = dac_ones->access(mpos+1)-1;
        }
      }
    }
    level++;
    while(level>=0){
      if(!baje){
        if (aux <=max_value && pos ==mi) return sum+aux; 
        aux1=BRR->select1(rank_rule);
        assert(BRR->access(aux1));
        aux2=1;
        rank_aux=aux1-BRR->rank1(aux1)+1; 
      }
      else{//backtrack
        aux1 = b_1[level];
        aux2 = b_2[level];
        while(aux2!=0 && BRR->access(aux1) ){
          aux2++;
          aux1++;
        }
        aux1++;
        aux2--;
        baje = false; 
        rank_aux=aux1-BRR->rank1(aux1); 
        if(BRR->access(aux1)) rank_aux++; 
      }
      if(aux2!=0){
        aux = GetField(symbols_new,nbits,rank_aux);
        rank_rule=aux-max_value;
        if(aux <2){ 
          curr_l = 1;
          conocida=true;
          curr_o = aux;
        } 
        else{  
          conocida = B_rules->access(aux); 
          if(conocida){
            mpos = B_rules->rank1(aux)-1;
            //curr_l =  v_lengths[mpos];
            //curr_o =  v_ones[mpos];
            curr_l = dac_lengths->access(mpos+1);
            curr_o = dac_ones->access(mpos+1)-1;
          }
        }
      }
      while(aux2!=0 && conocida && pos+curr_l<=mi)
      {
        if (BRR->access(aux1))
          aux2++;
        else{
          aux2--;
          pos+= curr_l;
          sum+=curr_o;
          rank_aux++;
          if(aux2==0) break;
          aux = GetField(symbols_new,nbits,rank_aux);
          rank_rule=aux-max_value;
          if(aux <2){ 
            curr_l = 1;
            conocida=true;
            curr_o = aux;
          } 
          else{  
            conocida = B_rules->access(aux); 
            if(conocida){
              mpos = B_rules->rank1(aux)-1;
              //curr_l =  v_lengths[mpos];
              //curr_o =  v_ones[mpos];
              curr_l = dac_lengths->access(mpos+1);
              curr_o = dac_ones->access(mpos+1)-1;
            }
          }
        }
        aux1++;
      }
      if(aux2==0){
        level-=1;
        baje = true; 
      }
      else{
        b_1[level] = aux1;
        b_2[level] = aux2;
        level ++;
      }
    }
    mj++;
    baje = false;
  }
  exit(-1);

  return 0;
}

uint BPE::C_ARB_GNS_SOME_rank(uint _mi) const{
  ulong aux,aux1,aux2,rank_aux;
  ulong curr_l;
  ulong pos;

  ulong curr_o;
  ulong sum = 0;

  uint mi = _mi;
  uint mj = gn_symb[mi/s];
  pos = s*(mi/s)-gn_off[mi/s]; 
  sum = gn_ones[mi/s];

  int level=-1;
  bool baje = false;
  ulong rank_rule;
  bool conocida=true;
  ulong mpos;
  while(true){
	
    aux = GetField(symbols,nbits,mj); //lread from C 
    rank_rule=aux-max_value;
    conocida = B_rules->access(aux);//assert(!(aux<2) OR conocida=false)
    if(aux <2){ 
      curr_l = 1;
      conocida=true;
      curr_o = aux;
    } 
    else{  
      conocida = B_rules->access(aux); 
      if(conocida){
        mpos = B_rules->rank1(aux)-1;
        curr_l =  v_lengths[mpos];
        curr_o =  v_ones[mpos];
      }
    }

    while(conocida && pos + (curr_l) <= mi )
    {
      pos+= curr_l; 
      sum+=curr_o;
      mj++;
      aux = GetField(symbols,nbits,mj);
      rank_rule=aux-max_value;
      if(aux <2){ 
        curr_l = 1;
        conocida=true;
        curr_o = aux;
      } 
      else{  
        conocida = B_rules->access(aux); 
        if(conocida){
          mpos = B_rules->rank1(aux)-1;
			fprintf(stderr,"hola");
          curr_l =  v_lengths[mpos];
			fprintf(stderr,"hola");
          curr_o =  v_ones[mpos];
        }
      }
    }
    level++;
    while(level>=0){
      if(!baje){
        if (aux <=max_value && pos ==mi) return sum+aux; 
        aux1=BRR->select1(rank_rule);
        assert(BRR->access(aux1));
        aux2=1;
        rank_aux=aux1-BRR->rank1(aux1)+1; 
      }
      else{//backtrack
        aux1 = b_1[level];
        aux2 = b_2[level];
        while(aux2!=0 && BRR->access(aux1) ){
          aux2++;
          aux1++;
        }
        aux1++;
        aux2--;
        baje = false; 
        rank_aux=aux1-BRR->rank1(aux1); 
        if(BRR->access(aux1)) rank_aux++; 
      }
      if(aux2!=0){
        aux = GetField(symbols_new,nbits,rank_aux);
        rank_rule=aux-max_value;
        if(aux <2){ 
          curr_l = 1;
          conocida=true;
          curr_o = aux;
        } 
        else{  
          conocida = B_rules->access(aux); 
          if(conocida){
            mpos = B_rules->rank1(aux)-1;
            curr_l =  v_lengths[mpos];
            curr_o =  v_ones[mpos];
          }
        }
      }
      while(aux2!=0 && conocida && pos+curr_l<=mi)
      {
        if (BRR->access(aux1))
          aux2++;
        else{
          aux2--;
          pos+= curr_l;
          sum+=curr_o;
          rank_aux++;
          if(aux2==0) break;
          aux = GetField(symbols_new,nbits,rank_aux);
          rank_rule=aux-max_value;
          if(aux <2){ 
            curr_l = 1;
            conocida=true;
            curr_o = aux;
          } 
          else{  
            conocida = B_rules->access(aux); 
            if(conocida){
              mpos = B_rules->rank1(aux)-1;
              curr_l =  v_lengths[mpos];
              curr_o =  v_ones[mpos];
            }
          }
        }
        aux1++;
      }
      if(aux2==0){
        level-=1;
        baje = true; 
      }
      else{
        b_1[level] = aux1;
        b_2[level] = aux2;
        level ++;
      }
    }
    mj++;
    baje = false;
  }
  exit(-1);

  return 0;
}

uint BPE::C_ARB_GNS_rank(uint _mi) const{
  ulong aux,aux1,aux2,rank_aux;
  ulong curr_l;
  ulong curr_o;
  ulong pos = 0;
  ulong sum = 0;
  int mj = 0;
  

  ulong mi = _mi;
  mj = gn_symb[mi/s];
  sum = gn_ones[mi/s];
  pos = s*(mi/s)-gn_off[mi/s]; 

  int level=-1;
  bool baje = false;
  ulong rank_rule;
  while(true){
    aux = GetField(symbols,nbits,mj); //lread from C 
    rank_rule=aux-max_value;
    if(aux <2){ 
      curr_l = 1;
      curr_o = aux;
    }
    else{  
      curr_l =  s_lengths[aux-max_value-1];
      curr_o =  s_ones[aux-max_value-1];
    }
    while(pos + (curr_l) <= mi && curr_l < largo_max ){
      pos+= curr_l; 
      sum+=curr_o;
      mj++;
      aux = GetField(symbols,nbits,mj);
      rank_rule=aux-max_value;
      if(aux <2){ 
        curr_l = 1;
        curr_o = aux;
      }
      else{
        curr_l =  s_lengths[aux-max_value-1];
        curr_o =  s_ones[aux-max_value-1];
      }
    }
    level++;
    while(level>=0){
      if(!baje){
        if (aux <=max_value && pos ==mi) return sum+aux; 
        aux1=BRR->select1(rank_rule);
        assert(BRR->access(aux1));
        aux2=1;
        rank_aux=aux1-rank_rule+1; 
      }
      else{//backtrack
        aux1 = b_1[level];
        aux2 = b_2[level];
        while(aux2!=0 && BRR->access(aux1) ){
          aux2++;
          aux1++;
        }
        aux1++;
        aux2--;
        baje = false; 
        rank_aux=aux1-BRR->rank1(aux1);  
        if(BRR->access(aux1)) rank_aux++; 
      }
      aux = GetField(symbols_new,nbits,rank_aux);
      rank_rule=aux-max_value;
      if(aux<2){
        curr_l=1;
        curr_o = aux;
      }
      else{
        curr_l =  s_lengths[aux-max_value-1];
        curr_o =  s_ones[aux-max_value-1];
      }
      while(aux2!=0 && pos+curr_l<=mi && curr_l < largo_max ){
        if (BRR->access(aux1))
          aux2++;
        else{
          aux2--;
          pos+= curr_l;
          sum+=curr_o;
          rank_aux++;
          aux = GetField(symbols_new,nbits,rank_aux);
          rank_rule=aux-max_value;
          if(aux<2){
            curr_l=1;
            curr_o = aux;
          }
          else{
            curr_l =  s_lengths[aux-max_value-1];
            curr_o =  s_ones[aux-max_value-1];
          }
        }
        aux1++;
      }
      if(aux2==0){
        level-=1;
        baje = true; 
      }
      else{
        b_1[level] = aux1;
        b_2[level] = aux2;
        level ++;
      }
    }
    mj++;
    baje = false;
  }
  exit(-1);
  return 0;

  return 0;
}


