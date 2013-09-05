#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <iostream>
#include <sstream>
#include "bpe.h"
#include <BitString.h>

//double ticks = 100.0;

using namespace std;

bool manual = false;

int main (int argc, char **argv)
{ 
	uint len = 10;
	
	BitString a(len);
  

	a.setBit(2,1);
	a.setBit(5,1);
	a.setBit(10,1);
  
	ulong *seq = new ulong[len];
    for(int i=0;i<len;i++)
        seq[i] = a.getBit(i);

    BPE *bpe = new BPE(seq,len, true,-1,len/4,1);
    bpe->build_GN_samples(128);
	int mp = bpe->getMaxProf();
	fprintf(stderr,"%d\n", mp);
    bpe->re_set_SOME_rules2(mp, mp);
	for(int i = 0; i < len; ++i)
	    printf("access(%d) = %d\n",i, bpe->access(i));


}
