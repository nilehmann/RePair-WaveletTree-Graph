#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <iostream>
#include <sstream>
#include "bpe.h"

//double ticks = 100.0;

using namespace std;

bool manual = false;

int main (int argc, char **argv)
{ 
	int p[] = {0x035F, 0xFF00};
	ulong *p2 = new ulong[2];
	p2[0] = p[0];
	p2[1] = p[1];

    BPE *bpe = new BPE(p2,128 , true,-2,64);//DV

	fprintf(stderr,"%d\n", bpe->rank1(1));
  }
