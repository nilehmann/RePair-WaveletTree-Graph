#include "Pairs.h"
#include "bpe.h"

Pairs::Pairs()
{
	pairs.inithashtable();
}


void
Pairs::reloadFirst(ulong * s_prev){
  pairs.reloadFirst(s_prev);
}

HASHREC **
Pairs::get_hash()
{
	return pairs.get();
}

HASHREC * 
Pairs::delete_pair(int w0, int w1)
{
  return pairs.delete_key(w0,w1);
}

HASHREC * 
Pairs::insert(int w0, int w1, Record * rec)
{
	return pairs.hashinsert(w0, w1, rec);
}

HASHREC *
Pairs::search(int w0, int w1)
{		
  return pairs.hashsearch(w0, w1);
}

Pairs::~Pairs()
{
	//printf("liberando pairs\n");
  //free(pairs);
}
