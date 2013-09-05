#ifndef PAIRS_H_
#define PAIRS_H_

#include <iostream>
#include "constants.h"
#include "HashTable.h"

using namespace std;

class Pairs
{
public:
	Pairs();

  void reloadFirst(ulong * s_prev);	
	HASHREC ** get_hash();
	HASHREC * delete_pair(int w0, int w1);
	HASHREC * insert(int w0, int w1, Record * rec);
	HASHREC * search(int w0, int w1);

	virtual ~Pairs();
	
protected:
	HashTable pairs;

	int count;
};

#endif /*PAIRS_H_*/
