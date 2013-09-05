// Memory management

// #include "basics.h" included later to avoid macro recursion for malloc
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>

int NullFreq = 1 << (8*sizeof(int)-1);
/* original one
*/
void *Malloc (int n)

   { void *p;
     if (n == 0) return NULL;
     p = (void*) malloc (n);
     if (p == NULL)
        { fprintf (stderr,"Could not allocate %i bytes\n",n);
          exit(1);
        }
     fflush(stdout);
     return p;
   }

//invento para valgrind:
/*
void *Malloc (int n)

   { 
     void *r;
     void *s;
     
     r = (void*) malloc (1000*n);
     free(r);

     void *p;
     if (n == 0) return NULL;
     p = (void*) malloc (n);
     if (p == NULL)
        { fprintf (stderr,"Could not allocate %i bytes\n",n);
          exit(1);
        }
     fflush(stdout);

     
     s = (void*) malloc (1000*n);
     free(s);
     
     return p;
   }
*/
void Free (void *p)

   { 
     fflush(stdout);
     if (p) free (p);
   }

void *Realloc (void *p, int n)

   { if (p == NULL) return Malloc (n);
     if (n == 0) { Free(p); return NULL; }
     p = (void*) realloc (p,n);
     if (p == NULL)
        { fprintf (stderr,"Could not allocate %i bytes\n",n);
          exit(1);
        }
     return p;
}
#include "maloc.h"
