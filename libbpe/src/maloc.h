
#ifndef MEMSINCLUDED
#define MEMSINCLUDED
  // Memory management

#define malloc(n) Malloc(n)
#define free(p) Free(p)
#define realloc(p,n) Realloc(p,n)

void *Malloc (int n);
void Free (void *p);
void *Realloc (void *p, int n);

#endif
