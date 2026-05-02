#include "decl_modelsq.h"

//// Class IDENT05_MODELSQ
IDENT05_MODELSQ::IDENT05_MODELSQ(vectorop yvalues, int dsize, int dn) :
        size(dsize),
	n(dn),
	y(vectorop(dsize))
{
  for (int i=0; i<dsize; i++) {
     y.values[i]=yvalues(i);
  }
}

IDENT05_MODELSQ::IDENT05_MODELSQ(const IDENT05_MODELSQ &source) :
	size(source.size),
	n(source.n),
	y(vectorop(source.size))
{
  for (int i=0; i<size; i++) {
     y.values[i]=source.y.values[i];
  }
}

IDENT05_MODELSQ &IDENT05_MODELSQ::operator=(const IDENT05_MODELSQ &source)
{
  size=source.size;
  n=source.n;
  for (int i=0; i<size; i++) {
     y.values[i]=source.y.values[i];

  }
}

IDENT05_MODELSQ::~IDENT05_MODELSQ()
{
}


