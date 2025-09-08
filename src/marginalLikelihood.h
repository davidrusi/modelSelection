#ifndef MARGLHOOD_H
#define MARGLHOOD_H 1

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include "crossprodmat.h"
#include "Polynomial.h"


/*
 * Function Prototypes
 */
double marginalLikelihood(arma::SpMat<short> *model, lmObject *lm);

#endif /* MARGLHOOD_H */
