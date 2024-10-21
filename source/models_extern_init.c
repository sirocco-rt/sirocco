#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "atomic.h"
#include "sirocco.h"
#include "models.h"


int ncomps;                     // The number of components that have been read.  
int nmods_tot;                  // The total number of models that have been read in 

ModelPtr mods;

struct ModSum comp[NCOMPS];
