#ifndef ForestFJRho_h
#define ForestFJRho_h

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestFJRho {
public :
   ForestFJRho(){};
   ~ForestFJRho(){};

   // Declaration of leaf types
   std::vector<double>  *etaMin = 0;
   std::vector<double>  *etaMax = 0;
   std::vector<double>  *rho    = 0;
   std::vector<double>  *rhom   = 0;

   // List of branches
   TBranch        *b_etaMin;   //!
   TBranch        *b_etaMax;   //!
   TBranch        *b_rho;      //!
   TBranch        *b_rhom;     //!
};
#endif
