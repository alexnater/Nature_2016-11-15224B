//---------------------------------------------------------------------------

#include <iostream>
using namespace std;

#include "TException.h"
//---------------------------------------------------------------------------
/** Output */
ostream&
operator << (ostream& os, TException& E) {
   switch(E.getIntensity()){
      default:
      case _WARNING      : os << endl << "WARNING: ";      break;
      case _ERROR        : os << endl << "\nERROR: ";        break;
      case _FATAL_ERROR  : os << endl << "\nFATAL ERROR: ";  break;
   }
   os << E.getMessage() << endl;
   return os;
}


