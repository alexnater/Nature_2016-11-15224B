//---------------------------------------------------------------------------

#ifndef TExceptionH
#define TExceptionH


#include <sstream>

using namespace std;

//---------------------------------------------------------------------------
enum EXCEPTION {_WARNING, _ERROR, _FATAL_ERROR};

class TException{
//   friend istream& operator >> (istream& is, TException& E);
   friend ostream& operator << (ostream& os, TException& E);
   private:
      std::string message;
      EXCEPTION intensity; //0: warning; 1: error; 2: fatal error

   public:
      TException(){
         message="No error message!";
         intensity=_WARNING;
      }

      TException(std::string error){
    	  message = error;
    	  intensity=_WARNING;
      }

      TException(std::string error, int i){
    	  message = error;
    	  intensity=(EXCEPTION) i;
      }

      TException(std::string error, EXCEPTION i){
    	  message = error;
    	  intensity = i;
      }

      template<typename T>
      TException(std::string errorPartOne, T something, std::string errorPartTwo, EXCEPTION i){
    	  ostringstream tos;
    	  tos << something;
		  message = errorPartOne + tos.str() + errorPartTwo;
		  intensity = i;
      }
      ~TException(){}

      EXCEPTION getIntensity(){return intensity;}
      std::string getMessage(){return message;}
};

#endif
