
#include "stringFunctions.h"

//-----------------------------------------------------------------------
// casting
//-----------------------------------------------------------------------
std::string toString(const int & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const long & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const std::vector<int>::size_type & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const unsigned int & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};


std::string toString(const float & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};

std::string toString(const double & input){
	std::ostringstream tos;
	tos << input;
	return tos.str();
};


int stringToInt(const std::string & s){
	return atoi(s.c_str());
};

long stringToLong(const std::string & s){
	return atol(s.c_str());
};

double stringToDouble(const std::string & s){
	return atof(s.c_str());
};

float stringToFloat(const std::string & s){
	return atof(s.c_str());
};
int stringToIntCheck(const std::string & s){
	return (int) stringToLongCheck(s);
};

long stringToLongCheck(const std::string & s){
	char** ptr=NULL;
	long i=strtol(s.c_str(), ptr, 10);
	if(ptr!=NULL) throw TException("String '" +s+"' is not a number!", _FATAL_ERROR);
	return i;
};

double stringToDoubleCheck(const std::string & s){
	char** ptr=NULL;
	double i=strtod(s.c_str(), ptr);
	if(ptr!=NULL) throw TException("String '" +s+"' is not a number!", _FATAL_ERROR);
	return i;
};

float stringToFloatCheck(const std::string & s){
	return (float) stringToDoubleCheck(s);
};


//-----------------------------------------------------------------------
//check
//-----------------------------------------------------------------------
bool stringContains(std::string & haystack, std::string needle){
	if(haystack.find(needle)==std::string::npos) return false;
	else return true;
};
bool stringContainsAny(std::string & haystack, std::string needle){
	if(haystack.find_first_of(needle)==std::string::npos) return false;
	else return true;
};
bool stringContains(std::string & haystack, char needle){
	if(haystack.find_first_of(needle)==std::string::npos) return false;
	else return true;
};
bool stringContainsOnly(std::string & haystack, std::string needle){
	if(haystack.find_first_not_of(needle)==std::string::npos) return true;
	else return false;
};

bool stringContainsLetters(std::string & haystack){
	return stringContainsAny(haystack, "abcdefghijklmnopqrstuvxyzABCDEFGHIJKLMNOPQRSTUVWXYZäöüÄÖÜàéèÀÉÈ");
};
bool stringContainsNumbers(std::string & haystack){
	return stringContainsAny(haystack, "1234567890");
};
bool allEntriesAreUnique(vector<std::string> vec){
	vector<std::string>::iterator it_second;
	for(vector<std::string>::iterator it=vec.begin(); it!= --vec.end(); ++it){
		it_second=it;
		++it_second;
		for(;it_second!=vec.end(); ++it){
			if((*it).compare(*it_second)==0) return false;
		}
	}
	return true;
}
std::string getFirstNonUniqueString(vector<std::string> vec){
	vector<std::string>::iterator it_second;
	for(vector<std::string>::iterator it=vec.begin(); it!= --vec.end(); ++it){
		it_second=it;
		++it_second;
		for(;it_second!=vec.end(); ++it_second){
			if((*it).compare(*it_second)==0) return *it;
		}
	}
	return "";
}

//-----------------------------------------------------------------------
//modify
//-----------------------------------------------------------------------
void eraseAllOccurences(std::string & s, std::string delim){
	std::string::size_type l=s.find(delim);
	while(l!=std::string::npos){
		s.erase(l, 1);
		l=s.find_first_of(delim);
	}
};

void eraseAllOccurencesAny(std::string & s, std::string delim){
	std::string::size_type l=s.find_first_of(delim);
	while(l!=std::string::npos){
		s.erase(l, 1);
		l=s.find_first_of(delim);
	}
};

void eraseAllWhiteSpaces(std::string & s){
	eraseAllOccurencesAny(s, " \t\f\v\n\r");
};

//-----------------------------------------------------------------------
//extract before
//-----------------------------------------------------------------------
std::string extractBefore(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBefore(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.find_first_of(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeAnyOf(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find_first_of(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeDoubleSlash(std::string & s){
	return extractBefore(s, "//");
};

std::string extractBeforeWhiteSpace(std::string & s){
	return extractBeforeAnyOf(s, " \t\f\v\n\r");
};

std::string extractUntil(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractUntil(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeLast(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractBeforeLast(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l);
		s.erase(0, l);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractUntilLast(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractUntilLast(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(0,l+1);
		s.erase(0, l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

//-----------------------------------------------------------------------
//extract after
//-----------------------------------------------------------------------
std::string extractAfter(std::string & s, std::string delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractAfter(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.find(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractAfterLast(std::string & s, std::string delim) {
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractAfterLast(std::string & s, char delim){
	std::string ret="";
	std::string::size_type l=s.rfind(delim);
	if(l!=std::string::npos){
		ret=s.substr(l+1);
		s.erase(l+1);
	} else {
		ret=s;
		s.clear();
	}
	return ret;
};

std::string extractPath(std::string & s){
	return extractUntilLast(s, "/");
};

void trimString(std::string & s){
	trimString(s, " \t\f\v\n\r");
};

void trimString(std::string & s, std::string what){
	//from beginning
	std::string::size_type l=s.find_first_not_of(what);
	if(l==std::string::npos){
		s.clear();
	} else {
		s.erase(0, l);
		//from end
		l=s.find_last_not_of(what);
		if(l!=std::string::npos)
			s.erase(l+1);
	}
};

void concatenateString(std::vector<std::string> & vec, std::string & s){
	s.clear();
	for(std::vector<std::string>::iterator it=vec.begin(); it!=vec.end(); ++it){
		s+=(*it);
	}
};
void concatenateString(std::vector<std::string> & vec, std::string & s, int from){
	s.clear();
	std::vector<std::string>::iterator it=vec.begin();
	it+=from;
	for(; it!=vec.end(); ++it){
		s+=(*it);
	}
}
void concatenateString(std::vector<std::string> & vec, std::string & s, std::string delim){
	s.clear();
	std::vector<std::string>::iterator it=vec.begin();
	s=*it;
	++it;
	for(; it!=vec.end(); ++it){
		s+=delim+(*it);
	}
};
//-----------------------------------------------------------------------
//split into vector
//-----------------------------------------------------------------------
void fillVectorFromString(std::string s, std::vector<std::string> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find(delim);
		while(l!=std::string::npos){
			vec.push_back(s.substr(0,l));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromStringAny(std::string s, std::vector<std::string> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(s.substr(0,l));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<std::string> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(tmp);
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<std::string> & vec){
	fillVectorFromStringAny(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<std::string> & vec){
	fillVectorFromStringAnySkipEmpty(s, vec, " \t\f\v\n\r");
}

void fillVectorFromString(std::string s, std::vector<std::string> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(s.substr(0,l));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(s.substr(0,l));
	}
};

void fillVectorFromString(std::string s, std::vector<float> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDouble(s.substr(0)));
	}
};

void fillVectorFromString(std::string s, std::vector<double> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDouble(s.substr(0)));
	}
};

void fillVectorFromStringAny(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDouble(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);		}
		vec.push_back(stringToDouble(s.substr(0)));
	}
};

void fillVectorFromStringAnyCheck(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToDoubleCheck(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(stringToDoubleCheck(s.substr(0)));
	}
};

void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(stringToDouble(tmp));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0);
		if(!tmp.empty()) vec.push_back(stringToDouble(s.substr(0)));
	}
};

void fillVectorFromStringAnySkipEmptyCheck(std::string s, std::vector<double> & vec, std::string delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		std::string tmp;
		while(l!=std::string::npos){
			tmp=s.substr(0,l);
			trimString(tmp);
			if(!tmp.empty()) vec.push_back(stringToDouble(tmp));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		tmp=s.substr(0);
		if(!tmp.empty()) vec.push_back(stringToDoubleCheck(s.substr(0)));
	}
};

void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAny(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceCheck(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAnyCheck(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAnySkipEmpty(s, vec, " \t\f\v\n\r");
};
void fillVectorFromStringWhiteSpaceSkipEmptyCheck(const std::string & s, std::vector<double> & vec){
	fillVectorFromStringAnySkipEmptyCheck(s, vec, " \t\f\v\n\r");
};

void fillVectorFromString(std::string s, std::vector<int> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(atoi(s.substr(0,l).c_str()));
	}
};

void fillVectorFromString(std::string s, std::vector<long> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToLong(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(atol(s.substr(0,l).c_str()));
	}
};


void fillVectorFromString(std::string s, std::vector<bool> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::string::size_type l=s.find_first_of(delim);
		while(l!=std::string::npos){
			vec.push_back(stringToInt(s.substr(0,l)));
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back(atoi(s.substr(0,l).c_str()));
	}
};

bool fillSequenceFromString(std::string s, std::vector<int> & vec, char delim){
	vec.clear();
	if(!s.empty()){
		std::vector<std::string> temp;
		fillVectorFromString(s, temp, delim);
		for(std::vector<std::string>::iterator it=temp.begin(); it!=temp.end(); ++it){
		   //If sequence fill sequence...
			std::string::size_type pos=it->find_first_of('-');
		   if(pos != std::string::npos){
			   int first=atoi((it->substr(0,pos).c_str()));
			   int second=atoi((it->substr(pos+1).c_str()));
			   if(second>first){
				   for(int j=first; j<=second; ++j) vec.push_back(j);
			   } else return false;
		   }
		   //if number, put back.
		   else {
			   int num=atoi(it->c_str());
			   vec.push_back(num);
		   }
		}
		return true;
	} else return false;
};

//-----------------------------------------------------------------------
//read from file
//-----------------------------------------------------------------------
void readHeaderAndValues(std::string & filename, std::vector<std::string> & header, std::vector<double> & values){
	std::string line;
	header.clear();
	values.clear();
	//open file stream
	std::ifstream is (filename.c_str()); // opening the file for reading
	if(!is) throw TException("The file '" + filename + "' could not be opened!", _FATAL_ERROR);

	//read header
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, header);

	//read observed Data
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, values);

	//close stream
	is.close();

	if(header.size() != values.size())
		throw TException("Number of values does not match number of header names in file '", filename, "'!", _FATAL_ERROR);
}
void readHeaderAndValuesUnique(std::string & filename, std::vector<std::string> & header, std::vector<double> & values){
	std::string line;
	header.clear();
	values.clear();
	//open file stream
	std::ifstream is (filename.c_str()); // opening the file for reading
	if(!is) throw TException("The file '" + filename + "' could not be opened!", _FATAL_ERROR);
	//read header
	getline(is, line);
	trimString(line);
	if(line.empty()) throw TException("The file '"+filename+"' appears to be empty!", _FATAL_ERROR);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, header);
	std::string s=getFirstNonUniqueString(header);
	if(!s.empty()) throw TException("Entry '"+s+"' is listed multiple times in header of file '"+ filename +"'!", _FATAL_ERROR);
	//read observed Data
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, values);
	//close stream
	is.close();

	if(header.size() != values.size())
		throw TException("Number of values does not match number of header names in file '", filename, "'!", _FATAL_ERROR);
}

void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromString(line, vec, delim);
};

void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, char delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromString(line, vec, delim);
};

void fillVectorFromLineAny(std::ifstream & is, std::vector<std::string> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringAny(line, vec, delim);
};

void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<std::string> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpace(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<std::string> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
};

void fillVectorFromLine(std::ifstream & is, std::vector<double> & vec, char delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromString(line, vec, delim);
};

void fillVectorFromLineAny(std::ifstream & is, std::vector<double> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringAny(line, vec, delim);
};
void fillVectorFromLineAnyCheck(std::ifstream & is, std::vector<double> & vec, std::string delim){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringAnyCheck(line, vec, delim);
}

void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<double> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpace(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<double> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmpty(line, vec);
};
void fillVectorFromLineWhiteSpaceSkipEmptyCheck(std::ifstream & is, std::vector<double> & vec){
	std::string line;
	getline(is, line);
	trimString(line);
	fillVectorFromStringWhiteSpaceSkipEmptyCheck(line, vec);
};

//-----------------------------------------------------------------------
//read from file
//-----------------------------------------------------------------------
std::string stringReplace(char needle, std::string replace, std::string & haystack){
	std::string s="";
	std::string::size_type l=haystack.find_first_of(needle);
	while(l!=std::string::npos){
		s=s+haystack.substr(0,l)+replace;
		haystack.erase(0, l+1);
		l=haystack.find_first_of(needle);
	}
	s=s+haystack;
	return s;
};








