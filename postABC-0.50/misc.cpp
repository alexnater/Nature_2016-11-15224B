#include <vector>
#include <string>

using namespace std;



int splitStringtoVector(std::string s, std::vector<std::string> &vec, std::string delim){
	vec.clear();
	if( !s.empty() ){
		std::string::size_type l=s.find(delim);
		while(l!=std::string::npos){
			vec.push_back( s.substr(0,l) );
			s.erase(0, l+1);
			l=s.find_first_of(delim);
		}
		vec.push_back( s.substr(0,l) );
		}
	return 1;
	}

