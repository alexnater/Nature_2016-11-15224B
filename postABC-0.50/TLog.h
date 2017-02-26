/*
 * TLog.h
 *
 *  Created on: Oct 17, 2010
 *      Author: wegmannd
 */

#ifndef TLOG_H_
#define TLOG_H_

#include <iostream>
#include "TException.h"
#include <fstream>

using namespace std;

class TLog{
private:
	bool isFile;
	std::string filename;
	bool isVerbose;
	int numIndent;
	std::string indent;
	bool printWarnings;

public:
	ofstream file;

	TLog(){
		isFile=false;
		isVerbose=true;
		printWarnings=true;
		numIndent=0;
		fillIndentString();
	};

	void close(){
		if(isFile) file.close();
		isFile=false;
	};

	~TLog(){ close();};

	void openFile(std::string Filename){
		filename=Filename;
		file.open(filename.c_str());
		if(!file) throw TException("Unable to open logfile '"+ filename +"'!", _FATAL_ERROR);
		isFile=true;
	};

	void setVerbose(bool Verbose){ isVerbose=Verbose;};
	bool verbose(){return isVerbose;};
	void suppressWarings(){printWarnings=false;};
	void showWarings(){printWarnings=true;};

	void fillIndentString(){
		indent="";
		for(int i=0; i<numIndent; ++i) indent+="   ";
		indent+="- ";
	};

	void addIndent(int n=1){
		numIndent+=n;
		fillIndentString();
	};

	void removeIndent(int n=1){
		numIndent-=n;
		if(numIndent<0) numIndent=0;
		fillIndentString();
	};

	template<typename T>
	void startIndent(T out){
		list(out);
		addIndent();
	};

	template<typename T, typename U, typename V>
	void startIndent(T first, U middle, V last){
		list(first, middle, last);
		addIndent();
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void startIndent(T first, U second, V third, W fourth, X fifth){
		list(first, second, third, fourth, fifth);
		addIndent();
	};

	template<typename T>
	void startIndentFlush(T out){
		listFlush(out);
		addIndent();
	};

	template<typename T>
	void endIndent(T out){
		list(out);
		removeIndent();
	};

	template<typename T, typename U, typename V>
	void endIndent(T first, U middle, V last){
		list(first, middle, last);
		removeIndent();
	};

	void endIndent(){
		removeIndent();
	};

	void newLine(){
		if(isFile) file << endl;
		cout << endl;
	};

	template<typename T>
	void write(T out){
		if(isFile) file << out << endl;
		if(isVerbose) cout << out << endl;
	};

	template<typename T, typename V>
	void write(T first, V last){
		if(isFile) file << first << last << endl;
		if(isVerbose) cout << first << last << endl;
	};

	template<typename T, typename U, typename V>
	void write(T first, U middle, V last){
		if(isFile) file << first << middle << last << endl;
		if(isVerbose) cout << first << middle << last << endl;
	};

	template<typename T>
	void list(T out){
		if(isFile) file << indent << out << endl;
		if(isVerbose) cout << indent << out << endl;
	};

	template<typename T, typename U>
	void list(T first, U last){
		if(isFile) file << indent << first << last << endl;
		if(isVerbose) cout << indent << first << last << endl;
	};

	template<typename T, typename U, typename V>
	void list(T first, U middle, V last){
		if(isFile) file << indent << first << middle << last << endl;
		if(isVerbose) cout << indent << first << middle << last << endl;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void list(T first, U second, V third, W fourth, X fifth){
		if(isFile) file << indent << first << second << third << fourth << fifth << endl;
		if(isVerbose) cout << indent << first << second << third << fourth << fifth << endl;
	};


	template<typename T>
	void conclude(T out){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile) file << temp << "-> " << out << endl;
		if(isVerbose) cout << temp << "-> " << out << endl;
	};

	template<typename T, typename U>
	void conclude(T first, U last){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile) file << temp << "-> " << first << last << endl;
		if(isVerbose) cout << temp << "-> " << first << last << endl;
	};

	template<typename T, typename U, typename V>
	void conclude(T first, U middle, V last){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile) file << temp << "-> " << first << middle << last << endl;
		if(isVerbose) cout << temp << "-> " << first << middle << last << endl;
	};

	template<typename T, typename U, typename V, typename W>
	void conclude(T first, U middle, V secondlast, W last){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile) file << temp << "-> " << first << middle << secondlast <<  last << endl;
		if(isVerbose) cout << temp << "-> " << first << middle << secondlast << last << endl;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void conclude(T first, U second, V third, W fourth, X fifth){
		std::string temp="";
		for(int i=0; i<=numIndent; ++i) temp+="   ";
		if(isFile) file << temp << "-> " << first << second << third << fourth << fifth << endl;
		if(isVerbose) cout << temp << "-> " << first << second << third << fourth << fifth << endl;
	};

	template<typename T>
	void overWrite(T out){
		if(isFile) file << '\xd'<< out << endl;
		if(isVerbose) cout << '\xd' << out << endl;
	};

	template<typename T, typename U, typename V>
	void overWrite(T first, U middle, V last){
		if(isFile) file << '\xd'<< first << middle << last << endl;
		if(isVerbose) cout << '\xd' << first << middle << last << endl;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void overWrite(T first, U second, V third, W fourth, X fifth){
		if(isFile) file << '\xd'<< first << second << third << fourth << fifth << endl;
		if(isVerbose) cout << '\xd' << first << second << third << fourth << fifth << endl;
	};

	template<typename T>
	void overList(T out){
		if(isFile) file << '\xd'<< indent << out << endl;
		if(isVerbose) cout << '\xd' << indent << out << endl;
	};

	template<typename T, typename U, typename V>
	void overList(T first, U middle, V last){
		if(isFile) file << '\xd'<< indent << first << middle << last << endl;
		if(isVerbose) cout << '\xd' << indent << first << middle << last << endl;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void overList(T first, U second, V third, W fourth, X fifth){
		if(isFile) file << '\xd'<< indent << first << second << third << fourth << fifth << endl;
		if(isVerbose) cout << '\xd' << indent << first << second << third << fourth << fifth << endl;
	};

	template<typename T>
	void writeFileOnly(T out){
		if(isFile) file << out << endl;
	};

	template<typename T>
	void listFileOnly(T out){
		if(isFile) file << indent << out << endl;
	};

	template<typename T>
	void listNoFile(T out){
		cout << indent << out << endl;
	};

	template<typename T>
	void warning(T out){
		if(printWarnings){
			if(isFile) file << "WARNING: " << out << endl;
			//cout << "WARNING: " << out << endl;
			cerr << "WARNING: " << out << endl;
		}
	};

	template<typename T, typename U, typename V>
	void warning(T first, U middle, V last){
		if(printWarnings){
			if(isFile) file << "WARNING: " << first << middle << last << endl;
			//cout << "WARNING: " << out << endl;
			cerr << "WARNING: " << first << middle << last << endl;
		}
	};

	template<typename T>
	void error(T out){
		newLine();
		if(isFile) file << "ERROR: " << out << endl;
		//cout << "ERROR: " << out << endl;
		cerr << "ERROR: " << out << endl;
	};

	template<typename T, typename U, typename V>
	void error(T first, U middle, V last){
		newLine();
		if(isFile) file << "ERROR: " << first << middle << last << endl;
		//cout << "ERROR: " << out << endl;
		cerr << "ERROR: " << first << middle << last << endl;
	};

	template<typename T>
	void add(T out){
		if(isFile) file << out;
		if(isVerbose) cout << out;
	};

	template<typename T>
	void flush(T out){
			if(isFile) file << out << std::flush;
			if(isVerbose) cout << out << std::flush;
	};

	template<typename T, typename U, typename V>
	void flush(T first, U middle, V last){
		if(isFile) file << first << middle << last << std::flush;
		if(isVerbose) cout << first << middle << last << std::flush;
	};

	template<typename T>
	void listFlush(T out){
			if(isFile) file << indent << out << std::flush;
			if(isVerbose) cout << indent << out << std::flush;
	};

	template<typename T, typename U, typename V>
	void listFlush(T first, U middle, V last){
			if(isFile) file << indent << first << middle << last << std::flush;
			if(isVerbose) cout << indent << first << middle << last << std::flush;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void listFlush(T first, U second, V third, W fourth, X fifth){
			if(isFile) file << indent << first << second << third << fourth << fifth << std::flush;
			if(isVerbose) cout << indent << first << second << third << fourth << fifth << std::flush;
	};


	template<typename T>
	void flushFileOnly(T out){
		if(isFile) file << out << std::flush;
	};

	template<typename T>
	void listFlushFileOnly(T out){
		if(isFile) file << indent << out << std::flush;
	};


	template<typename T>
	void overFlush(T out){
			if(isFile) file << '\xd' << out << std::flush;
			if(isVerbose) cout << '\xd' << out << std::flush;
	};

	template<typename T, typename U, typename V>
	void overFlush(T first, U middle, V last){
		if(isFile) file << '\xd' << first << middle << last << std::flush;
		if(isVerbose) cout << '\xd' << first << middle << last << std::flush;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void overFlush(T first, U second, V third, W fourth, X fifth){
		if(isFile) file << '\xd' << first << second << third << fourth << fifth << std::flush;
		if(isVerbose) cout << '\xd' << first << second << third << fourth << fifth << std::flush;
	};

	template<typename T>
	void listOverFlush(T out){
			if(isFile) file << '\xd' << indent << out << std::flush;
			if(isVerbose) cout << '\xd' << indent << out << std::flush;
	};

	template<typename T, typename U, typename V>
	void listOverFlush(T first, U middle, V last){
			if(isFile) file << '\xd' << indent << first << middle << last << std::flush;
			if(isVerbose) cout << '\xd' << indent << first << middle << last << std::flush;
	};

	template<typename T, typename U, typename V, typename W, typename X>
	void listOverFlush(T first, U second, V third, W fourth, X fifth){
			if(isFile) file << '\xd' << indent << first << second << third << fourth << fifth << std::flush;
			if(isVerbose) cout << '\xd' << indent << first << second << third << fourth << fifth << std::flush;
	};


	std::string getFilename(){
		if(isFile) return filename;
		else return "";
	};
};



#endif /* TLOG_H_ */
