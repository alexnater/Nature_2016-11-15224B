/*
 * stringFunctions.h
 *
 *  Created on: May 4, 2012
 *      Author: wegmannd
 */

#ifndef STRINGFUNCTIONS_H_
#define STRINGFUNCTIONS_H_

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <sstream>
#include "TException.h"
#include <stdio.h>

std::string toString(const int & input);
std::string toString(const long & input);
std::string toString(const std::vector<int>::size_type & input);
std::string toString(const unsigned int & input);
std::string toString(const float & input);
std::string toString(const double & input);

int stringToInt(const std::string & s);
long stringToLong(const std::string & s);
double stringToDouble(const std::string & s);
float stringToFloat(const std::string & s);

int stringToIntCheck(const std::string & s);
long stringToLongCheck(const std::string & s);
double stringToDoubleCheck(const std::string & s);
float stringToFloatCheck(const std::string & s);

//check
bool stringContains(std::string & haystack, std::string needle);
bool stringContainsAny(std::string & haystack, std::string needle);
bool stringContains(std::string & haystack, char needle);
bool stringContainsOnly(std::string & haystack, std::string needle);
bool stringContainsLetters(std::string & haystack);
bool stringContainsNumbers(std::string & haystack);
bool allEntriesAreUnique(vector<std::string> vec);
std::string getFirstNonUniqueString(vector<std::string> vec);

//modify
void eraseAllOccurences(std::string & s, std::string delim);
void eraseAllOccurencesAny(std::string & s, std::string delim);
void eraseAllWhiteSpaces(std::string & s);

//manipulations
std::string extractBefore(std::string & s, std::string delim);
std::string extractBefore(std::string & s, char delim);
std::string extractBeforeAnyOf(std::string & s, std::string delim);
std::string extractBeforeDoubleSlash(std::string & s);
std::string extractBeforeWhiteSpace(std::string & s);
std::string extractUntil(std::string & s, std::string delim);
std::string extractUntil(std::string & s, char delim);

std::string extractBeforeLast(std::string & s, std::string delim);
std::string extractBeforeLast(std::string & s, char delim);
std::string extractUntilLast(std::string & s, std::string delim);
std::string extractUntilLast(std::string & s, char delim);

std::string extractAfter(std::string & s, std::string delim);
std::string extractAfter(std::string & s, char delim);

std::string extractAfterLast(std::string & s, std::string delim);
std::string extractAfterLast(std::string & s, char delim);

std::string extractPath(std::string & s);

void trimString(std::string & s);
void trimString(std::string & s, std::string what);
void concatenateString(std::vector<std::string> & vec, std::string & s);
void concatenateString(std::vector<std::string> & vec, std::string & s, int from);
void concatenateString(std::vector<std::string> & vec, std::string & s, std::string delim);

//split into vector
void fillVectorFromString(std::string s, std::vector<std::string> & vec, std::string delim);
void fillVectorFromStringAny(std::string s, std::vector<std::string> & vec, std::string delim);
void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<std::string> & vec);
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<std::string> & vec);
void fillVectorFromString(std::string s, std::vector<std::string> & vec, char delim);
void fillVectorFromString(std::string s, std::vector<float> & vec, char delim);
void fillVectorFromString(std::string s, std::vector<double> & vec, char delim);
void fillVectorFromStringAny(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringAnyCheck(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringAnySkipEmpty(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringAnySkipEmptyCheck(std::string s, std::vector<double> & vec, std::string delim);
void fillVectorFromStringWhiteSpace(const std::string & s, std::vector<double> & vec);
void fillVectorFromStringWhiteSpaceCheck(const std::string & s, std::vector<double> & vec);
void fillVectorFromStringWhiteSpaceSkipEmpty(const std::string & s, std::vector<double> & vec);
void fillVectorFromStringWhiteSpaceSkipEmptyCheck(const std::string & s, std::vector<double> & vec);
void fillVectorFromString(std::string s, std::vector<int> & vec, char delim);
void fillVectorFromString(std::string s, std::vector<long> & vec, char delim);
void fillVectorFromString(std::string s, std::vector<bool> & vec, char delim);
bool fillSequenceFromString(std::string s, std::vector<int> & vec, char delim);

//read from file
void readHeaderAndValues(std::string & filename, std::vector<std::string> & header, std::vector<double> & values);
void readHeaderAndValuesUnique(std::string & filename, std::vector<std::string> & header, std::vector<double> & values);

void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, std::string delim);
void fillVectorFromLine(std::ifstream & is, std::vector<std::string> & vec, char delim);
void fillVectorFromLineAny(std::ifstream & is, std::vector<std::string> & vec, std::string delim);
void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<std::string> & vec);
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<std::string> & vec);

void fillVectorFromLine(std::ifstream & is, std::vector<double> & vec, char delim);
void fillVectorFromLineAny(std::ifstream & is, std::vector<double> & vec, std::string delim);
void fillVectorFromLineAnyCheck(std::ifstream & is, std::vector<double> & vec, std::string delim);
void fillVectorFromLineWhiteSpace(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmpty(std::ifstream & is, std::vector<double> & vec);
void fillVectorFromLineWhiteSpaceSkipEmptyCheck(std::ifstream & is, std::vector<double> & vec);
std::string stringReplace(char needle, std::string replace, std::string & haystack);

#endif /* STRINGFUNCTIONS_H_ */
