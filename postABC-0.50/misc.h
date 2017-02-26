#ifndef MISC_H
#define MISC_H


#include <vector>
#include <string>


int splitStringtoVector(std::string s, std::vector<std::string> &vec, std::string delim);

template <typename T>
T max(int n, T *list){
	T max=0;
	for (int i=0; i<n; i++){
		if (list[i]>max) max=list[i];
		}
	return max;
	}

template <typename T>
T sum(int n, T *list){
	T total=0;
	for (int i=0; i<n; i++) total+=list[i];
	return total;
	}

template <typename T>
double average(int n, T *list){
	double avg=0;
	for (int i=0; i<n; i++) avg+=(double)list[i];
	avg /= n;
	return avg;
	}

template <typename T>
T * allocation(int size1){
	T *pointer=new T [size1]();
	return pointer;
	}

template <typename T>
T ** allocation(int size1, int size2){
	T **pointer=new T *[size1];
	for (int i=0; i<size1; i++){
		pointer[i]=new T [size2]();
		}
	return pointer;
	}

template <typename T>
T *** allocation(int size1, int size2, int size3){
	T ***pointer=new T **[size1];
	for (int i=0; i<size1; i++){
		pointer[i]=new T *[size2];
		for (int j=0; j<size2; j++){
			pointer[i][j]=new T [size3]();
			}
		}
	return pointer;
	}


template <typename T>
int freemem(T **pointer, int size1){
	for(int i=0; i<size1; i++) delete [] pointer[i];
	delete [] pointer;
	return 1;
	}

template <typename T>
int freemem(T ***pointer, int size1, int size2){
	for(int i=0; i<size1; i++){
		for(int j=0; j<size2; j++) delete [] pointer[i][j];
		delete [] pointer[i];
		}
	delete [] pointer;
	return 1;
	}

template <typename T>
int freemem_c(T **pointer, int size1){
	for(int i=0; i<size1; i++) free(pointer[i]);
	delete [] pointer;
	return 1;
	}

template <typename T>
int freemem_c(T ***pointer, int size1, int size2){
	for(int i=0; i<size1; i++){
		for(int j=0; j<size2; j++) free(pointer[i][j]);
		free(pointer[i]);
		}
	delete [] pointer;
	return 1;
	}



#endif
