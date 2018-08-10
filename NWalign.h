#ifndef NWALIGN_H
#define NWALIGN_H

#include "iostream"
#include "string"
#include "stdlib.h"
using namespace std;

int blosum62[24][24];
int* ali = NULL;	// save the aligned information // index start from 0
int* a_order = NULL;
int* b_order = NULL;

int a_len = 0;
int b_len = 0;

int final_sco = 0;

void loadBLOSUM62();
void run_needleman_wunsch(string, string, int, int);
int* mapAAinSeq2AAOrderInBLOSUM62(string seq);
int** new2DIntArr(int row, int col);
void release2DIntArr(int n, int ** Arr);
void print2DIntArr(int** arr, int row, int col);
void printAliInfo(const int*, string, string);

#endif