#include "NWalign.h"

int main(int argc, char* args[]){
	
	if (argc < 3){
		cout << endl;
		cout << "query : Query sequence" << endl;
		cout << "templ : Template sequence" << endl;
		cout << endl << "e.g.,:" << endl;
		cout << "./NWalign GKVFLTNAFSINMLKEFPTTITIDKLDEEDFCLKLELRLEDGTLINAIGHDSTINLVNTL VQGAGVVGETPTIPPNTAYQYTSGTVLDTPFGIMYGTYGMVSESGEHFNAIIKPFRLATP" << endl;
		cout << endl;
		exit(0);
	}

	string a_seq(args[1]);
	string b_seq(args[2]);
	
	a_len = a_seq.length();
	b_len = b_seq.length();

	a_order = mapAAinSeq2AAOrderInBLOSUM62(a_seq);
	b_order = mapAAinSeq2AAOrderInBLOSUM62(b_seq);

	loadBLOSUM62();	
	run_needleman_wunsch(a_seq, b_seq, -11, -1);	

	printAliInfo(ali, a_seq, b_seq);

	delete[] ali; ali = NULL;
	delete[] b_order; b_order = NULL;
	delete[] a_order; a_order = NULL;
	return 1;
}

void run_needleman_wunsch(string a_seq, string b_seq, int gap_open, int gap_extend){
	int i = 0, j = 0, k = 0;

	int** S = new2DIntArr(a_len+1, b_len+1); // saving the score matrix
	int** D = new2DIntArr(a_len+1, b_len+1); // saving the direction 0: start, 1 : diag, 2 : horizontal, 3, vertical 
	int** H = new2DIntArr(a_len+1, b_len+1); // saving the horizontal matrix
	int** V = new2DIntArr(a_len+1, b_len+1); // saving the vertical matrix
	int** H_jump_steps = new2DIntArr(a_len+1, b_len+1);
	int** V_jump_steps = new2DIntArr(a_len+1, b_len+1);

	S[1][0] = gap_open;
	for (i = 2; i < a_len+1; i++){
		S[i][0] = S[i-1][0] + gap_extend;
	}
	for (i = 1; i < a_len+1; i++){
		H[i][0] = S[i][0];
	}

	S[0][1] = gap_open;
	for (j = 2; j < b_len+1; j++){
		S[0][j] = S[0][j-1] + gap_extend;
	}
	for (j = 1; j < b_len+1; j++){
		V[0][j] = S[0][j];
	}
	
	//run dynamical programming
	for (i = 1; i <= a_len; i++){
		for (j = 1; j <= b_len; j++){
			// Diagonal direction 
			int d_val = S[i-1][j-1] + blosum62[a_order[i-1]][b_order[j-1]];
			
			// Horizontal direction
			H_jump_steps[i][j] = 1;
			int h_val = S[i][j-1] + gap_open;
			int _h_val = H[i][j-1] + gap_extend;
			if (h_val < _h_val){
				h_val = _h_val;
				H_jump_steps[i][j] = H_jump_steps[i][j-1] + 1;
			}

			// vertical direction
			V_jump_steps[i][j] = 1;
			int v_val = S[i-1][j] + gap_open;
			int _v_val = V[i-1][j] + gap_extend;
			if (v_val < _v_val){
				v_val = _v_val;
				V_jump_steps[i][j] = V_jump_steps[i-1][j] + 1;
			}

			H[i][j] = h_val;
			V[i][j] = v_val;
		
			/*
			if (d_val >= h_val && d_val >= v_val){
				D[i][j] = 1;
				S[i][j] = d_val;
			}else if (h_val >= v_val){
				D[i][j] = 2;
				S[i][j] = h_val;
			}else{
				D[i][j] = 3;
				S[i][j] = v_val;
			} 
			*/

			if (d_val > h_val && d_val > v_val){
				D[i][j] = 1;
				S[i][j] = d_val;
			}else if (h_val > v_val){
				D[i][j] = 2;
				S[i][j] = h_val;
			}else{
				D[i][j] = 3;
				S[i][j] = v_val;
			}
		}//END FOR J
	}// END FOR I

	// tracking back
	ali = new int[a_len]; // ali[0] is useless.
	for (i = 0; i < a_len; i++)
		ali[i] = -1;

	i = a_len;
	j = b_len;
	while (i > 0 && j > 0){
		if (1 == D[i][j]){
			ali[i-1] = j-1; // index start from 0
			i--;
			j--;
		}else if (2 == D[i][j]){
			j = j-H_jump_steps[i][j];
			j = j>0 ? j : 0; 
		}else if (3 == D[i][j]){
			i = i-V_jump_steps[i][j];
			i = i>0 ? i : 0; 
		}
	}
	
	final_sco = S[a_len][b_len];
	
	release2DIntArr(a_len+1, V_jump_steps);
	release2DIntArr(a_len+1, H_jump_steps);
	release2DIntArr(a_len+1, V);
	release2DIntArr(a_len+1, H);
	release2DIntArr(a_len+1, D);
	release2DIntArr(a_len+1, S);
}

//The return array stores the amino acid order numbers of query seq.
//Please release the return value when it is useless.
int* mapAAinSeq2AAOrderInBLOSUM62(string seq){
	string AA_order = "*ARNDCQEGHILKMFPSTWYVBZX";	       // uppercase amino acide order in the BLAST's scoring matrix (e.g.,Blosum62). 
	string aa_order = "*arndcqegjilkmfpstwyvbzx";	       // lowercase amino acide order in the BLAST's scoring matrix (e.g.,Blosum62). 
	
	unsigned int len = seq.length();
	
	int* ans = new int[len];
	for (unsigned int i = 0; i < len; i++){
		int ind = AA_order.find(seq[i]);
		if (ind == -1)
			ind = aa_order.find(seq[i]);

		ans[i] = ind;
	}

	return ans;
}

void printAliInfo(const int* ali, string a_seq, string b_seq){
	int identity_num=0;
	int aligned_num=0;
	
	int i = 0;
	for(i = 0; i < a_len; i++){
		if(ali[i] != -1){
			int j = ali[i];
	    	aligned_num++;
		    
			if(a_order[i] == b_order[j])
				identity_num++;
		}
	}
	cout << endl;
	cout << "First (query) sequence length     : " << a_len << endl;
	cout << "Second (template) sequence length : " << b_len << endl;
	cout << "Number of identical pairs : " << identity_num << endl;
	cout << "Number of aligned pairs : " << aligned_num << endl;
	cout << "Sequence identity : " << (1.0*identity_num / a_len) << " (" << identity_num << "/" << a_len << ")" << endl;
	cout << endl;

	char* aligned_a_seq = new char[a_len+b_len+1];
	char* aligned_tags = new char[a_len+b_len+1];	// aligned tag
	char* aligned_b_seq = new char[a_len+b_len+1];

	int k = 0;
	int pre_ali_val = -1;
	for (i = 0; i < a_len; i++){
		if (i > 0 && -1 != ali[i-1])
			pre_ali_val = ali[i-1];

		if (-1 == ali[i]){
			aligned_a_seq[k] = a_seq[i];
			aligned_tags[k] = ' ';
			aligned_b_seq[k] = '-';

			k++;
		}else{
			int bias = ali[i] - pre_ali_val;
			for (int j = 0; j < bias-1; j++){
				aligned_a_seq[k] = '-';
				aligned_tags[k] = ' ';
				aligned_b_seq[k] = b_seq[ pre_ali_val+j+1 ];

				k++;
			}
			
			if (a_order[i] == b_order[ ali[i] ]){
				aligned_a_seq[k] = a_seq[i];
				aligned_tags[k] = ':';
				aligned_b_seq[k] = b_seq[ ali[i] ];
			}else{
				aligned_a_seq[k] = a_seq[i];
				aligned_tags[k] = ' ';
				aligned_b_seq[k] = b_seq[ ali[i] ];
			}
			k++;
		}
	}
	
	if (i > 0 && -1 != ali[a_len-1])
			pre_ali_val = ali[a_len-1];

	for (i = pre_ali_val+1; i < b_len; i++){
			aligned_a_seq[k] = '-';
			aligned_tags[k] = ' ';
			aligned_b_seq[k] = b_seq[i];

			k++;
	}

	aligned_a_seq[k] = '\0';
	aligned_tags[k] = '\0';
	aligned_b_seq[k] = '\0';


	cout << aligned_a_seq << endl;
	cout << aligned_tags << endl;
	cout << aligned_b_seq << endl;
	cout << endl;

	delete[] aligned_a_seq; aligned_a_seq = NULL;
	delete[] aligned_tags; aligned_tags = NULL;
	delete[] aligned_b_seq; aligned_b_seq = NULL;
}

void print2DIntArr(int** arr, int row, int col){
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			cout << arr[i][j] << "\t";
		}
		cout << endl;
	}
}

int** new2DIntArr(int row, int col){
	int **ans=new int*[row];
	for(int i=0;i<row;i++){
		ans[i]=new int[col];
		for(int j=0; j<col; j++)
			ans[i][j] = 0;
	}
	
	return ans;
}

void release2DIntArr(int n, int ** Arr){
	for(int i = 0; i < n; i++){
		delete [] Arr[i];
	}
	delete [] Arr;
	Arr = NULL;
}

void loadBLOSUM62(){
		blosum62[1][1]=4;
		blosum62[1][2]=-1;
		blosum62[1][3]=-2;
		blosum62[1][4]=-2;
		blosum62[1][5]=0;
		blosum62[1][6]=-1;
		blosum62[1][7]=-1;
		blosum62[1][8]=0;
		blosum62[1][9]=-2;
		blosum62[1][10]=-1;
		blosum62[1][11]=-1;
		blosum62[1][12]=-1;
		blosum62[1][13]=-1;
		blosum62[1][14]=-2;
		blosum62[1][15]=-1;
		blosum62[1][16]=1;
		blosum62[1][17]=0;
		blosum62[1][18]=-3;
		blosum62[1][19]=-2;
		blosum62[1][20]=0;
		blosum62[1][21]=-2;
		blosum62[1][22]=-1;
		blosum62[1][23]=0;
		blosum62[2][1]=-1;
		blosum62[2][2]=5;
		blosum62[2][3]=0;
		blosum62[2][4]=-2;
		blosum62[2][5]=-3;
		blosum62[2][6]=1;
		blosum62[2][7]=0;
		blosum62[2][8]=-2;
		blosum62[2][9]=0;
		blosum62[2][10]=-3;
		blosum62[2][11]=-2;
		blosum62[2][12]=2;
		blosum62[2][13]=-1;
		blosum62[2][14]=-3;
		blosum62[2][15]=-2;
		blosum62[2][16]=-1;
		blosum62[2][17]=-1;
		blosum62[2][18]=-3;
		blosum62[2][19]=-2;
		blosum62[2][20]=-3;
		blosum62[2][21]=-1;
		blosum62[2][22]=0;
		blosum62[2][23]=-1;
		blosum62[3][1]=-2;
		blosum62[3][2]=0;
		blosum62[3][3]=6;
		blosum62[3][4]=1;
		blosum62[3][5]=-3;
		blosum62[3][6]=0;
		blosum62[3][7]=0;
		blosum62[3][8]=0;
		blosum62[3][9]=1;
		blosum62[3][10]=-3;
		blosum62[3][11]=-3;
		blosum62[3][12]=0;
		blosum62[3][13]=-2;
		blosum62[3][14]=-3;
		blosum62[3][15]=-2;
		blosum62[3][16]=1;
		blosum62[3][17]=0;
		blosum62[3][18]=-4;
		blosum62[3][19]=-2;
		blosum62[3][20]=-3;
		blosum62[3][21]=3;
		blosum62[3][22]=0;
		blosum62[3][23]=-1;
		blosum62[4][1]=-2;
		blosum62[4][2]=-2;
		blosum62[4][3]=1;
		blosum62[4][4]=6;
		blosum62[4][5]=-3;
		blosum62[4][6]=0;
		blosum62[4][7]=2;
		blosum62[4][8]=-1;
		blosum62[4][9]=-1;
		blosum62[4][10]=-3;
		blosum62[4][11]=-4;
		blosum62[4][12]=-1;
		blosum62[4][13]=-3;
		blosum62[4][14]=-3;
		blosum62[4][15]=-1;
		blosum62[4][16]=0;
		blosum62[4][17]=-1;
		blosum62[4][18]=-4;
		blosum62[4][19]=-3;
		blosum62[4][20]=-3;
		blosum62[4][21]=4;
		blosum62[4][22]=1;
		blosum62[4][23]=-1;
		blosum62[5][1]=0;
		blosum62[5][2]=-3;
		blosum62[5][3]=-3;
		blosum62[5][4]=-3;
		blosum62[5][5]=9;
		blosum62[5][6]=-3;
		blosum62[5][7]=-4;
		blosum62[5][8]=-3;
		blosum62[5][9]=-3;
		blosum62[5][10]=-1;
		blosum62[5][11]=-1;
		blosum62[5][12]=-3;
		blosum62[5][13]=-1;
		blosum62[5][14]=-2;
		blosum62[5][15]=-3;
		blosum62[5][16]=-1;
		blosum62[5][17]=-1;
		blosum62[5][18]=-2;
		blosum62[5][19]=-2;
		blosum62[5][20]=-1;
		blosum62[5][21]=-3;
		blosum62[5][22]=-3;
		blosum62[5][23]=-2;
		blosum62[6][1]=-1;
		blosum62[6][2]=1;
		blosum62[6][3]=0;
		blosum62[6][4]=0;
		blosum62[6][5]=-3;
		blosum62[6][6]=5;
		blosum62[6][7]=2;
		blosum62[6][8]=-2;
		blosum62[6][9]=0;
		blosum62[6][10]=-3;
		blosum62[6][11]=-2;
		blosum62[6][12]=1;
		blosum62[6][13]=0;
		blosum62[6][14]=-3;
		blosum62[6][15]=-1;
		blosum62[6][16]=0;
		blosum62[6][17]=-1;
		blosum62[6][18]=-2;
		blosum62[6][19]=-1;
		blosum62[6][20]=-2;
		blosum62[6][21]=0;
		blosum62[6][22]=3;
		blosum62[6][23]=-1;
		blosum62[7][1]=-1;
		blosum62[7][2]=0;
		blosum62[7][3]=0;
		blosum62[7][4]=2;
		blosum62[7][5]=-4;
		blosum62[7][6]=2;
		blosum62[7][7]=5;
		blosum62[7][8]=-2;
		blosum62[7][9]=0;
		blosum62[7][10]=-3;
		blosum62[7][11]=-3;
		blosum62[7][12]=1;
		blosum62[7][13]=-2;
		blosum62[7][14]=-3;
		blosum62[7][15]=-1;
		blosum62[7][16]=0;
		blosum62[7][17]=-1;
		blosum62[7][18]=-3;
		blosum62[7][19]=-2;
		blosum62[7][20]=-2;
		blosum62[7][21]=1;
		blosum62[7][22]=4;
		blosum62[7][23]=-1;
		blosum62[8][1]=0;
		blosum62[8][2]=-2;
		blosum62[8][3]=0;
		blosum62[8][4]=-1;
		blosum62[8][5]=-3;
		blosum62[8][6]=-2;
		blosum62[8][7]=-2;
		blosum62[8][8]=6;
		blosum62[8][9]=-2;
		blosum62[8][10]=-4;
		blosum62[8][11]=-4;
		blosum62[8][12]=-2;
		blosum62[8][13]=-3;
		blosum62[8][14]=-3;
		blosum62[8][15]=-2;
		blosum62[8][16]=0;
		blosum62[8][17]=-2;
		blosum62[8][18]=-2;
		blosum62[8][19]=-3;
		blosum62[8][20]=-3;
		blosum62[8][21]=-1;
		blosum62[8][22]=-2;
		blosum62[8][23]=-1;
		blosum62[9][1]=-2;
		blosum62[9][2]=0;
		blosum62[9][3]=1;
		blosum62[9][4]=-1;
		blosum62[9][5]=-3;
		blosum62[9][6]=0;
		blosum62[9][7]=0;
		blosum62[9][8]=-2;
		blosum62[9][9]=8;
		blosum62[9][10]=-3;
		blosum62[9][11]=-3;
		blosum62[9][12]=-1;
		blosum62[9][13]=-2;
		blosum62[9][14]=-1;
		blosum62[9][15]=-2;
		blosum62[9][16]=-1;
		blosum62[9][17]=-2;
		blosum62[9][18]=-2;
		blosum62[9][19]=2;
		blosum62[9][20]=-3;
		blosum62[9][21]=0;
		blosum62[9][22]=0;
		blosum62[9][23]=-1;
		blosum62[10][1]=-1;
		blosum62[10][2]=-3;
		blosum62[10][3]=-3;
		blosum62[10][4]=-3;
		blosum62[10][5]=-1;
		blosum62[10][6]=-3;
		blosum62[10][7]=-3;
		blosum62[10][8]=-4;
		blosum62[10][9]=-3;
		blosum62[10][10]=4;
		blosum62[10][11]=2;
		blosum62[10][12]=-3;
		blosum62[10][13]=1;
		blosum62[10][14]=0;
		blosum62[10][15]=-3;
		blosum62[10][16]=-2;
		blosum62[10][17]=-1;
		blosum62[10][18]=-3;
		blosum62[10][19]=-1;
		blosum62[10][20]=3;
		blosum62[10][21]=-3;
		blosum62[10][22]=-3;
		blosum62[10][23]=-1;
		blosum62[11][1]=-1;
		blosum62[11][2]=-2;
		blosum62[11][3]=-3;
		blosum62[11][4]=-4;
		blosum62[11][5]=-1;
		blosum62[11][6]=-2;
		blosum62[11][7]=-3;
		blosum62[11][8]=-4;
		blosum62[11][9]=-3;
		blosum62[11][10]=2;
		blosum62[11][11]=4;
		blosum62[11][12]=-2;
		blosum62[11][13]=2;
		blosum62[11][14]=0;
		blosum62[11][15]=-3;
		blosum62[11][16]=-2;
		blosum62[11][17]=-1;
		blosum62[11][18]=-2;
		blosum62[11][19]=-1;
		blosum62[11][20]=1;
		blosum62[11][21]=-4;
		blosum62[11][22]=-3;
		blosum62[11][23]=-1;
		blosum62[12][1]=-1;
		blosum62[12][2]=2;
		blosum62[12][3]=0;
		blosum62[12][4]=-1;
		blosum62[12][5]=-3;
		blosum62[12][6]=1;
		blosum62[12][7]=1;
		blosum62[12][8]=-2;
		blosum62[12][9]=-1;
		blosum62[12][10]=-3;
		blosum62[12][11]=-2;
		blosum62[12][12]=5;
		blosum62[12][13]=-1;
		blosum62[12][14]=-3;
		blosum62[12][15]=-1;
		blosum62[12][16]=0;
		blosum62[12][17]=-1;
		blosum62[12][18]=-3;
		blosum62[12][19]=-2;
		blosum62[12][20]=-2;
		blosum62[12][21]=0;
		blosum62[12][22]=1;
		blosum62[12][23]=-1;
		blosum62[13][1]=-1;
		blosum62[13][2]=-1;
		blosum62[13][3]=-2;
		blosum62[13][4]=-3;
		blosum62[13][5]=-1;
		blosum62[13][6]=0;
		blosum62[13][7]=-2;
		blosum62[13][8]=-3;
		blosum62[13][9]=-2;
		blosum62[13][10]=1;
		blosum62[13][11]=2;
		blosum62[13][12]=-1;
		blosum62[13][13]=5;
		blosum62[13][14]=0;
		blosum62[13][15]=-2;
		blosum62[13][16]=-1;
		blosum62[13][17]=-1;
		blosum62[13][18]=-1;
		blosum62[13][19]=-1;
		blosum62[13][20]=1;
		blosum62[13][21]=-3;
		blosum62[13][22]=-1;
		blosum62[13][23]=-1;
		blosum62[14][1]=-2;
		blosum62[14][2]=-3;
		blosum62[14][3]=-3;
		blosum62[14][4]=-3;
		blosum62[14][5]=-2;
		blosum62[14][6]=-3;
		blosum62[14][7]=-3;
		blosum62[14][8]=-3;
		blosum62[14][9]=-1;
		blosum62[14][10]=0;
		blosum62[14][11]=0;
		blosum62[14][12]=-3;
		blosum62[14][13]=0;
		blosum62[14][14]=6;
		blosum62[14][15]=-4;
		blosum62[14][16]=-2;
		blosum62[14][17]=-2;
		blosum62[14][18]=1;
		blosum62[14][19]=3;
		blosum62[14][20]=-1;
		blosum62[14][21]=-3;
		blosum62[14][22]=-3;
		blosum62[14][23]=-1;
		blosum62[15][1]=-1;
		blosum62[15][2]=-2;
		blosum62[15][3]=-2;
		blosum62[15][4]=-1;
		blosum62[15][5]=-3;
		blosum62[15][6]=-1;
		blosum62[15][7]=-1;
		blosum62[15][8]=-2;
		blosum62[15][9]=-2;
		blosum62[15][10]=-3;
		blosum62[15][11]=-3;
		blosum62[15][12]=-1;
		blosum62[15][13]=-2;
		blosum62[15][14]=-4;
		blosum62[15][15]=7;
		blosum62[15][16]=-1;
		blosum62[15][17]=-1;
		blosum62[15][18]=-4;
		blosum62[15][19]=-3;
		blosum62[15][20]=-2;
		blosum62[15][21]=-2;
		blosum62[15][22]=-1;
		blosum62[15][23]=-2;
		blosum62[16][1]=1;
		blosum62[16][2]=-1;
		blosum62[16][3]=1;
		blosum62[16][4]=0;
		blosum62[16][5]=-1;
		blosum62[16][6]=0;
		blosum62[16][7]=0;
		blosum62[16][8]=0;
		blosum62[16][9]=-1;
		blosum62[16][10]=-2;
		blosum62[16][11]=-2;
		blosum62[16][12]=0;
		blosum62[16][13]=-1;
		blosum62[16][14]=-2;
		blosum62[16][15]=-1;
		blosum62[16][16]=4;
		blosum62[16][17]=1;
		blosum62[16][18]=-3;
		blosum62[16][19]=-2;
		blosum62[16][20]=-2;
		blosum62[16][21]=0;
		blosum62[16][22]=0;
		blosum62[16][23]=0;
		blosum62[17][1]=0;
		blosum62[17][2]=-1;
		blosum62[17][3]=0;
		blosum62[17][4]=-1;
		blosum62[17][5]=-1;
		blosum62[17][6]=-1;
		blosum62[17][7]=-1;
		blosum62[17][8]=-2;
		blosum62[17][9]=-2;
		blosum62[17][10]=-1;
		blosum62[17][11]=-1;
		blosum62[17][12]=-1;
		blosum62[17][13]=-1;
		blosum62[17][14]=-2;
		blosum62[17][15]=-1;
		blosum62[17][16]=1;
		blosum62[17][17]=5;
		blosum62[17][18]=-2;
		blosum62[17][19]=-2;
		blosum62[17][20]=0;
		blosum62[17][21]=-1;
		blosum62[17][22]=-1;
		blosum62[17][23]=0;
		blosum62[18][1]=-3;
		blosum62[18][2]=-3;
		blosum62[18][3]=-4;
		blosum62[18][4]=-4;
		blosum62[18][5]=-2;
		blosum62[18][6]=-2;
		blosum62[18][7]=-3;
		blosum62[18][8]=-2;
		blosum62[18][9]=-2;
		blosum62[18][10]=-3;
		blosum62[18][11]=-2;
		blosum62[18][12]=-3;
		blosum62[18][13]=-1;
		blosum62[18][14]=1;
		blosum62[18][15]=-4;
		blosum62[18][16]=-3;
		blosum62[18][17]=-2;
		blosum62[18][18]=11;
		blosum62[18][19]=2;
		blosum62[18][20]=-3;
		blosum62[18][21]=-4;
		blosum62[18][22]=-3;
		blosum62[18][23]=-2;
		blosum62[19][1]=-2;
		blosum62[19][2]=-2;
		blosum62[19][3]=-2;
		blosum62[19][4]=-3;
		blosum62[19][5]=-2;
		blosum62[19][6]=-1;
		blosum62[19][7]=-2;
		blosum62[19][8]=-3;
		blosum62[19][9]=2;
		blosum62[19][10]=-1;
		blosum62[19][11]=-1;
		blosum62[19][12]=-2;
		blosum62[19][13]=-1;
		blosum62[19][14]=3;
		blosum62[19][15]=-3;
		blosum62[19][16]=-2;
		blosum62[19][17]=-2;
		blosum62[19][18]=2;
		blosum62[19][19]=7;
		blosum62[19][20]=-1;
		blosum62[19][21]=-3;
		blosum62[19][22]=-2;
		blosum62[19][23]=-1;
		blosum62[20][1]=0;
		blosum62[20][2]=-3;
		blosum62[20][3]=-3;
		blosum62[20][4]=-3;
		blosum62[20][5]=-1;
		blosum62[20][6]=-2;
		blosum62[20][7]=-2;
		blosum62[20][8]=-3;
		blosum62[20][9]=-3;
		blosum62[20][10]=3;
		blosum62[20][11]=1;
		blosum62[20][12]=-2;
		blosum62[20][13]=1;
		blosum62[20][14]=-1;
		blosum62[20][15]=-2;
		blosum62[20][16]=-2;
		blosum62[20][17]=0;
		blosum62[20][18]=-3;
		blosum62[20][19]=-1;
		blosum62[20][20]=4;
		blosum62[20][21]=-3;
		blosum62[20][22]=-2;
		blosum62[20][23]=-1;
		blosum62[21][1]=-2;
		blosum62[21][2]=-1;
		blosum62[21][3]=3;
		blosum62[21][4]=4;
		blosum62[21][5]=-3;
		blosum62[21][6]=0;
		blosum62[21][7]=1;
		blosum62[21][8]=-1;
		blosum62[21][9]=0;
		blosum62[21][10]=-3;
		blosum62[21][11]=-4;
		blosum62[21][12]=0;
		blosum62[21][13]=-3;
		blosum62[21][14]=-3;
		blosum62[21][15]=-2;
		blosum62[21][16]=0;
		blosum62[21][17]=-1;
		blosum62[21][18]=-4;
		blosum62[21][19]=-3;
		blosum62[21][20]=-3;
		blosum62[21][21]=4;
		blosum62[21][22]=1;
		blosum62[21][23]=-1;
		blosum62[22][1]=-1;
		blosum62[22][2]=0;
		blosum62[22][3]=0;
		blosum62[22][4]=1;
		blosum62[22][5]=-3;
		blosum62[22][6]=3;
		blosum62[22][7]=4;
		blosum62[22][8]=-2;
		blosum62[22][9]=0;
		blosum62[22][10]=-3;
		blosum62[22][11]=-3;
		blosum62[22][12]=1;
		blosum62[22][13]=-1;
		blosum62[22][14]=-3;
		blosum62[22][15]=-1;
		blosum62[22][16]=0;
		blosum62[22][17]=-1;
		blosum62[22][18]=-3;
		blosum62[22][19]=-2;
		blosum62[22][20]=-2;
		blosum62[22][21]=1;
		blosum62[22][22]=4;
		blosum62[22][23]=-1;
		blosum62[23][1]=0;
		blosum62[23][2]=-1;
		blosum62[23][3]=-1;
		blosum62[23][4]=-1;
		blosum62[23][5]=-2;
		blosum62[23][6]=-1;
		blosum62[23][7]=-1;
		blosum62[23][8]=-1;
		blosum62[23][9]=-1;
		blosum62[23][10]=-1;
		blosum62[23][11]=-1;
		blosum62[23][12]=-1;
		blosum62[23][13]=-1;
		blosum62[23][14]=-1;
		blosum62[23][15]=-2;
		blosum62[23][16]=0;
		blosum62[23][17]=0;
		blosum62[23][18]=-2;
		blosum62[23][19]=-1;
		blosum62[23][20]=-1;
		blosum62[23][21]=-1;
		blosum62[23][22]=-1;
		blosum62[23][23]=-1;
}