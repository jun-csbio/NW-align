#NW-align 

By using with BLOSUM62 matrix and Needleman-Wunsch dynamical programming, this program can be used to align two protein sequences.

This program, which is wrote by using C/C++ lauguage, can implement the similar function of the web page (https://zhanglab.ccmb.med.umich.edu/NW-align). 

On 10,000 protein sequence pairs, which are randomly selected from PDB, we compared the running time of this program and the fortran version of NW-align, which is downloaded from https://zhanglab.ccmb.med.umich.edu/NW-align. The comparison results show this C/C++ program (average running time is 0.0032 s) can run faster than the fortran version (average running time is 0.0075 s). 

You can download the C/C++ source code and compile it by yourself, as follows:

      >g++ -static -O3 -ffast-math -lm -o NWalign NWalign.cpp
      
      or
      
      >g++ -o NWalign NWalign.cpp

Here, an example of usage is given as:

      >./NWalign
        
            query : Query sequence
            templ : Template sequence
            
            e.g.,:
            ./NWalign GKVFLTNAFSINMLKEFPTTITIDKLDEEDFCLKLELRLEDGTLINAIGHDSTINLVNTL VQGAGVVGETPTIPPNTAYQYTSGTVLDTPFGIMYGTYGMVSESGEHFNAIIKPFRLATP
      
      or
      
      >./NWalign GKVFLTNAFSINMLKEFPTTITIDKLDEEDFCLKLELRLEDGTLINAIGHDSTINLVNTL VQGAGVVGETPTIPPNTAYQYTSGTVLDTPFGIMYGTYGMVSESGEHFNAIIKPFRLATP
      
            First (query) sequence length     : 60
            Second (template) sequence length : 60
            Number of identical pairs : 12
            Number of aligned pairs : 53
            Sequence identity : 0.2 (12/60)
      
            ----GKVFLTNAFSINMLKEFPTTITIDKLDEEDFCLKLE---LRLEDGTLINAIGHDSTINLVNTL
                : :  :     :           :      :           : :   :::          : 
            VQGAGVVGETPTIPPNTAYQYTSGTVLDT----PFGIMYGTYGMVSESGEHFNAIIKPFRLA---TP



Contact: Jun Hu (jun_cs@126.com)
