#NW-align 

By using with BLOSUM62 matrix and Needleman-Wunsch dynamical programming, this program can be used to align two protein sequences.

This program, which is wrote by using C/C++ lauguage. 

You can download the C/C++ source code and compile it by yourself, as follows:

      >g++ -o NWalign NWalign.cpp
      
      or
      
      >g++ -static -O3 -ffast-math -lm -o NWalign NWalign.cpp

Here, an example is given as:

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
