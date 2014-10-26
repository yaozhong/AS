this is the developed version for weighted EM algorithms for RNA-seq
the current version is very slow and far from being used.

We will refine the code from the following aspects:
(1). multi-core process for each transcripts
(2). other optimization 

2014/10/26 successfully applied mclapply to transcripts inside of each chromesome,
	       the speed is aroud 3~4x depends the core numbers and computer
		   e.g., in chromesome 21, K562 data, 1-core:9min, 7-core:3min)
		   the speed also depeends on the balance of splitted data. 

"to apply mclapply" the data should be totally indpendent, no shared varible access.

