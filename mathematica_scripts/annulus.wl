#!/usr/local/bin/MathematicaScript -script


expmatrix3={{1.0*10^(-2),0,0},{1.7*10^(-2),1.4*10^(-2),0},{4.5*10^(-2),5.3*10^(-2),1.0*10^(-1)}}; (*m^2-1eV*)
compmatrix={{(1-expmatrix3[[1,1]]),0,0},{expmatrix3[[2,1]],(1-expmatrix3[[2,2]]),0},{expmatrix3[[3,1]],expmatrix3[[3,2]],(1-expmatrix3[[3,3]])}}


Annulus /:
	Random`DistributionVector[
		Annulus[min_?Positive, max_?Positive], n_Integer, prec_?Positive] :=
		RandomReal[{min,max}, n, WorkingPrecision -> prec] * Transpose[{Cos[#], Sin[#]} & [RandomReal[ 2 Pi, n, WorkingPrecision -> prec ]]];

ep={1,2,3};

For[i=1, i < 4, i++,
	gene = RandomReal[Annulus[compmatrix[[i,i]],1]];
	z = gene[[1]] + I gene[[2]];
	ep[[i]] = z;
];

Print[ep];
