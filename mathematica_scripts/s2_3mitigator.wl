#!/usr/local/bin/MathematicaScript -script


expmatrix1={{1.3*10^(-3),0,0},{6.8*10^(-4),2.2*10^(-4),0},{2.7*10^(-3),1.2*10^(-3),2.8*10^(-3)}};  (*m<EW*)
expmatrix2={{2.4*10^(-2),0,0},{2.5*10^(-2),2.2*10^(-2),0},{6.9*10^(-2),1.2*10^(-2),1.0*10^(-1)}}; (*m^2>100eV*)
expmatrix3={{1.0*10^(-2),0,0},{1.7*10^(-2),1.4*10^(-2),0},{4.5*10^(-2),5.3*10^(-2),1.0*10^(-1)}}; (*m^2-1eV*)

resultmatrix = Table[Subscript[x, i, j], {i, 3}, {j, 6}];
resultmatrix[[1,1]] = "e1";
resultmatrix[[1,2]] = "e2";
resultmatrix[[1,3]] = "e3";
resultmatrix[[1,4]] = "s1";
resultmatrix[[1,5]] = "s2";
resultmatrix[[1,6]] = "s3";

(*==========Settings===============*)

expmat = 3; (* 1 for m<EW, 2 for m^2>100eV, 3 for m^2 - 1eV *)

out = 0; (* if 1 output will be printed*)

randommatricesquantity = 10^9;

decraesequantity = 1000;

loopbreak = 2; (* if number of matrices found is equal to this, then go sjump lower *)

imaginary = 1 (*if 1, table of imaginary numbers with power ip will be added to generatef real eigenvalues *)

ip = -4;   	(*power of imaginary numbers in table imag*)

singval3 = 0.999; (* starting point of s3 decraese *)

singval2 = 0.999; (* starting point of s2 decraese *)

sjump = 0.001; (* s decraese quantity *)

SeedRandom[];   (*diferrent random number generator each time based on system time*) 
(*===============================*)

Print["Settings:"];
Print["expmat: ", expmat, ", random matrices: ", randommatricesquantity, ", 'found' boundry: ", loopbreak, ", imaginary numbers power: ", ip ];
Print["s2: ", singval2, " s3: ", singval3];

If[expmat == 1,
	compmatrix={{(1-expmatrix1[[1,1]]),0,0},{expmatrix1[[2,1]],(1-expmatrix1[[2,2]]),0},{expmatrix1[[3,1]],expmatrix1[[3,2]],(1-expmatrix1[[3,3]])}},
	If[expmat == 2,
		compmatrix={{(1-expmatrix2[[1,1]]),0,0},{expmatrix2[[2,1]],(1-expmatrix2[[2,2]]),0},{expmatrix2[[3,1]],expmatrix2[[3,2]],(1-expmatrix2[[3,3]])}},
		If[expmat == 3,
			compmatrix={{(1-expmatrix3[[1,1]]),0,0},{expmatrix3[[2,1]],(1-expmatrix3[[2,2]]),0},{expmatrix3[[3,1]],expmatrix3[[3,2]],(1-expmatrix3[[3,3]])}}
		];
	];
];


control = 3;


	
	pos = 0;
	neg = 0;

	Do[
		s = {1.0,singval2,singval3};	
       		ep = {RandomReal[{compmatrix[[1,1]],1}] , RandomReal[{compmatrix[[2,2]],1}] , RandomReal[{compmatrix[[3,3]],1}]};  (*Produce random eigenvalues*)
	       	imag = { I 10^(ip), I 10^(ip), I 10^(ip)};	(*produces imaginary table*)
		If[ imaginary == 1,
       			e = ep + imag,	(*adding of random eigenvalues and imaginary tables*)
        		e = ep; 
		];
        		n=Length[e];
	        	s=Sort[s];
        		a=DiagonalMatrix[s];

        		nzeros = n - Total[Sign[s]];
	        	sumcorrection = 0;
        	        (*If[Max[Or[Abs[e]<0, s!=Abs[s]]],Print["e s must be nonneg"],*)
        		temp = Min[Log[Apply[Times,Sort[Abs[e]]] / Apply[Times,s]]];
        		If[temp < (-2*^-16)*n*Max[s], neg++,
        			For[i = 1,i < n, i++,
        				posn = i;
		        		For[j = i,j < n, j++,
        				If[a[[j,j]] <= Re[e[[i]]],
        					posn = j
        				];
	        		];
        			j = posn;				
             			correction = Max[{a[[j,j]]-Re[e[[i]]], 0}];		
	        		(*If[ correction > 0,
					Print["========================="];
        				Print[j];
        				Print[{j, a[[j,j]], Re[e[[i]]]}];
	    				For[l=1, l<n+1, l++,
    						a[[l,l]] = Abs[a[[l,l]]];
    					];
        				Print[MatrixForm[a]];
        			];*)
	        		a[[j,j]] = Min[{a[[j,j]] , Re[e[[i]]]}];
        					
        			sumcorrection = correction + sumcorrection;
        			correction = Max[{Re[e[[i]]]-a[[j+1,j+1]], 0}];
        			(*If[ correction > 0,
					(*Print["========================="];
        				Print["j+1"];
     					Print[{a[[j+1,j+1]], e[[i]]}];*)
        			];*)
	        		a[[j+1,j+1]] = Max[a[[j+1,j+1]], Re[e[[i]]]];
	        		sumcorrection = correction + sumcorrection;
	        		perm = Append[Append[Range[1,i-1],posn],posn+1];
	        		For[j = i, j<n+1, j++,
	        			If[And[j != posn, j != posn+1],
	        				perm = Append[perm,j];
	        			];	
	        		]; 
	        		a = a[[perm,perm]]; 
	        		If[ Re[e[[i]]] != 0,
	        			ee = e[[i]];
	        			s1 = a[[i+1,i+1]];
	        			s2 = a[[i,i]];
	        			U = IdentityMatrix[2];
	        			If[(s2^2 - s1^2) != 0,
	        				cu = Sqrt[((s1^2-ee^2)/(s1^2 - s2^2))];	
	        				su = Sqrt[((ee^2-s2^2)/(s1^2 - s2^2))];
	        				U = {{-cu, su},{su, cu}};
	        			];
	        		       	a[[{i,i+1}]]=U.a[[{i,i+1}]];
	        			If[ Abs[a[[i,i+1]]] > Abs[a[[i,i]]], 
	        				a[[All, {i,i+1}]]=a[[All,{i+1,i}]];
	        			];
	        			tangent = a[[i,i+1]]/a[[i,i]];
	        			cosine = 1/Sqrt[1+tangent^2];
	        			sine = cosine*tangent;
	        			V = {{cosine, sine},{-sine, cosine}};
	        			a[[All,{i,i+1}]]=a[[All,{i,i+1}]].Transpose[V];
	        			a[[i,i+1]]=0;			
	        			a[[All,i]]=Sign[a[[i,i]]]*a[[All,i]];
	        			a[[i+1,i+1]]=Abs[a[[i+1,i+1]]],
	        		
	        			If[nzeros == 1,
	        				If[Apply[Times,Sign[e[[i+1 ;; n]]]]==0,
	        					V={{0,1},{1,0}};
	        					a[[All,{i,i+1}]]=a[[All,{i,i+1}]].V;
        						a[[i+1,i+1]]=Abs[a[[i+1,i+1]]];
	        					nzeros=nzeros+1,	
	        					f = Apply[Times,e[[i+1 ;; n]]]/Apply[Times,Diagonal[a[[i+2 ;; n,i+2 ;; n]]]];
	        					a[[i+1,i]] = Sqrt[a[[i+2,i+1]]^2 - f^2];
	        					a[[i+1,i+1]] = f;
	        				];
	        			];
	        			nzeros=nzeros - 1;
	        		];
	        		indx = Ordering[Diagonal[a[[i+1 ;; n, i+1 ;; n]]]];
	        		perm = Join[Range[1,i],(indx+i*ConstantArray[1,Dimensions[indx]])];
	        		a = a[[perm,perm]];			
	        			
	        	];
	        	(*If[sumcorrection > 2*^-16*n*Max[s],
	        		Print["Warning", sumcorrection];
	        	];*)
	          	
	    		For[k=1,k<n+1,k++,
	    			a[[k,k]] = Abs[a[[k,k]]];
				
	    		];
			
			If[ out == 1,	
	        		Print["abs a:"];
	    			Print[MatrixForm[Abs[a]]];
	    			Print["compmatrix:"];
	    			Print[MatrixForm[compmatrix]];
				Print["========================================"];
			];

			If[sumcorrection < 2*^-16*n*Max[s] && correction == 0, 
				(*Print["s2: ",singval2," s3: ",singval3];
				Print["Control: ", control];*)	
	    			If[control == 3,	(*if control is equal 3:*)

					If[compmatrix[[1,1]] <= a[[1,1]] <= 1 && compmatrix[[2,2]] <= a[[2,2]] <= 1 && compmatrix[[3,3]] <= Abs[a[[3,3]]] <= 1 && Abs[a[[2,1]]] <= compmatrix[[2,1]] && Abs[a[[3,1]]] <= compmatrix[[3,1]] && Abs[a[[3,2]]] <= compmatrix[[3,2]],
						control = 2;
		  				singval3 = singval3 - sjump;
		    				resultmatrix[[2]] = Join[ep,s];
						resultmatrix[[3]] = Join[Diagonal[a],SingularValueList[a]];
						Print["==================================================================="];
						Print[MatrixForm[resultmatrix]];
						Print["-------------------------------------------------------------------"];
						Print[MatrixForm[a]];
				
		    			],	(*if control is not equal to 2:*)
	
					If[compmatrix[[1,1]] <= a[[1,1]] <= 1 && compmatrix[[2,2]] <= a[[2,2]] <= 1 && compmatrix[[3,3]] <= Abs[a[[3,3]]] <= 1 && Abs[a[[2,1]]] <= compmatrix[[2,1]] && Abs[a[[3,1]]] <= compmatrix[[3,1]] && Abs[a[[3,2]]] <= compmatrix[[3,2]],
						control = 3;
						pos++;
		  				singval2 = singval2 - sjump;
		    				resultmatrix[[2]] = Join[ep,s];
						resultmatrix[[3]] = Join[Diagonal[a],SingularValueList[a]];
						Print["==================================================================="];
						Print[MatrixForm[resultmatrix]];
						Print["-------------------------------------------------------------------"];
						Print[MatrixForm[a]];
				
		    			];
				];	(*If control*)
			]; (*If sumcorrection*)

	        ];	(*If temp < *)        	
        	(*]; (*if nonneg*)*)

,{randommatricesquantity}];

Print["Could not find any matrix for s2: ", singval2, " s3: ", singval3, " - break"];
	
   


