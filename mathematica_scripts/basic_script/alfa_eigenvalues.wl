#!/bin/bash

math << FILE


Do[
		Print[" "];
	(*	s={1.0003,0.999,0.9964};*)
		s={0.111,1.0,1.0};
		Print["Choosen s: ",s];
	(*	e = {0.9972, 0.9998, 0.9987};*)
		e={0.999237, 0.999979, 0.987385};
		Print["Choosen e: ", e];
		Print[" "];
		n=Length[e];
		s=Sort[s];
		


		a=DiagonalMatrix[s];
		nzeros = n - Total[Sign[s]];
		sumcorrection = 0;
		If[Max[Or[e!=Abs[e], s!=Abs[s]]],Print["e s must be nonneg"],
			temp = Min[Log[Apply[Times,Sort[e]] / Apply[Times,s]]];
			If[temp < (-2*^-16)*n*Max[s], Print["majorization not satisfied"],
				For[i = 1,i < n, i++,
					posn = i;
					For[j = i,j < n, j++,
						If[a[[j,j]] <= e[[i]],
							posn = j
						];
					];
					j = posn;
					

					correction = Max[{a[[j,j]]-e[[i]], 0}];
							
					If[ correction > 0,
						Print[j];
						Print[{j, a[[j,j]], e[[i]]}];
						Print[MatrixForm[a]];
					];
					a[[j,j]] = Min[{a[[j,j]] , e[[i]]}];
					
					sumcorrection = correction + sumcorrection;
					correction = Max[{e[[i]]-a[[j+1,j+1]], 0}];
					If[ correction > 0,
						Print["j+1"];
						Print[{a[[j+1,j+1]], e[[i]]}];
					];
					a[[j+1,j+1]] = Max[a[[j+1,j+1]], e[[i]]];
					sumcorrection = correction + sumcorrection;
					perm = Append[Append[Range[1,i-1],posn],posn+1];
						For[j = i, j<n+1, j++,
							If[And[j != posn, j != posn+1],
							perm = Append[perm,j];
							];	
						]; 
					a = a[[perm,perm]]; 
					If[ e[[i]] != 0,
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
				(*Print["a"];*)
				(*Print[MatrixForm[a]];*)
				(*Print["----------------------------------------------------------"];*)
			];
			If[sumcorrection > 2*^-16*n*Max[s],
				Print["Warning", sumcorrection];
			];
				
		];	(*If temp < *)

	
	Print["a:"];
	Print[MatrixForm[a]];
	Print["========================================"];
	s=SingularValueList[a];
	Print[s];
	e=Eigenvalues[a];
	Print[e];
	]; (*if nonneg*)
]; (* first Do *)

FILE
