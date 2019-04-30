#!/bin/bash

math << FILE








Do[

		s={1.0003 ,0.999, 0.9964};
		Print["Choosen s: ",s];
		e = {0.9972, 0.9998, 0.9987};
		Print["Choosen e: ", e];

		n=Length[e];
		s=Sort[s];
		Print[(-2*^-16)*n*Max[s]];


		a=DiagonalMatrix[s];
		nzeros = n - Total[Sign[s]];
		Print["nzeros: ", nzeros];
		sumcorrection = 0;
		Print["sumcorrection: ", sumcorrection];
		If[Max[Or[e!=Abs[e], s!=Abs[s]]],Print["e s must be nonneg"],
			temp = Min[Log[Apply[Times,Sort[e]] / Apply[Times,s]]];
			Print["temp: ", temp];
			If[temp < (-2*^-16)*n*Max[s], Print["majorization not satisfied"], 
				For[i = 1,i < n, i++,
					posn = i;
					Print["posn, pętla i: ", posn];
					For[j = i,j < n, j++,
						If[a[[j,j]] <= e[[i]],
							posn = j
						];
						Print["posn, pętla j: " ,posn];
					];
					j = posn;
					Print["j = posn: ", j];

					correction = Max[{a[[j,j]]-e[[i]], 0}];
					Print["correction: ", correction];
		
					If[ correction > 0,
						Print[j];
						Print[{j, a[[j,j]], e[[i]]}];
						Print[a];
					];
					a[[j,j]] = Min[{a[[j,j]] , e[[i]]}];
					Print["a"];
					Print[MatrixForm[a]];	
					sumcorrection = correction + sumcorrection;
					Print["sumcorrection ", sumcorrection];
					correction = Max[{e[[i]]-a[[j+1,j+1]], 0}];
					Print["correction: ", correction];
					If[ correction > 0,
						Print["j+1"];
						Print[{a[[j+1,j+1]], e[[i]]}];
					];
					a[[j+1,j+1]] = Max[a[[j+1,j+1]], e[[i]]];  (*błąd*)
					Print["a"];
					Print[MatrixForm[a]];
					sumcorrection = correction + sumcorrection;
					Print["sumcorrection: ", sumcorrection];
					perm = Append[Append[Range[1,i-1],posn],posn+1];
					Print["perm: ", perm];
						For[j = i, j<n+1, j++,
							If[And[j != posn, j != posn+1],
							perm = Append[perm,j];
							Print["appended perm: ", perm];
							];	
						]; 
					a = a[[perm,perm]];
					Print["a"];
					Print[MatrixForm[a]]; 
					If[ e[[i]] != 0,
						ee = e[[i]];
						Print["ee: ", ee];
						s1 = a[[i+1,i+1]];
						Print["s1:", s1];
						s2 = a[[i,i]];
						Print["s2: ", s2];
						U = IdentityMatrix[2];
						Print["U"];
						Print[MatrixForm[U]];
						If[(s2^2 - s1^2) != 0,
							cu = Sqrt[((s1^2-ee^2)/(s1^2 - s2^2))];	
							Print["cu: ", cu];
							su = Sqrt[((ee^2-s2^2)/(s1^2 - s2^2))];
							Print["su: ", su];
							U = {{-cu, su},{su, cu}};
							Print["U"];
							Print[MatrixForm[U]];
						];
				        	a[[{i,i+1}]]=U.a[[{i,i+1}]];
						Print["a"];
						Print[MatrixForm[a]];
						If[ Abs[a[[i,i+1]]] > Abs[a[[i,i]]], 
							a[[All, {i,i+1}]]=a[[All,{i+1,i}]];	
							Print["a"];
							Print[MatrixForm[a]];
						];
						tangent = a[[i,i+1]]/a[[i,i]];
						Print["tangent"];
						Print[tangent];
						cosine = 1/Sqrt[1+tangent^2];
						Print["cosine: ", cosine];
						sine = cosine*tangent;
						Print["sine:", sine];
						V = {{cosine, sine},{-sine, cosine}};			
						Print["V"];
						Print[MatrixForm[V]];
						a[[All,{i,i+1}]]=a[[All,{i,i+1}]].Transpose[V];
						Print["a"]
						Print[MatrixForm[a]];
						a[[i,i+1]]=0;
						Print["a"];
						Print[MatrixForm[a]];				
						a[[All,i]]=Sign[a[[i,i]]]*a[[All,i]];	
						Print["a"];
						Print[MatrixForm[a]];
						a[[i+1,i+1]]=Abs[a[[i+1,i+1]]];
						Print["a"];
						Print[MatrixForm[a]],
			
						If[nzeros ==1,
							If[Apply[Times,Sign[e[[i+1 ;; n]]]]==0,
								V={{0,1},{1,0}};
								Print["V"];
								Print[MatrixForm[V]];
								a[[All,{i,i+1}]]=a[[All,{i,i+1}]].V;
								Print["a"];
								Print[MatrixForm[a]];
								a[[i+1,i+1]]=Abs[a[[i+1,i+1]]];
								Print["a"];
								Print[MatrixForm[a]];
								nzeros=nzeros+1,	
								Print["nzeros: ", nzeros];
								f = Apply[Times,e[[i+1 ;; n]]]/Apply[Times,Diagonal[a[[i+2 ;; n,i+2 ;; n]]]];
								Print["f: ", f];
								a[[i+1,i]] = Sqrt[a[[i+2,i+1]]^2 - f^2];
								Print["a"];
								Print[MatrixForm[a]];
								a[[i+1,i+1]] = f;
								Print["a"];
								Print[MatrixForm[a]];
							];
						];
						nzeros=nzeros-1;
						Print["nzeros: ", nzeros];
					];
				indx = Ordering[Diagonal[a[[i+1 ;; n, i+1 ;; n]]]];
				Print["indx: ", indx];
				perm = Join[Range[1,i],(indx+i*ConstantArray[1,Dimensions[indx]])];
				Print["perm: ", perm];
				a = a[[perm,perm]];			
				Print["last a"];
				Print[MatrixForm[a]];
				Print["=========================================================="];
			];
			If[sumcorrection > 2*^-16*n*Max[s],
				Print["Warning", sumcorrection];
			];
				
		];	(*If temp < *)

	]; (* If nonneg*)

]; (* first Do *)

FILE
