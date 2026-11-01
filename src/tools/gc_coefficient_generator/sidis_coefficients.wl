(* ::Package:: *)

(* ::Text:: *)
(*WARNING: This package uses many common variable names with the assumption that they have not been defined. The code here does not follow the usual Mathematica practices regarding encapsulation.*)
(**)
(*It is therefore highly recommended to run this package in a fresh kernel.*)


(* ::Section:: *)
(*Basic utilities*)


parallelize[expr_List]:=ParallelMap[#&,Unevaluated/@Unevaluated@expr,Method->"FinestGrained"]
SetAttributes[parallelize,HoldAll];


cartesian[list1_,list2_]:=Outer[Times,list1,list2]//Flatten


f1[ft_]:=ft/2
f2[ft_,fl_]:=ft+fl


(* ::Section:: *)
(*Front-end*)


(* ::Subsection:: *)
(*Parsing FORM*)


FORMreplacement1=StringReplace[{
	"id "->"",
	"[1-xi]"->"onemxi",
	"[1-zeta]"->"onemzeta",
	"[1-xi-zeta]"->"onemximzeta",
	"[zeta-xi]"->"zetamxi",
	"[xi-zeta]"->"ximzeta",
	"[1+xi]"->"onepxi",
	"[1+zeta]"->"onepzeta",
	"(1 + xi)"->"onepxi",
	"(1 + zeta)"->"onepzeta"
}];


FORMreplacement2=StringReplace[{
	"Li2(" ~~ Shortest[expr__] ~~ ")" :> "PolyLog[2," ~~ expr ~~ "]",
	"Li3(" ~~ Shortest[expr__] ~~ ")" :> "PolyLog[3," ~~ expr ~~ "]",

	"ln("~~Shortest[expr__]~~")":>"Log["~~expr~~"]",

	"ArcTan("~~Shortest[expr__]~~")":>"ArcTan["~~expr~~"]",
	"InvTanInt("~~Shortest[expr__]~~")":>"InvTanInt["~~expr~~"]",

	"T("~~Shortest[expr__]~~")":>"T["~~expr~~"]",

	"Dd("~~Shortest[expr__]~~")":>"Dd["~~expr~~"]",
	"Dn("~~idx_~~","~~Shortest[expr__]~~")":>"Dn["~~idx~~","~~expr~~"]",

	"pi"->"\[Pi]",
	"Pi"->"\[Pi]"
}];


FORMreplacement3=StringReplace[{
	"sqrtxz1"->"sqrt1",
	"sqrtxz2"->"sqrt2",
	"sqrtxz3"->"sqrt3",
	"zeta3"->"\[Zeta]3",
	"x"->"xi",
	"z"->"zeta",
	"rln2"->"Log[2]"
}];


FORMreplacement=FORMreplacement2@*FORMreplacement1@*FORMreplacement3;


(* ::Subsection:: *)
(*Coefficient convolution*)


getExpression[coefficientExpression[expr_,coeffs_,regions_]]:=expr
getCoeffs[coefficientExpression[expr_,coeffs_,regions_]]:=coeffs
getRegions[coefficientExpression[expr_,coeffs_,regions_]]:=regions


getCoefficientVariable[coefficient[var_,expr_]]:=var
getCoefficientExpr[coefficient[var_,expr_]]:=expr

SetAttributes[#,Listable]&/@{getCoefficientVariable,getCoefficientExpr};

getRegionVariable[region[var_,expr_]]:=var
getRegionExpr[region[var_,expr_]]:=expr

SetAttributes[#,Listable]&/@{getRegionVariable,getRegionExpr};


splitByCoefficients[expr_,coeffs_,vars_]:=Module[{coeffExprs,coefficientObjs,newExpr,mask,filter},
	coeffExprs=Simplify[Coefficient[expr,coeffs]];
	mask=Boole[#=!=0]&/@coeffExprs;

	filter[coefficient[v_,e_]]:=e=!=0;

	coefficientObjs=Select[Thread[coefficient[vars,coeffExprs]],filter];

	newExpr=(mask*coeffs) . vars;
	coefficientExpression[newExpr,coefficientObjs,{}]
]


splitByCoefficients[expr_]:=Module[{coeffs},
	coeffs=cartesian[{pdf[x]/(1-x),pdf[x/\[Xi]],pdfsub[0],pdfsub[1],pdfsub[2],pdfsub[3]},{ff[z]/(1-z),ff[z/\[Zeta]],ffsub[0],ffsub[1],ffsub[2],ffsub[3]}];
	splitByCoefficients[expr,coeffs,Table[Symbol["c"<>ToString[i]],{i,1,Length[coeffs]}]]
]


evaluate[expr_]:=((splitByCoefficients[(Expand[expr*pdf[x/\[Xi]]*ff[z/\[Zeta]]]//.{
	Dd[onemxi]pdf[x/\[Xi]]->pdf[x]/(1-x),
	Dd[onemzeta]ff[z/\[Zeta]]->ff[z]/(1-z),
	Dn[k_,onemxi]pdf[x/\[Xi]]:>(pdfsub[k]),
	Dn[k_,onemzeta]ff[z/\[Zeta]]:>(ffsub[k])
})/.{
	xi->\[Xi],zeta->\[Zeta],onemxi->1-\[Xi],onemzeta->1-\[Zeta],onemximzeta->1-\[Xi]-\[Zeta],zetamxi->\[Zeta]-\[Xi],ximzeta->\[Xi]-\[Zeta],onepxi->1+\[Xi],onepzeta->1+\[Zeta],lnonepxi->Log[1+\[Xi]],lnonepzeta->Log[1+\[Zeta]]
}])/.{
	pdfsub[k_]:>(pdf[x/\[Xi]]-pdf[x]) Log[1-\[Xi]]^k/(1-\[Xi])+pdf[x] Log[1-x]^(k+1)/((1-x)(k+1)),ffsub[k_]:>(ff[z/\[Zeta]]-ff[z]) Log[1-\[Zeta]]^k/(1-\[Zeta])+ff[z] Log[1-z]^(k+1)/((1-z)(k+1))
})


(* ::Section:: *)
(*Middle-end*)


(* ::Subsection:: *)
(*Region extraction*)


extractRegions[expr_,regs_]:=Module[{coeffAssoc,filteredAssoc,coeffs,filteredRegs,variables,newExpr},
	coeffAssoc=Association@@Thread[regs->(Coefficient[expr,T[#]]&/@regs)];
	filteredAssoc=Select[coeffAssoc,#=!=0&];
	{filteredRegs,coeffs}=Comap[{Keys,Values},filteredAssoc];
	variables=(Symbol["T"<>ToString[#]])&/@filteredRegs;
	newExpr=Plus@@variables+expr/.Thread[(T/@filteredRegs)->0];
	{newExpr,Thread[region[variables,coeffs]]}
]


extractRegions[coefficientExpression[expr_,coeffs_,__],regs_:{r1,r2,t1,t2,u1,u2,u3,u4}]:=Module[{extracted,extractedExpr,extractedRules},
	{extractedExpr,extractedRules}=extractRegions[#,regs]&/@getCoefficientExpr[coeffs]\[Transpose];
	coefficientExpression[expr,Thread[coefficient[getCoefficientVariable[coeffs],extractedExpr]],Flatten[extractedRules]]
]


(* ::Subsection:: *)
(*Common transcendental extraction*)


collapseCommonTranscendentals[expr_]:=Module[{
		collectExpressions,collectPolylogs,collectLogs,collectArcTans,collectInvTanInts,generateCommonCollapseRules,
		logs,polylogs,arctans,invtanints,
		logRules,polylogRules,arctanRules,invtanintRules
	},

	collectExpressions[e_,type_]:=DeleteDuplicates[Cases[e,type,Infinity]];

	collectPolylogs[e_]:=collectExpressions[e,_PolyLog];
	collectLogs[e_]:=collectExpressions[e,_Log];
	collectArcTans[e_]:=collectExpressions[e,_ArcTan];
	collectInvTanInts[e_]:=collectExpressions[e,_InvTanInt];

	generateCommonCollapseRules[terms_,prefix_]:=Module[{vars},
		vars=Symbol/@Table[prefix<>"v"<>ToString[i],{i,1,Length[terms]}];
		Thread[terms->vars]
	];

	logs=collectLogs[expr];
	polylogs=collectPolylogs[expr];
	arctans=collectArcTans[expr];
	invtanints=collectInvTanInts[expr];

	logRules=generateCommonCollapseRules[logs,"log"];
	polylogRules=generateCommonCollapseRules[polylogs,"polylog"];
	arctanRules=generateCommonCollapseRules[arctans,"arctan"];
	invtanintRules=generateCommonCollapseRules[invtanints,"invtanint"];

	Join[logRules,polylogRules,arctanRules,invtanintRules]
]


collapseCommonTranscendentals[exprs_?ListQ]:=Module[{rules,coeffs,regions},
	coeffs=getCoeffs/@exprs;
	regions=getRegions/@exprs;
	rules=collapseCommonTranscendentals[{coeffs,regions}];
	computation[exprs/.rules,rules,{}]
]/;MatchQ[exprs,{_coefficientExpression..}]


(* ::Subsection:: *)
(*Common subexpression elimination*)


(* ::Input:: *)
(*(**)
(*Original: Stonewall Ballard *)
(*	http://stoney.sb.org/wordpress/2009/06/converting-symbolic-mathematica-expressions-to-c-code/*)
(*Modification: Eric Cousineau (eacousineau@gmail.com)*)
(*	. Reformatted code, placed it in a package*)
(*	. Fixed a small bug allowing it to return expressions with no subexpressions*)
(*	. Changed variable numbering to any number of digits*)
(*	. Fixed small bug in listSubExps[] handling if it returns no results from Scan[]*)
(**)*)


listSubExps[exp_,minLeafCount_]:=Module[{res},
Select[
Sort[
res=Reap[
Scan[
If[!AtomQ[#],Sow[#]]&,exp,Infinity
]
];
If[Length[res[[2]]]>0,res[[2,1]],{}],
LeafCount[#1]<LeafCount[#2]&
]
,LeafCount[#]>=minLeafCount&
]
];

replaceRepeats[expr_,minLeafCount_,acceptance_]:=Module[{ceqns, eqns, rules},
ceqns=Reap[
Module[{vnum=0, modExpr=expr, subExprs=listSubExps[expr,minLeafCount], subExpr, mapping,originalCount},
originalCount=Length[subExprs];
Monitor[While[subExprs!={},
subExpr = First[subExprs];
subExprs = Rest[subExprs];
If[Count[subExprs, subExpr]>0&&acceptance[subExpr],
mapping = subExpr->Symbol["vv"<>IntegerString[vnum++,10]];
Sow[mapping];
modExpr=modExpr//.mapping;
subExprs=DeleteCases[subExprs,subExpr]//.mapping
]
],
GraphicsRow[{ProgressIndicator[originalCount-Length[subExprs],{0,originalCount}],ToString[originalCount-Length[subExprs]]<>" / "<>ToString[originalCount]}]];
modExpr
]
];
eqns = ceqns[[1]];
rules = ceqns[[2]];
rules = If[Length[rules]>0,rules[[1]],{}];
Return[{eqns, rules}];
]

foldOneUseVariables[ceqns:{eqns_,rules_}]:=Module[{newRules,updates},
{newRules,updates}=Reap[
Select[rules,
If[Count[ceqns,#[[2]],Infinity]==2,Sow[#[[2]]->#[[1]]];False,True]&
]
];
{eqns,newRules}//.Flatten[updates]
]

renumberVariables[ceqns:{_,rules_}]:=Module[{vnum=0},
ceqns/.Map[#[[2]]->Symbol["v"<>IntegerString[vnum++]]&,rules]
]

hoistCommonSubexpressions[expr_,minLeafCount_,acceptance_:Function[True]]:=
(renumberVariables@*foldOneUseVariables@*replaceRepeats)[expr,minLeafCount,acceptance]


cse[computation[coefficients_,transcendentals_,_],minLeafCount_,acceptance_:Function[True]]:=Module[{coeffs,regions,expr,rules},
	coeffs=getCoeffs/@coefficients;
	regions=getRegions/@coefficients;

	{expr,rules}=hoistCommonSubexpressions[{coeffs,regions},minLeafCount,acceptance];

	computation[coefficients//.rules,transcendentals,rules]
]


(* ::Subsection:: *)
(*Intermediate representation*)


(* ::Subsubsection:: *)
(*IR utilities*)


applyReplace[f_,g_,expr_]:=f@@Thread[(g@*List)@@expr]


(* ::Subsubsection:: *)
(*IR generation*)


SetAttributes[intermediateRepresentation,Listable]
intermediateRepresentation[x_->y_]:=intermediateRepresentation[x]->intermediateRepresentation[y]
intermediateRepresentation[computation[coeffs_,transcendentals_,cse_,coeffExprs_,regions_]]:=computation[intermediateRepresentation[coeffs],intermediateRepresentation[transcendentals],intermediateRepresentation[cse],intermediateRepresentation[coeffExprs],intermediateRepresentation[regions]]
intermediateRepresentation[coefficientExpression[expr_,coeffs_,regions_]]:=coefficientExpression[intermediateRepresentation[expr],intermediateRepresentation[coeffs],intermediateRepresentation[regions]]
intermediateRepresentation[coefficient[v_,e_]]:=coefficient[intermediateRepresentation[v],intermediateRepresentation[e]]
intermediateRepresentation[region[v_,e_]]:=region[intermediateRepresentation[v],intermediateRepresentation[e]]


regionVariableMatch=StringMatchQ[RegularExpression["T(r1|r2|t1|t2|u1|u2|u3|u4)t\\d+"]];

intermediateRepresentation[expr_]:=applyReplace[IRsum,intermediateRepresentation,expr]/;Head[expr]==Plus
intermediateRepresentation[expr_]:=applyReplace[IRproduct,intermediateRepresentation,expr]/;Head[expr]==Times
intermediateRepresentation[Power[x_,n_?IntegerQ]]:=IRpow[intermediateRepresentation[x],IRinteger[n]]
intermediateRepresentation[Rational[a_?NumberQ,b_?NumberQ]]:=IRrational[IRnumber[a],IRnumber[b]]
intermediateRepresentation[Rational[a_,b_]]:=IRfraction[intermediateRepresentation[a],intermediateRepresentation[b]]

intermediateRepresentation[pdf[x]]:=IRspecial[pdf]
intermediateRepresentation[pdf[x/\[Xi]]]:=IRspecial[pdfhat]
intermediateRepresentation[ff[z]]:=IRspecial[ff]
intermediateRepresentation[ff[z/\[Zeta]]]:=IRspecial[ffhat]

intermediateRepresentation[Log[expr_]]:=IRlog[intermediateRepresentation[expr]]
intermediateRepresentation[PolyLog[n_?IntegerQ,expr_]]:=IRpolylog[IRinteger[n],intermediateRepresentation[expr]]
intermediateRepresentation[ArcTan[expr_]]:=IRarctan[intermediateRepresentation[expr]]
intermediateRepresentation[InvTanInt[expr_]]:=IRinvtanint[intermediateRepresentation[expr]]

intermediateRepresentation[reg_]:=IRT[T[Symbol[StringTake[ToString[reg],{2,3}]]],reg]/;regionVariableMatch[ToString[reg]];

intermediateRepresentation[\[Pi]]:=IRconstant[IRpi]
intermediateRepresentation[\[Zeta]3]:=IRconstant[IR\[Zeta]3]
intermediateRepresentation[x]:=IRspecial[x]
intermediateRepresentation[z]:=IRspecial[z]
intermediateRepresentation[\[Xi]]:=IRspecial[\[Xi]]
intermediateRepresentation[\[Zeta]]:=IRspecial[\[Zeta]]
intermediateRepresentation[x_]:=IRvariable[x]/;Head[x]==Symbol&&!regionVariableMatch[ToString[x]];

intermediateRepresentation[x_?NumberQ]:=IRnumber[x]

intermediateRepresentation[lnonepxi]:=intermediateRepresentation[Log[1+\[Xi]]]
intermediateRepresentation[lnonepzeta]:=intermediateRepresentation[Log[1+\[Zeta]]]


(* ::Subsubsection:: *)
(*IR simplification*)


IRproduct[IRpow[a_,n_],IRpow[b_,n_],c___]:=IRproduct[IRpow[IRproduct[a,b],n],c]


(* ::Section:: *)
(*Back-end*)


(* ::Subsection:: *)
(*Code generation*)


addParenthesis[expr_]:="("<>expr<>")"
unary[expr_]:=addParenthesis[codegen[expr]]
binary[op_,left_,right_]:="("<>codegen[left]<>op<>codegen[right]<>")"
nary[op_,exprs__]:=addParenthesis[StringJoin[Riffle[codegen/@{exprs},op]]]
function[f_,expr_]:=f<>unary[expr]
templateFunction[f_,templ_,expr_]:=f<>"<"<>codegen[templ]<>">"<>unary[expr]


codegen[IRvariable[x_]]:=SymbolName[x]
codegen[IRnumber[x_]]:=ToString[N[x]]
codegen[IRinteger[x_]]:=ToString[x]
codegen[IRrational[a_,b_]]:=binary["/",a,b]

codegen[IRsum[x__]]:=nary["+",x]
codegen[IRproduct[x__]]:=nary["*",x]
codegen[IRpow[x_,n_]]:=templateFunction["math::power",n,x]

codegen[IRlog[expr_]]:=function["std::log",expr]
codegen[IRpolylog[n_,expr_]]:=templateFunction["math::polylogT",n,expr]
codegen[IRarctan[expr_]]:=function["std::atan",expr]
codegen[IRinvtanint[expr_]]:=function["math::invtanint",expr]
codegen[IRT[_,v_]]:=codegen[IRvariable[v]]
codegen[T[expr_]]:=function["T",expr]

codegen[IRspecial[x]]:="p.pdf_momentum_fraction"
codegen[IRspecial[z]]:="p.z"
codegen[IRspecial[\[Xi]]]:="p.xi"
codegen[IRspecial[\[Zeta]]]:="p.zeta"
codegen[IRspecial[pdf]]:="p.pdf"
codegen[IRspecial[pdfhat]]:="p.pdf_hat"
codegen[IRspecial[ff]]:="p.ff"
codegen[IRspecial[ffhat]]:="p.ff_hat"

codegen[IRconstant[IRpi]]:="std::numbers::pi"
codegen[IRconstant[IR\[Zeta]3]]:="constants::zeta3"

codegen[x_]:=SymbolName[x]/;Head[x]==Symbol


(* ::Subsection:: *)
(*Destructuring*)


generateVariableName[v_]:=ToString[v]
generateVariableName[IRvariable[v_]]:=generateVariableName[v]
generateVariableName[IRT[T[_],v_]]:=generateVariableName[v]


generateVariableDeclaration[region[IRT[T[v1_,v2_]],e_]]:="const double "<>generateVariableName[v2]<>" = "<>codegen[T[v1]]<>" ? "<>codegen[e]<> " : 0.0;\n"
generateVariableDeclaration[coefficient[v_,e_]]:="const double "<>generateVariableName[v]<>" = "<>codegen[e]<>";\n"
generateVariableDeclaration[e_->v_,declarationPrefix_]:=declarationPrefix<>generateVariableName[v]<>" = "<>codegen[e]<>";\n"


generateReturnStatement[expr_]:="return "<>codegen[expr]<>";"


generateVariableBlock[expr_?ListQ]:=StringJoin[generateVariableDeclaration/@expr]


generateCode[filename_,coefficientExpression[expr_,coeffs_,regions_]]:=file[filename,generateVariableBlock[regions]<>generateVariableBlock[coeffs]<>generateReturnStatement[expr]]


generateCode[
coefficientFilenames_,updateFilename_,internalDeclarationFilename_,variableDeclarationFilename_,pto_,computation[coeffs_,transcendentals_,cse_,coeffFactors_,regions_]]:=Module[{generateCoefficientFile,coefficientFiles,updateFile,rules,internalDeclarationFile,variableDeclarationFile},

	generateCoefficientFile[filename_,coefficientExpression[expr_,cvars_,tvars_]]:=Module[{varDeclaration},
		varDeclaration[var_]:="const double "<>generateVariableName[var]<>" = p.c."<>pto<>"."<>generateVariableName[var]<>";\n";
		file[
			filename,
			StringJoin[varDeclaration/@cvars]<>generateReturnStatement[expr]
		]
	];

	coefficientFiles=MapThread[generateCoefficientFile,{coefficientFilenames,coeffs}];
	rules=Join[transcendentals,cse];

	updateFile=file[
		updateFilename,
		StringJoin[MapThread[generateVariableName[#1]<>" = "<>codegen[#2]<>";\n"&,{Values[rules],Keys[rules]}]]
		<>
		StringJoin[MapThread[generateVariableName[#1]<>" = "<>codegen[T[Symbol[StringTake[generateVariableName[#1],{2,3}]]]]<>" ? "<>codegen[#2]<> " : 0.0;\n"&,{getRegionVariable[regions],getRegionExpr[regions]}]]
		<>
		StringJoin[MapThread[generateVariableName[#1]<>" = "<>codegen[#2]<>";\n"&,{getCoefficientVariable[coeffFactors],getCoefficientExpr[coeffFactors]}]]
	];

	internalDeclarationFile=file[internalDeclarationFilename,StringJoin["[[maybe_unused]] const double "<>generateVariableName[#]<>" = p.c."<>pto<>"."<>generateVariableName[#]<>";\n"&/@Values[rules]]];
	variableDeclarationFile=file[variableDeclarationFilename,StringJoin["double "<>generateVariableName[#]<>";\n"&/@Join[Values[rules],getRegionVariable[regions],getCoefficientVariable[coeffFactors]]]];
	
	directory[
		Join[
			coefficientFiles,
			{updateFile,internalDeclarationFile,variableDeclarationFile}
		]
	]
]


(* ::Section:: *)
(*Driver*)


optimize[expr_,minLeafCount_]:=liftCTFromOptimized[
	cse[
		collapseCommonTranscendentals[
			makeCoefficientVariableNamesUnique[extractRegions/@expr]
		],
		minLeafCount,
		Function[subExpr,
			And[
				!MemberQ[{coefficient,List},Head[subExpr]],
				FreeQ[subExpr,x]
			]
		]
	]
]


makeCoefficientVariableNamesUnique[list_?ListQ]:=Module[{counter=0,renameVariables},
	renameVariables[coefficientExpression[expr_,coeffs_,regions_]]:=Module[{readCVariable,readTVariable,cvars,tvars,newCVars,newTVars,rules},
		readCVariable[coefficient[v_,e_]]:=v;
		cvars=readCVariable/@coeffs;
		newCVars=Symbol[ToString[#]<>"t"<>ToString[counter]]&/@cvars;

		readTVariable[region[v_,e_]]:=v;
		tvars=readTVariable/@regions;
		newTVars=Symbol[ToString[#]<>"t"<>ToString[counter]]&/@tvars;

		rules=Join[Thread[cvars->newCVars],Thread[tvars->newTVars]];
		counter++;
		coefficientExpression[expr,coeffs,regions]/.rules
	];
	renameVariables/@list
]


liftCTFromOptimized[computation[coefficients_,transcendentals_,cse_]]:=Module[{coeffs,regions,stripCoefficients,newCoefficients},
	coeffs=Flatten[getCoeffs/@coefficients];
	regions=Flatten[getRegions/@coefficients];
	stripCoefficients[coefficientExpression[expr_,cs_,rgs_]]:=coefficientExpression[expr,getCoefficientVariable[cs],getRegionVariable[rgs]];
	newCoefficients=stripCoefficients/@coefficients;
	computation[newCoefficients,transcendentals,cse,coeffs,regions]
]


(* ::Section:: *)
(*Exporting*)


export[directory[files_]]:=export/@files
export[file[filename_,content_]]:=Export[filename,content,"Text"]


(* ::Section:: *)
(*Reference tables*)


coefficientToRules[coefficient[var_,expr_]]:=var->expr


undoCoefficientExpression[coefficientExpression[expr_,coeffs_,regions_]]:=expr/.(coefficientToRules/@coeffs)
undoCoefficientExpression[expr_]:=expr


prepareForIntegration[expr_]:=expr/.{\[Zeta]3->Zeta[3],poly2->1+2\[Xi]+\[Xi]^2-4\[Xi] \[Zeta],sqrt1->Sqrt[1-2\[Zeta]+\[Zeta]^2+4\[Xi] \[Zeta]],sqrt2->Sqrt[1+2\[Xi]+\[Xi]^2-4\[Xi] \[Zeta]],sqrt3->Sqrt[\[Xi]/\[Zeta]],T[r1]->Boole[\[Zeta]<1-\[Xi]],T[r2]->Boole[\[Zeta]>1-\[Xi]],T[t1]->Boole[\[Zeta]>\[Xi]],T[t2]->Boole[\[Zeta]<\[Xi]],T[u1]->Boole[\[Zeta]<1-\[Xi]&&\[Zeta]<\[Xi]],T[u2]->Boole[\[Zeta]>1-\[Xi]&&\[Zeta]<\[Xi]],T[u3]->Boole[\[Zeta]<1-\[Xi]&&\[Zeta]>\[Xi]],T[u4]->Boole[\[Zeta]>1-\[Xi]&&\[Zeta]>\[Xi]],InvTanInt->(x|->1/2 I (PolyLog[2,-I Abs[x]]-PolyLog[2,I Abs[x]]) Sign[x])}


makeIntegrand[expr_,pdfF_,ffF_]:=prepareForIntegration[undoCoefficientExpression[expr]]/.{pdf->pdfF,ff->ffF}


calculateIntegralWith[expr_,pdfF_,ffF_]:=Integrate[
	makeIntegrand[expr,pdfF,ffF],
	{\[Xi],x,1},{\[Zeta],z,1},
	Assumptions->{0<x<1,0<z<1},
	GenerateConditions->False
]


pdfF=x|->x^-2 (1-x)^3;
ffF=z|->z^3 (1-z)^(3/2);


setupIntegrals[Plus[a_,b__]]:=Plus[setupIntegrals[a],setupIntegrals[Plus[b]]]
setupIntegrals[expr_]:=setupIntegrals@*prepareForIntegration@*undoCoefficientExpression@expr/;Head[expr]==coefficientExpression
setupIntegrals[expr_Times]:=Module[{factors=List@@expr,factor,otherFactor,dependentVars,independentVars,dependentFactors,independentFactors,selectIntegral},

factor[\[Xi]]=Times@@Select[factors,!FreeQ[#,pdf|x]&];
factor[\[Zeta]]=Times@@Select[factors,!FreeQ[#,ff|z]&];
otherFactor=Times@@Select[factors,FreeQ[#,pdf|x|ff|z]&];

dependentVars=Select[{\[Xi],\[Zeta]},!FreeQ[otherFactor,#]&];
independentVars=Complement[{\[Xi],\[Zeta]},dependentVars];

dependentFactors=factor/@dependentVars;
independentFactors=factor/@independentVars;

selectIntegral[e_,{}]:=e;
selectIntegral[e_,{\[Xi]}]:=integral1D[e,{\[Xi],x,1}];
selectIntegral[e_,{\[Zeta]}]:=integral1D[e,{\[Zeta],z,1}];
selectIntegral[e_,{\[Xi],\[Zeta]}]:=integralCell[e,{\[Xi],x,1},{\[Zeta],z,1}];

selectIntegral[otherFactor (Times@@dependentFactors),dependentVars]Times@@(MapThread[basisFunction[#1,#2]&,{independentFactors,independentVars}])
]


computeCells[x_,z_,upperCorner_,conditions_]:=Module[{filteredConditions,\[Zeta]solutions,\[Zeta]bounds,\[Zeta]regions,\[Xi]solutions,\[Xi]interiorSolutions,\[Xi]bounds,cells},
	filteredConditions=Cases[conditions,_Less|_Greater,{0,Infinity}];
	\[Zeta]solutions=DeleteDuplicates[
		Flatten[
			SolveValues[Equal@@#,\[Zeta]]&/@filteredConditions
		]
	];
	\[Zeta]bounds=Join[{z,upperCorner},\[Zeta]solutions];
	\[Zeta]regions=Subsets[\[Zeta]bounds,{2}];
	\[Xi]solutions=Flatten[
		SolveValues[Equal@@#,\[Xi]]&/@\[Zeta]regions
	];
	\[Xi]interiorSolutions=Select[\[Xi]solutions,x<#<upperCorner&];
	\[Xi]bounds=Union[{x,upperCorner},\[Xi]interiorSolutions];
	cells[{\[Xi]a_,\[Xi]b_}]:=With[{at=(#/.{\[Xi]->Mean[{\[Xi]a,\[Xi]b}]})&},
		{{\[Xi]a,\[Xi]b},{##}}&@@@Partition[SortBy[Select[\[Zeta]bounds,z<=at[#]<=upperCorner&], at], 2, 1]
	];
	Join@@cells/@Partition[\[Xi]bounds,2,1]
]


applyCell[expr_,{{\[Xi]a_,\[Xi]b_},{\[Zeta]a_,\[Zeta]b_}}]:=With[{\[Xi]p=Mean[{\[Xi]a,\[Xi]b}],\[Zeta]p=Mean[{\[Zeta]a,\[Zeta]b}]},
	expr/.Boole[cond_]:>Boole[cond/.{\[Xi]->\[Xi]p,\[Zeta]->(\[Zeta]p/.\[Xi]->\[Xi]p)}]
]


vectorDERule[h_,tmax_,prec_ : MachinePrecision]:=Module[{k,t,s,abscissas,weights},
	k=Range[-(tmax/h),tmax/h];
	t=N[k h,prec];
	s=\[Pi]/2 Sinh[t];
	abscissas=1/2 (1+Tanh[s]);
	weights=\[Pi]/4 h Cosh[t]Sech[s]^2;
	{abscissas,weights}
]


quadratureValues[func_,rules_,zero_]:=Quiet[Replace[(func/.rules)+zero,Except[_?NumberQ]->0,{1}]]


quadrature1D[func_,{var_,a_,b_},{u_,w_}]:=With[{values=quadratureValues[func,var->a+(b-a)u,0u]},
	(b-a)w . values
]


quadrature2Dcell[func_,{{\[Xi]a_,\[Xi]b_},{\[Zeta]a_,\[Zeta]b_}},{u_,w_}]:=Module[{\[Xi]v=\[Xi]a+(\[Xi]b-\[Xi]a)u,\[Zeta]low,\[Zeta]high,\[Zeta]width,values},
	\[Zeta]low=(\[Zeta]a/.\[Xi]->\[Xi]v)+0\[Xi]v;
	\[Zeta]high=(\[Zeta]b/.\[Xi]->\[Xi]v);
	\[Zeta]width=\[Zeta]high-\[Zeta]low;
	values=quadratureValues[func,{
		\[Xi]->Flatten[Transpose[ConstantArray[\[Xi]v,Length[u]]]], 
		\[Zeta]->Flatten[\[Zeta]low+Outer[Times,\[Zeta]width,u]]
	},0 Flatten[Outer[Times,\[Xi]v,u]]];
	(\[Xi]b-\[Xi]a) Flatten[KroneckerProduct[\[Zeta]width w,w]] . values
]


basisIntegral[func_,limits_]:=basisIntegral[func,limits]=NIntegrate[
	func,limits,
	WorkingPrecision->30,PrecisionGoal->12,MaxRecursion->50
]


integrateCells[func_,{\[Xi], \[Xi]a_, \[Xi]b_},{\[Zeta], \[Zeta]a_, \[Zeta]b_},rule_,prec_]:=Module[{conditions},
	conditions=DeleteDuplicates[
		Cases[func, Boole[cond_] :> cond, {0, Infinity}]
	];
	Total[
		quadrature2Dcell[N[applyCell[func,#],prec],#,rule]&/@computeCells[\[Xi]a, \[Zeta]a, \[Xi]b, conditions]
	]
]


calculateIntegral[expr_,{xv_,zv_},{pdfF_,ffF_},rules_,prec_:MachinePrecision]:=Module[{rule=vectorDERule[1/16,9/4,prec]},
	N[
		Block[{$MinPrecision=If[prec===MachinePrecision,0,prec],$MaxPrecision=If[prec===MachinePrecision,Infinity,prec]},
			Re[expr/.{pdf->pdfF,ff->ffF}/.rules/.{x->xv,z->zv}/. {
					basisFunction[func_,\[Xi]]:>basisIntegral[func,{\[Xi],xv,1}],basisFunction[func_,\[Zeta]]:>basisIntegral[func,{\[Zeta],zv,1}],integral1D[func_,limits_]:>quadrature1D[N[func,prec],limits,rule],integralCell[func_,\[Xi]limits_,\[Zeta]limits_]:>integrateCells[func,\[Xi]limits,\[Zeta]limits,rule,prec]
				}
			]
		]
	]
]
