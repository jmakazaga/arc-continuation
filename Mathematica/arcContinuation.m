(* ::Package:: *)

Clear[ArcContinuation,SetArcContinuation,ReLinkArcContinuation,SetArcContinuationExample,
ArcContinuationDir, AuxiliaryRules,OSType,RealJacobian,RealRest,RealFindRoot];

{DoubleContinuation,SetDoubleContinuation,ReLinkDoubleContinuation,SetDoubleContinuationExample,DoubleJacobian,DoubleRest,DoubleFindRoot} = 
{ArcContinuation,SetArcContinuation,ReLinkArcContinuation,SetArcContinuationExample,RealJacobian,RealRest,RealFindRoot};



Options[SetArcContinuation] = 
{AuxiliaryRules->{}, ArcContinuationDir->"~/ik/homotopia_double",
               TargetDir->Automatic, OSType->"Linux", RealParameters->{}, IntegerParameters->{},
               HomogeneousVariables->False};

SetArcContinuation[rest_List, t_Symbol, vars_List, opts___Rule]:=
Module[{CCDir,ostype,auxrules,exampleargs, DDrulesall,DDrules0,DDrules,
        parsrules,auxvarsrules, varsrules,dauxvarsrules,allrules,denbora, 
      indaux,ind1,ind2,makefile,(*resta,restb,*)targetdir,fileconf,varsymbollist,makecmd},
(*denbora = SessionTime[];*)
Clear[DDrule,hvarsQ,simplrulesfcn,preal,aux,daux,var,cvar,auxvars,dauxvars,
      resaux,auxaux,dauxaux,jacaux,dir];
CCDir = ArcContinuationDir/.{opts}/.Options[SetArcContinuation];
targetdir = TargetDir/.{opts}/.Options[SetArcContinuation];
If[targetdir===Automatic, targetdir = CCDir <> "/temp"];
If[!DirectoryQ[targetdir],CreateDirectory[targetdir]];
ostype = OSType/.{opts}/.Options[SetArcContinuation];
hvarsQ = HomogeneousVariables/.{opts}/.Options[SetArcContinuation];

varsymbollist =Apply[Union,Map[Cases[#,_Symbol,Infinity,Heads->True]&,vars]];
Map[SetAttributes[#,NHoldAll]&,varsymbollist];

auxrules = AuxiliaryRules/.{opts}/.Options[SetArcContinuation];
resta = rest;
restb = rest;
splicevars = vars;
If[hvarsQ,
resta = Prepend[resta, 0.];
restb = Prepend[restb, (vars . Array[cvar,Length[vars],0])-1];
auxrules = Join[auxrules,
                Array[cvar[#]->vars[[#+1]]&,Length[vars],0]];
(*auxrules = Join[auxrules,
                Array[cvar[#]->Conjugate[vars[[#+1]]]&,Length[vars],0]];*)
  ];
pars = RealParameters/.{opts}/.Options[SetArcContinuation];
ipars = IntegerParameters/.{opts}/.Options[SetArcContinuation]; (* honetarako \[Integral]user_fcns.mc fitxategia prestatu gabe dago oraindik*)
If[ostype==="Mac", makefile = CCDir <> "/Makefile.mac",
				    makefile = CCDir <> "/Makefile"
                   (*If[ostype==="linux32", makefile = CCDir <> "/Makefile.linux32",
                                           makefile = CCDir <> "/Makefile.linux-eli"
                      ]*)
 ];

auxvars = Map[First,auxrules];

DDrule[a_,b_]:= DD[a,b]->(D[a/.auxrules,b]+Map[D[a/.auxrules,#]&,auxvars].Map[DD[#,b]&,auxvars]);
DDrule[cvar[j_],b_]:= DD[cvar[j],b]->0;

simplrulesfcn[{rules0_,rules1_}]:=
     Module[{rules00},
     rules00=Select[rules1,NumberQ[Last[#1]]||MemberQ[Last[#1],vars]||
                           MemberQ[Last[#1],auxvars]&];
     {Join[rules0,rules00],Complement[rules1,rules00]/. rules00}
           ];

(*Print["DDrules definitzera goaz. Length[auxvars]*Length[vars]=",Length[auxvars]*Length[vars]];
Print["Orain arteko denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)

DDrulesall=Flatten[Outer[DDrule,auxvars,vars]];

(*Print["DDrulesall definitu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)


(*{DDrules0,DDrules}=FixedPoint[simplrulesfcn,{{},DDrulesall}];*)
DDrules = DDrulesall;
DDrules0 ={};

(*Print["DDrules definitu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)

dauxvars = Map[First,DDrules];
parsrules = Thread[pars->Array[HoldForm[preal[[#]]]&,Length[pars],0]];
auxvarsrules= Thread[auxvars->Array[HoldForm[aux[[#]]]&,Length[auxvars],0]];
varsrules = Thread[vars->Array[HoldForm[var[[#]]]&,Length[vars],0]];
dauxvarsrules = Thread[dauxvars->Array[HoldForm[daux[[#]]]&,Length[dauxvars],0]];
allrules = Join[varsrules,parsrules,auxvarsrules,dauxvarsrules,
                {t -> HoldForm[t]}];
auxaux=Array[SequenceForm[CForm[HoldForm[aux[[#]]]],
 " = ",
CForm[auxrules[[#+1,2]]/.allrules],
";"]&,Length[auxvars],0];

resaux=Array[SequenceForm[CForm[HoldForm[res[[#]]]],
 "=",
CForm[resta[[#+1]]/.allrules],
";"]&,Length[resta],0];


(*Print["resaux definitu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)

dauxaux=Array[SequenceForm[CForm[HoldForm[daux[[#]]]],
 " = ",
CForm[DDrules[[#+1,2]]/.allrules],
";"]&,Length[dauxvars],0];

indaux=Flatten[Table[{j,i},{i,Length[vars]},{j,Length[vars]}],1];
ind1 = Map[First,indaux];
ind2 = Map[Last,indaux];

(*Print["indaux definitu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)


jacaux=Array[{SequenceForm["//",CForm[HoldForm[jac[[#]]]],
 " = ","D[res[",ind1[[#+1]]-1,"],","var[" ,ind2[[#+1]]-1,"]]"],
              SequenceForm[CForm[HoldForm[jac[[#]]]],
 " = ",CForm[D[restb[[ind1[[#+1]]]],vars[[ind2[[#+1]]]]]+
             Sum[D[restb[[ind1[[#+1]]]],auxvars[[j]] ]DD[auxvars[[j]],vars[[ind2[[#+1]] ]] ],
                 {j,Length[auxvars]}
                ]/.DDrules0/.allrules],";"]}&,
              Length[vars]^2,0]//Flatten;

(*Print["jacaux definitu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)

dir = Directory[];
SetDirectory[targetdir];
Splice[CCDir <>"/Mathematica/user_fcns.mc","user_fcns.c",FormatType->OutputForm];


fileconf = "./arcContinuation.conf";
If[!FileExistsQ[fileconf], Run[ "cp " <> CCDir <> "/arcContinuation.conf" <> " ./"]];
(*Run["rm emaitzak.m"];*)
makecmd = "echo \"VERSION="<>ToString[$VersionNumber]<>"0\" > Makefile";
Run[makecmd];
makecmd = "echo \"SYS="<>$SystemID<>"\" >> Makefile";
Run[makecmd];
makecmd = "echo \"MHDIR="<>CCDir<>"\" >> Makefile";
Run[makecmd];
Run[ "cat " <> makefile <> ">> Makefile"];
(*
Run[ "cp " <> makefile <> " Makefile"];
makecmd = "make "<>"VERSION="<>ToString[$VersionNumber]<>"0 SYS="<>$SystemID<>" MHDIR="<>CCDir;
Print[makecmd];
Run[makecmd];
*)
Run["make clean"]; 
Run["make"];



(*Print["konpilazioa bukatu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)

esteka=Install["mathArcContinuation"];
SetDirectory[dir];

(*Print["math-homotopiaren instalazioa bukatu degu"];
Print["Tarte honetako denbora=",SessionTime[]-denbora];
denbora = SessionTime[];*)

esteka
];




ReLinkArcContinuation[opts___Rule]:=
Module[{CCDir,targetdir,dir,esteka,makecmd},
CCDir = ArcContinuationDir/.{opts}/.Options[SetArcContinuation];
targetdir = TargetDir/.{opts}/.Options[SetArcContinuation];
If[targetdir===Automatic, targetdir = CCDir <> "/temp"];
If[!DirectoryQ[targetdir],CreateDirectory[targetdir]];
dir = Directory[];
SetDirectory[targetdir];
esteka=Install["mathArcContinuation"];
SetDirectory[dir];
esteka
];


Options[SetArcContinuationExample] = {ArcContinuationDir->"~/ownCloud/ik/github/arc-continuation", TargetDir->Automatic,
                          RealParameters->{}, IntegerParameters->{}};


SetArcContinuationExample[{var0__}, {beta0__},t0_?NumberQ, tend_?NumberQ, action_?NumberQ, opts___]:=
 (*Module[{CCDir,hvarsQ,(*vars0,rt0,rtend,*)dir,targetdir(*,rpars0,ipars0*)},
Clear[rpars0,ipars0];*)
 Block[{CCDir,hvarsQ,vars0,betas0,action0,rt0,rtend,dir,targetdir,rpars0,ipars0},
CCDir = ArcContinuationDir/.{opts}/.Options[SetArcContinuation];
targetdir = TargetDir/.{opts}/.Options[SetArcContinuation];
If[targetdir===Automatic, targetdir = CCDir <> "/temp"];
hvarsQ = HomogeneousVariables/.{opts}/.Options[SetArcContinuation];
vars0 = N[{var0}]+ 0.;
If[hvarsQ,vars0=vars0/Norm[vars0]];
betas0=N[{beta0}]+ 0.;
rt0 = Re[N[t0]];
rtend = Re[N[tend]];
action0=Round[action];
rpars0 = Re[N[RealParameters/.{opts}/.Options[SetArcContinuationExample]]];
ipars0 = IntegerParameters/.{opts}/.Options[SetArcContinuationExample];
dir = Directory[];
SetDirectory[targetdir];
Run["make terminalArcContinuation"];
SetDirectory[dir];
Splice[CCDir <>"/Mathematica/arcContinuation.minit",targetdir <>"/arcContinuation.init",
      FormatType->OutputForm]
];



(* ::Input:: *)
(**)


RealRest[{var0__},t0_?NumberQ,opts___]:= 
 Module[{vars0,rt0,rpars0,ipars0},
vars0 = N[{var0}];
rt0 = Re[N[t0]];
rpars0 = Re[N[RealParameters/.{opts}/.Options[ArcContinuation]]];
ipars0 = Round[IntegerParameters/.{opts}/.Options[ArcContinuation]];
realRest[vars0,rt0,rpars0,ipars0]//First];


RealJacobian[{var0__},t0_?NumberQ,opts___]:= 
 Module[{vars0,rt0,rpars0,ipars0},
vars0 = N[{var0}];
rt0 = Re[N[t0]];
rpars0 = Re[N[RealParameters/.{opts}/.Options[ArcContinuation]]];
ipars0 = Round[IntegerParameters/.{opts}/.Options[ArcContinuation]];
Transpose[realJacobian[vars0,rt0,rpars0,ipars0]]
];


Options[RealFindRoot] = {RealParameters->{},IntegerParameters->{},RelativePrecisionGoal->10, AbsolutePrecisionGoal->16, 
                            NumericalJacobian->True};

RealFindRoot[{var0__},t0_?NumberQ,opts___]:= 
 Module[{vars0,rt0,rpars0,ipars0,acc,prec,njacQ},
vars0 = N[{var0}];
rt0 = Re[N[t0]];
rpars0 = Re[N[RealParameters/.{opts}/.Options[RealFindRoot]]];
ipars0 = Round[IntegerParameters/.{opts}/.Options[RealFindRoot]];
acc = N[RelativePrecisionGoal/.{opts}/.Options[RealFindRoot]];
prec = N[AbsolutePrecisionGoal/.{opts}/.Options[RealFindRoot]];
njacQ = N[NumericalJacobian/.{opts}/.Options[RealFindRoot]];
realFindroot[vars0,rt0,rpars0,ipars0,acc,prec,If[njacQ,1,0,0]]
];

Options[ArcContinuation] = {RealParameters->{},IntegerParameters->{}};



ArcContinuation[{var0__},{beta0__},t0_?NumberQ,tend_?NumberQ, action_?NumberQ,opts___]:=
 Module[{vars0,betas0,rt0,rtend,rpars0,action0,ipars0,hvarsQ},
hvarsQ = HomogeneousVariables/.{opts}/.Options[SetArcContinuation];
vars0 = N[{var0}];
If[hvarsQ,vars0=vars0/Norm[vars0]];
betas0 = N[{beta0}];
rt0 = Re[N[t0]];
rtend = Re[N[tend]];
action0 = Round[action];
rpars0 = Re[N[RealParameters/.{opts}/.Options[ArcContinuation]]];
ipars0 = Round[IntegerParameters/.{opts}/.Options[ArcContinuation]];
arcContinuation[vars0,betas0,rt0,rtend,action0,rpars0,ipars0]
];



