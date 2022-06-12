(* ::Package:: *)

BeginPackage["Rotation`"]


PolarForm::usage="";
sph::usage="";
ToQuaternion::usage="";
ToMatrix::usage="";
CreateQuaternion::usage="";
UnesiMkvaterniona::usage="";
Rotation::usage="";
RotationAveraging::usage="";
InterpolationTwoQuat::usage="";
SimulationConsensus::usage="";
SimulationInterpolation::usage="";
VisualizationConsensus::usage="";
VisualizationInterpolation::usage="";
DistribSync::usage="";
SimulationDistribSync::usage="";
VisualizationDistribSync::usage="";


Begin["`Private`"]


PolarForm[q_]:={ArcCos[q [[1]]],{q [[2]]/Sin[ArcCos[q [[1]]]],q [[3]]/Sin[ArcCos[q [[1]]]],q [[4]]/Sin[ArcCos[q [[1]]]]}}


sph[\[Theta]_,\[Phi]_]:={ Cos[\[Phi]] Sin[\[Theta]],Sin[\[Phi]] Sin[\[Theta]],Cos[\[Theta]]}


ToQuaternion[matrix_] := {Re[matrix[[1]][[1]]], Im[matrix[[1]][[1]]], Re[matrix[[1]][[2]]],Im[matrix[[1]][[2]]]}


ToMatrix[quat_]:={{quat[[1]]+I quat[[2]],quat[[3]]+I quat[[4]]},{-quat[[3]]+I quat[[4]],quat[[1]]-I quat[[2]]}}


CreateQuaternion[rotac_,phi_,theta_]:={Cos[rotac],sph[phi,theta] [[1]]* Sin[rotac],sph[phi,theta][[2]]* Sin[rotac],sph[phi,theta][[3]]* Sin[rotac]}


UnesiMkvaterniona[n_]:=Table[CreateQuaternion[Input[],Input[],Input[]],{i,1,n}]


Rotation[lista_]:=Module[{list=lista,broj,g},
broj=Length[list];
Do[Subscript[g,i]=Graphics3D[GeometricTransformation[ExampleData[{"Geometry3D","SpaceShuttle"},"GraphicsComplex"],RotationTransform[2PolarForm[list[[i]]][[1]],PolarForm[list[[i]]][[2]]]]],{i,1,broj}];
Return[Grid[ArrayReshape[Table[{Subscript[g,i],Style[Framed[Chop[list[[i]][[1]]]+Chop[list[[i]][[2]]]"\[DoubleStruckI]"+Chop[list[[i]][[3]]]"\[DoubleStruckJ]"+Chop[list[[i]][[4]]]"\[DoubleStruckK]"],10,Blue]},{i,1,broj}],{broj,2}]]]]


RotationAveraging[ListUnitQuat_,Coupling_]:=Module[{K=Coupling,eqns},
Clear[n];
n=Length[ListUnitQuat];
eqns=Table[{Subscript[x1,i]'[t]==-2Subscript[x1,i][t] Subscript[x2,i][t] (K /(2n) Sum[Subscript[x2,j][t],{j,1,n}])-2Subscript[x1,i][t] Subscript[x3,i][t] (K /(2n) Sum[Subscript[x3,j][t],{j,1,n}])-2Subscript[x1,i][t] Subscript[x4,i][t] (K /(2n) Sum[Subscript[x4,j][t],{j,1,n}])+(Subscript[x1,i][t]^2-Subscript[x2,i][t]^2-Subscript[x3,i][t]^2-Subscript[x4,i][t]^2-1)(-(K /(2n))Sum[Subscript[x1,j][t],{j,1,n}]),Subscript[x2,i]'[t]==2Subscript[x1,i][t] Subscript[x2,i][t] (-(K /(2n))Sum[Subscript[x1,j][t],{j,1,n}])-2Subscript[x2,i][t] Subscript[x3,i][t] (K /(2n) Sum[Subscript[x3,j][t],{j,1,n}])-2Subscript[x2,i][t] Subscript[x4,i][t] (K /(2n) Sum[Subscript[x4,j][t],{j,1,n}])+(Subscript[x1,i][t]^2-Subscript[x2,i][t]^2+Subscript[x3,i][t]^2+Subscript[x4,i][t]^2+1)(K /(2n) Sum[Subscript[x2,j][t],{j,1,n}]),
Subscript[x3,i]'[t]==2Subscript[x1,i][t] Subscript[x3,i][t] (-(K /(2n))Sum[Subscript[x1,j][t],{j,1,n}])-2Subscript[x2,i][t] Subscript[x3,i][t] (K /(2n) Sum[Subscript[x2,j][t],{j,1,n}])-2Subscript[x3,i][t] Subscript[x4,i][t] (K/(2n) Sum[Subscript[x4,j][t],{j,1,n}])+(Subscript[x1,i][t]^2+Subscript[x2,i][t]^2-Subscript[x3,i][t]^2+Subscript[x4,i][t]^2+1)(K /(2n) Sum[Subscript[x3,j][t],{j,1,n}]),
Subscript[x4,i]'[t]==2Subscript[x1,i][t] Subscript[x4,i][t] (-(K /(2n))Sum[Subscript[x1,j][t],{j,1,n}])-2Subscript[x2,i][t] Subscript[x4,i][t] (K /(2n) Sum[Subscript[x2,j][t],{j,1,n}])-2Subscript[x3,i][t] Subscript[x4,i][t] (K /(2n) Sum[Subscript[x3,j][t],{j,1,n}])+(Subscript[x1,i][t]^2+Subscript[x2,i][t]^2+Subscript[x3,i][t]^2-Subscript[x4,i][t]^2+1)(K /(2n) Sum[Subscript[x4,j][t],{j,1,n}]),
Subscript[x1,i][0]==ListUnitQuat[[i]][[1]],Subscript[x2,i][0]==ListUnitQuat[[i]][[2]],Subscript[x3,i][0]==ListUnitQuat[[i]][[3]],Subscript[x4,i][0]==ListUnitQuat[[i]][[4]]},{i,1,n}];
sol=NDSolve[eqns,Flatten[{Table[Subscript[x1,i][t],{i,1,n}],Table[Subscript[x2,i][t],{i,1,n}],Table[Subscript[x3,i][t],{i,1,n}],Table[Subscript[x4,i][t],{i,1,n}]}],{t,0,100},Method->{"EquationSimplification"->"Residual"}];
Do[Subscript[Quat,i][t_]:=Evaluate[{{Subscript[x1,i][t]/.sol[[1]],Subscript[x2,i][t]/.sol[[1]],Subscript[x3,i][t]/.sol[[1]],Subscript[x4,i][t]/.sol[[1]]}}],{i,1,n}];
Return[Subscript[Quat,1][100]]]


DistribSync[brojkvat_,coupling_,timecoup_]:=Module[{K=coupling,aa1,aa2,aa3,q1,q2,q3,q4,qq1,qq2,qq3,qq4,w,eqns3,eqns4},
tc=timecoup;
nn=brojkvat;
Do[Subscript[Q,i_][t_]:={{Subscript[q1,i][t]+I Subscript[q2,i][t],-Subscript[q3,i][t]+I Subscript[q4,i][t]},{Subscript[q3,i][t]+I Subscript[q4,i][t],Subscript[q1,i][t]-I Subscript[q2,i][t]}},{i,1,nn}];
Do[aa1=RandomReal[];aa2=RandomReal[];aa3=RandomReal[];
Subscript[w,i]={{aa1 I,aa2 +aa3  I},{-aa2 +aa3  I,-aa1  I}},{i,1,nn}];
eqns3=Table[{Subscript[Q,i]'[t]==Subscript[w,i].Subscript[Q,i][t],Subscript[Q,i][0]=={{1,0},{0,1}}},{i,1,nn}];
sol3=NDSolve[eqns3,Flatten[{Table[Subscript[q1,i][t],{i,1,nn}],Table[Subscript[q2,i][t],{i,1,nn}],Table[Subscript[q3,i][t],{i,1,nn}],Table[Subscript[q4,i][t],{i,1,nn}]}],{t,0.01,20}];
Do[Subscript[fun2,i][t_]:=Evaluate[Subscript[Q,i][t]/.sol3],{i,1,nn}];
Do[Subscript[QQ,i_][t_]:={{Subscript[qq1,i][t]+I Subscript[qq2,i][t],-Subscript[qq3,i][t]+I Subscript[qq4,i][t]},{Subscript[qq3,i][t]+I Subscript[qq4,i][t],Subscript[qq1,i][t]-I Subscript[qq2,i][t]}},{i,1,nn}];
eqns4=Table[{Subscript[QQ,i]'[t]==Subscript[QQ,i][t].(-(K/(2nn))Sum[ConjugateTranspose[Subscript[QQ,j][t]],{j,1,nn}]).Subscript[QQ,i][t]+Subscript[w,i].Subscript[QQ,i][t]-ConjugateTranspose[-(K/(2nn))Sum[ConjugateTranspose[Subscript[QQ,j][t]],{j,1,nn}]],Subscript[QQ,i][0]==Subscript[fun2,i][tc][[1]]},{i,1,nn}];
sol4=NDSolve[eqns4,Flatten[{Table[Subscript[qq1,i][t],{i,1,nn}],Table[Subscript[qq2,i][t],{i,1,nn}],Table[Subscript[qq3,i][t],{i,1,nn}],Table[Subscript[qq4,i][t],{i,1,nn}]}],{t,0,20}];
Do[Subscript[fun3,i][t_]:=Evaluate[Subscript[QQ,i][t]/.sol4],{i,1,nn}];]


InterpolationTwoQuat[quatA_,quatB_,coupling_]:=Module[{K=coupling,q1,q2,q3,q4,t,eqns2},
k1=quatA;
k2=quatB;
P[t_]:={{q1[t]+I q2[t],-q3[t]+I q4[t]},{q3[t]+I q4[t],q1[t]-I q2[t]}};
eqns2={P'[t]==P[t].(-(K/2)ConjugateTranspose[ToMatrix[k2]+P[t]]).P[t]-ConjugateTranspose[-(K/2)ConjugateTranspose[ToMatrix[k2]+P[t]]],P[0]==ToMatrix[k1]};
sol2=NDSolve[eqns2,Flatten[{q1,q2,q3,q4}],{t,0,20}];
fun[t_]:=Evaluate[P[t]/.sol2];]


SimulationConsensus[time_]:=Module[{tt=time,g},
Do[Subscript[g,i]=Graphics3D[GeometricTransformation[ExampleData[{"Geometry3D","SpaceShuttle"},"GraphicsComplex"],RotationTransform[2PolarForm[Subscript[Quat,i][tt][[1]]][[1]],PolarForm[Subscript[Quat,i][tt][[1]]][[2]]]
],SphericalRegion->True],{i,1,n}];
Grid[ArrayReshape[Table[Subscript[g,i],{i,1,n}],{1,n}]]]


SimulationInterpolation[time_]:=Module[{tt=time,g1,g2},g1=Graphics3D[GeometricTransformation[ExampleData[{"Geometry3D","SpaceShuttle"},"GraphicsComplex"],RotationTransform[2PolarForm[ToQuaternion[fun[tt][[1]]]][[1]],PolarForm[ToQuaternion[fun[tt][[1]]]][[2]]]],SphericalRegion->True];g2=Graphics3D[GeometricTransformation[ExampleData[{"Geometry3D","Cow"},"GraphicsComplex"],RotationTransform[2PolarForm[k2][[1]],PolarForm[k2][[2]]]],SphericalRegion->True];
Grid[{{g1,g2}}]
];


SimulationDistribSync[time_]:=Module[{tt=time,gg},
Do[Subscript[gg,i]=Show[Graphics3D[If[tt<=tc,GeometricTransformation[ExampleData[{"Geometry3D","SpaceShuttle"},"GraphicsComplex"],RotationTransform[2PolarForm[ToQuaternion[Subscript[fun2,i][tt][[1]]]][[1]],PolarForm[ToQuaternion[Subscript[fun2,i][tt][[1]]]][[2]]]],SphericalRegion->True]],
Graphics3D[If[tt>=tc,GeometricTransformation[ExampleData[{"Geometry3D","SpaceShuttle"},"GraphicsComplex"],RotationTransform[2PolarForm[ToQuaternion[Subscript[fun3,i][tt-tc][[1]]]][[1]],PolarForm[ToQuaternion[Subscript[fun3,i][tt-tc][[1]]]][[2]]]]],SphericalRegion->True],SphericalRegion->True],{i,1,nn}];
Grid[ArrayReshape[Table[Subscript[gg,i],{i,1,nn}],{1,nn}]]]


VisualizationConsensus:=Manipulate[Quiet@SimulationConsensus[time],Item["Rotation average:",Alignment->Center],Item[Style[Framed[Subscript[Quat,1][100][[1]][[1]]+Subscript[Quat,1][100][[1]][[2]]"\[DoubleStruckI]"+Subscript[Quat,1][100][[1]][[3]]"\[DoubleStruckJ]"+Subscript[Quat,1][100][[1]][[4]]"\[DoubleStruckK]"],10,Blue],Alignment->Center],{{time,0,"time"},0,10,0.1,Appearance->"Labeled"},SaveDefinitions->True]


VisualizationInterpolation:=Manipulate[Quiet@SimulationInterpolation[time],Item["Quaternion A:",Alignment->Center],Item[Style[Framed[k1[[1]]+k1[[2]]"\[DoubleStruckI]"+k1[[3]]"\[DoubleStruckJ]"+k1[[4]]"\[DoubleStruckK]"],10,Blue],Alignment->Center],Item["Quaternion B:",Alignment->Center],Item[Style[Framed[k2[[1]]+k2[[2]]"\[DoubleStruckI]"+k2[[3]]"\[DoubleStruckJ]"+k2[[4]]"\[DoubleStruckK]"],10,Blue],Alignment->Center],{{time,0,"time"},0,20,0.1,Appearance->"Labeled"},SaveDefinitions->True]


VisualizationDistribSync:=Manipulate[Quiet@SimulationDistribSync[time],{{time,0.01,"time"},0,8,0.1,Appearance->"Labeled"},SaveDefinitions->True]


End[ ]


EndPackage[ ]
