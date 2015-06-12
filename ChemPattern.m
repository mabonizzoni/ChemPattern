(* ::Package:: *)

(* ::Section:: *)
(*A note on dataset normalization*)


(* ::Text:: *)
(*It is typically suggested that the dataset be normalized/rescaled before running LDA on it. However, normalization of all values to the interval (0,1) actually throws away information about the different dynamic ranges of different variables. If any normalization must be done, just in order to deal with numbers of the same order of magnitude, then Standardize is the most appropriate.*)
(**)
(*However, if the LDA is run on the normalized dataset, the eigenvectors cannot be properly used to transform new non-normalized data!! Also, new non-normalized data cannot be normalized in the same way that the original training set was, unless we record and transfer the maximum and minimum values used in the normalization of the original training set. This is probably a bad idea and a serious complication, at least at the initial level. *)
(**)
(*In the function below normalization is turned OFF by default, i.e. running LDA with no switch will assume that the data should not be normalized, and that the dataset has headers. If this assumption is incorrect, and in fact the dataset doesn't have headers, then the first row of the dataset will be discarded. Although this assumption may lead to wasting a sample replicate from the first group, it is still safer than assuming that no headers are present: in fact, in that case the headers will hang the calculation and one has to quit the kernel to recover.*)


(* ::Section:: *)
(*Functional implementation of LDA*)


(* ::Subsection:: *)
(*Things that need to be done:*)


(* ::Text:: *)
(*check whether to normalize full data set column-wise*)
(*generate a series of subsets by using the labels in the first column to gather data by its class*)
(*calculate within class scatter for each class (Swi) and sum them up to give total within-class scatter (Sw): \!\(\*UnderscriptBox[\(\[Sum]\), \(sets\)]\)\!\(\*UnderscriptBox[\(\[Sum]\), \(samples\ in\ a\ set\)]\)(x-Subscript[\!\(\*OverscriptBox[\(x\), \(_\)]\), set])\[CenterDot](x-Subscript[\!\(\*OverscriptBox[\(x\), \(_\)]\), set])^T*)
(*calculate between class scatter (Sb): \!\(\*UnderscriptBox[\(\[Sum]\), \(sets\)]\)(Subscript[\!\(\*OverscriptBox[\(x\), \(_\)]\), set]-Subscript[\[Mu], global])\[CenterDot](Subscript[\!\(\*OverscriptBox[\(x\), \(_\)]\), set]-Subscript[\[Mu], global])^T*)
(*solve eigensystem of matrix Sw^-1\[CenterDot]Sb*)
(*plot transformed dataset in 2D with tooltips*)
(*return eigensystem*)


(* ::Subsection:: *)
(*Implementation*)


Remove[lda];
Options[lda]={applyfunc->Identity,output->"2DL",ellipsoidcolor->Automatic,swapaxes->False};
lda::usage="lda[dataset] carries out Linear Discriminant Analysis on dataset and returns the transformed data as factor scores (default), or other numerical / graphical results. Each row of dataset should contain a sample; the first column contains the class identifier for that sample.\nOptions:\napplyfunc (Identity (default), Standardize, Rescale, ...)\noutput (\"scores\", \"vartable\", \"eigenvectors\", \"eigensystem\", \"2D\" (= 2D score plot), \"2DL\" (= 2D score and loading plots, default), \"3D\", \"3DL\")\nswapaxes (default = {False,False}).";
lda::outputoptions="The value `1` is not a valid plotting option. Valid options are: \"eigenvectors\" (default), \"eigensystem\", \"vartable\", \"2D\", \"3D\".";
lda::swapaxesnotboolean="The value `1` is not a valid swapaxes option. Use only True / False or combinations thereof.";
lda::swapaxeslength="The value `1` given for the swapaxes option is not valid. Acceptable values are a single True / False, to be applied to all axes, or a list containing as many True/False values as there are axes in the requested plot.";
lda::ellcolor="The value `1` for the ellipsoidcolor option is not valid. Acceptable values are: Automatic, True, False. Automatic settings were used.";


lda[matrix_/;MatrixQ[matrix],OptionsPattern[]]:=Module[
{(* the definition below is important, althouhg it looks silly!*)
(* If one used matrix directly in the function body, the VALUE of matrix would be substituted *)
(* everywhere BEFORE any further evaluation. the name of the pattern is not a proper variable!! *)
(* This is the reason why the standardization in place did not work before *)
dataset=matrix,

(* other local variables*)
hascolumnheaders,columnheaderlist,classlist,
partitioneddata,transformeddata,labeledtransformed,partitionedscores,
grandmean,clustermeans,Swi,Sw,Sb,
eigenvals,eigenvecs,
plotdata,readyforplot,
ellipsoids2D,coloredellipsoids2D,
ellipsoids3D,coloredellipsoids3D
},

(* if column headers are present, extract them and assign them to columnheaderlist. *)
(* if not, create a generic columnheaderlist and add it to the top of the dataset *)
(* this simplifies manipulation later on because we can assume the presence of column headers *)
(* and don't need to repeatedly check, no matter whether the original data had them or not *)
(* Remember that from now on dataset[[1,1]] contains no intereseting data *)
If[VectorQ[dataset[[1,2;;]],NumberQ]==True,
hascolumnheaders=False;Prepend[columnheaderlist=Join[{""},ToString/@Array["var",Last@Dimensions@dataset-1]]][dataset],
hascolumnheaders=True;columnheaderlist=dataset[[1,2;;]]
];

(* extract the class list from the first column in the dataset; remember that the first row is the column headers *)
classlist=DeleteDuplicates[dataset[[2;;,1]]];

(* Apply the applyfunc function to the numerical part of the dataaset *)
(* applyfunction's default is Identity, i.e. do nothing; if standardization is desired, *)
(* an appropriate standardizing function (e.g. Standardize itself) can be passed in *)
dataset[[2;;,2;;]]=OptionValue[applyfunc][dataset[[2;;,2;;]]];

(* Row headers are used as class labels to partition the data into the user-specified classes *)
(* The first row of dataset contains the column headers so it is ignored here *)
(* The class labels are then removed before assigning the list of partitioned datasets to partitioneddata *)
partitioneddata=GatherBy[dataset[[2;;]],First][[All,All,2;;]];

(* Calculate the within-cluster scatter by applying the Swi function to each cluster, and summing up each cluster's contribution *)
Swi[set_]:=Total[Map[Transpose[{#-Mean[set]}] . {#-Mean[set]}&,set]];
Sw=Total[Map[Swi,partitioneddata]];

(* Calculate the between-cluster scatter *)
(* In the case of a standardized dataset, the grand mean, i.e. the column-wise mean over the entire dataset should be exactly zero *)(* However it will typically calculate out to very small non-zero numbers at machine precision. To avoid accumulating inaccuracies, *)
(* the grand mean is chopped. If the dataset is not standardized, chopping has no significant effect *)
grandmean=Chop@Mean[dataset[[2;;,2;;]]];
clustermeans=Map[Mean,partitioneddata];
Sb=Total[Map[(Transpose[{#-grandmean}] . {#-grandmean}&),clustermeans]];

(* Maximize the Sb/Sw ratio by solving the equivalent eigenproblem *)
(* The Check[] wrapper around Inverse stops computation and returns if the inverse cannot be computed, e.g. for singular matrices *)
{eigenvals,eigenvecs}=Chop@Eigensystem[Check[Inverse[Sw] . Sb,Abort[],{Inverse::sing}]];

(* Calculate data scores along the LDA dimensions *)
transformeddata=Chop[dataset[[2;;,2;;]] . Transpose@eigenvecs];

(* Add back the column and row headers to the transformed data *)
labeledtransformed=Prepend[dataset[[1]]][Transpose@Insert[Transpose@transformeddata,dataset[[2;;,1]],1]];

Which[

(********************)
(* Numerical output *)
(********************)

OptionValue[output]=="scores",(* Default option: return the transformed data as labeled scores, e.g. for external plotting *)
Return[labeledtransformed],

OptionValue[output]=="vartable",(* Return a formatted table of the contributions of each variable to the first three factors *)
Return[
Style[
TableForm[Transpose@Round[100eigenvecs[[1;;3]]^2,1],TableHeadings->{columnheaderlist, {"F1","F2","F3"}},TableAlignments->Right],
FontFamily->"Arial",FontSize->14
]
],

OptionValue[output]=="eigenvectors",(* Return eigenvector matrix *)
Return[eigenvecs],

OptionValue[output]=="eigensystem", (* Return the list: {eigenvalues,eigenvectors} *)
Return[{eigenvals,eigenvecs}],

(***************)
(* 2D PLOTTING *)
(***************)

OptionValue[output]=="2D"||OptionValue[output]=="2DL",(* 2D plot of results was requested *)
plotdata=labeledtransformed[[All,1;;3]];

(* Swap values along x or y axes if requested; default is to do nothing *)
Module[
{xflip,yflip},
Switch[Length@OptionValue[swapaxes],

0, (* Atomic expression, i.e. a single value was passed *)
Switch[OptionValue[swapaxes],
True,xflip=yflip=-1,(* swap both axes *)
False,xflip=yflip=1,(* don't swap any axis *)
_,(* incorrect option; throw error and return Null *)Message[lda::swapaxesnotboolean,OptionValue[swapaxes]];Return[]
],

1,(* a list with one element; this is ambiguous and may be a syntax error on the part of the user; throw error *)
Message[lda::swapaxeslength,OptionValue[swapaxes]];Return[],

2,(* a list of two values*)
If[Not[BooleanQ@OptionValue[swapaxes][[1]]&&BooleanQ@OptionValue[swapaxes][[2]]],Message[lda::swapaxesnotboolean,OptionValue[swapaxes]];Return[]];
xflip=If[OptionValue[swapaxes][[1]]===True,-1,1];
yflip=If[OptionValue[swapaxes][[2]]===True,-1,1],

_,(* too many parameters for a 2D plot; possibly ambiguous *)
Message[lda::swapaxeslength,OptionValue[swapaxes]];Return[]
];

plotdata=plotdata/.List[class_?StringQ,x_?NumberQ,y_?NumberQ]->List[class,xflip x,yflip y]
];

(* Ellipsoids: in a non-correlated binormal distribution, 90% of the points lie within 2.15 standard deviations of the mean; 95% within 2.45 stdev; 99% within 3 stdev; 99.5% within 3.25 stdev *)
(* The ellipsoids are expressed as a function of the covariance; so n times StDev = n^2 times Covariance *)
(* To plot 95% confidence ellipsoids we need to stay within 2.45 stdev = (2.45)^2 covariance = ca. 6 covariance *)
partitionedscores=GatherBy[plotdata[[2;;]],First][[All,All,2;;3]];
ellipsoids2D={Opacity[0],EdgeForm[{Darker@Gray}],Ellipsoid[Mean[#],6Covariance[#]]}&/@partitionedscores;
coloredellipsoids2D=MapThread[
{Opacity[0],EdgeForm[{#2,AbsoluteThickness[2]}],Ellipsoid[Mean[#1],6Covariance[#1]]}&,
{
partitionedscores,
ColorData[97,"ColorList"][[1;;First@Dimensions@partitionedscores]]
}];

readyforplot=MapThread[Tooltip,{partitionedscores,classlist}];

Return[
If[OptionValue[output]=="2DL",GraphicsRow[#,ImageSize->Scaled[0.6]]&,Show[#[[1]],ImageSize->Scaled[0.3]]&]@
List[
(* 2D score plot *)
ListPlot[readyforplot,
Frame->True,FrameStyle->Directive[Black,FontSize->15],FrameLabel->{
Style["Factor 1 ("<>ToString[Round[100eigenvals[[1]]/Total@eigenvals,0.1]]<>"%)",FontSize->16,Blue],
Style["Factor 2 ("<>ToString[Round[100eigenvals[[2]]/Total@eigenvals,0.1]]<>"%)",FontSize->16,Red]
},
AspectRatio->1,PlotRange->All,PlotRangePadding->Scaled[0.10],
Epilog->Which[
OptionValue[ellipsoidcolor]===Automatic||OptionValue[ellipsoidcolor]===False,ellipsoids2D,
OptionValue[ellipsoidcolor]===True,coloredellipsoids2D,
True,Message[lda::ellcolor,OptionValue[ellipsoidcolor]];ellipsoids2D
]
],
If[OptionValue[output]=="2DL",
(* 2D loading plot *)
ListPlot[
MapThread[Labeled[100#1,Style[#2<>" "<>ToString@Round[100#1,1],Medium]]&,{Transpose[eigenvecs[[1;;2]]^2],columnheaderlist}],
PlotStyle->Directive[Black,PointSize[0.025]],
(* The aspect ratio and plotrange definitions below make the plot square, while still adapting the plot range to the values being plotted *)
AspectRatio->1,PlotRange->{{0,105Max[Transpose[eigenvecs[[1;;2]]^2]]},{0,105Max[Transpose[eigenvecs[[1;;2]]^2]]}},PlotRangePadding->Scaled[.05],
AxesOrigin->{0,0},
Frame->{True,True,False,False},FrameStyle->Directive[Black,FontSize->15],FrameLabel->{
Style["Contrib. to F1 (%)",FontSize->16,Blue],
Style["Contrib. to F2 (%)",FontSize->16,Red]
}
],
(* no 2D loading plot: add "nothing" *)
Unevaluated@Sequence[]
]
]
],

(***************)
(* 3D PLOTTING *)
(***************)

OptionValue[output]=="3D"||OptionValue[output]=="3DL",(* 3D plot of results was requested *)
plotdata=labeledtransformed[[All,1;;4]];

(* Swap values along x or y axes if requested; default is to do nothing *)
Module[
{xflip,yflip,zflip},
Switch[Length@OptionValue[swapaxes],

0, (* Atomic expression, i.e. a single value was passed *)
Switch[OptionValue[swapaxes],
True,xflip=yflip=zflip=-1,(* swap both axes *)
False,xflip=yflip=zflip=1,(* don't swap any axis *)
_,(* incorrect option; throw error and return Null *)Message[lda::swapaxesnotboolean,OptionValue[swapaxes]];Return[]
],

1,(* a list with one element; this is ambiguous and may be a syntax error on the part of the user; throw error *)
Message[lda::swapaxeslength,OptionValue[swapaxes]];Return[],

2,(* a list with two elements; not enough for a 3D plot; throw error *)
Message[lda::swapaxeslength,OptionValue[swapaxes]];Return[],

3,(* a list of three values*)
If[Not[BooleanQ@OptionValue[swapaxes][[1]]&&BooleanQ@OptionValue[swapaxes][[2]]&&BooleanQ@OptionValue[swapaxes][[3]]],Message[lda::swapaxesnotboolean,OptionValue[swapaxes]];Return[]];
xflip=If[OptionValue[swapaxes][[1]]===True,-1,1];
yflip=If[OptionValue[swapaxes][[2]]===True,-1,1];
zflip=If[OptionValue[swapaxes][[3]]===True,-1,1],

_,(* too many parameters for a 3D plot; possibly ambiguous *)
Message[lda::swapaxeslength,OptionValue[swapaxes]];Return[]
];

plotdata=plotdata/.List[class_?StringQ,x_?NumberQ,y_?NumberQ,z_?NumberQ]->List[class,xflip x,yflip y,zflip z]
];

(* Ellipsoids: see above in plot 2D for discussion of width of ellipsoids *)
partitionedscores=GatherBy[plotdata[[2;;]],First][[All,All,2;;4]];
ellipsoids3D={Opacity[0.1,Black],Ellipsoid[Mean[#],6Covariance[#]]}&/@partitionedscores;
coloredellipsoids3D=MapThread[
{Opacity[0.2,#2],Ellipsoid[Mean[#1],6Covariance[#1]]}&,
{
partitionedscores,
ColorData[97,"ColorList"][[1;;First@Dimensions@partitionedscores]]
}];

readyforplot=MapThread[Tooltip[{#3,Point[#1]},#2]&,{partitionedscores,classlist,ColorData[97,"ColorList"][[1;;First@Dimensions@partitionedscores]]}];

Return[
If[OptionValue[output]=="3DL",GraphicsRow[#,ImageSize->Scaled[1]]&,Show[#[[1]],ImageSize->Scaled[0.6]]&]@
List[
(* 3D score plot *)
Graphics3D[{
{PointSize->0.006,readyforplot},
Which[
OptionValue[ellipsoidcolor]===Automatic||OptionValue[ellipsoidcolor]===True,coloredellipsoids3D,
OptionValue[ellipsoidcolor]===False,ellipsoids3D,
True,Message[lda::ellcolor,OptionValue[ellipsoidcolor]];coloredellipsoids3D
]
},
Axes->True,AxesLabel->{
Style["Factor 1 ("<>ToString[Round[100eigenvals[[1]]/Total@eigenvals,0.1]]<>"%)",FontSize->Scaled[0.04],FontFamily->"Arial",Blue],Style["Factor 2 ("<>ToString[Round[100eigenvals[[2]]/Total@eigenvals,0.1]]<>"%)",FontSize->Scaled[0.04],FontFamily->"Arial",Red],
Style["Factor 3 ("<>ToString[Round[100eigenvals[[3]]/Total@eigenvals,0.1]]<>"%)",FontSize->Scaled[0.04],FontFamily->"Arial",Darker@Green]
},
PlotRange->All,PlotRangePadding->Scaled[0.05],BoxRatios->{1, 1, 1},Lighting->"Neutral",RotationAction->"Clip"
],
If[OptionValue[output]=="3DL",
(* 3D loading plot *)
Graphics3D[
MapThread[Tooltip[Style[Point[#1],Red,PointSize[0.02]],Style[#2<>" "<>ToString@Round[100#1,1],Medium]]&,{Transpose[eigenvecs[[1;;3]]^2],columnheaderlist}],
PlotRange->{{0,1},{0,1},{0,1}},PlotRangePadding->Scaled[.05],
Axes->True,AxesLabel->{
Style["Contrib. to F1",FontSize->Scaled[0.04],FontFamily->"Arial",Blue],
Style["Contrib. to F2",FontSize->Scaled[0.04],FontFamily->"Arial",Red],
Style["Contrib. to F3",FontSize->Scaled[0.04],FontFamily->"Arial",Darker@Green]
},
RotationAction->"Clip"
],
(* no 3D loading plot: add "nothing" *)
Unevaluated@Sequence[]
]
]
],

(************************************)
(* Default: incorrect output option *)
(************************************)

True,(* default option: if this is reached, the selected output option is incorrect, so throw error and return Null *)
Message[lda::outputoptions,OptionValue[output]];Return[]
];
]


(* This alternative definition of lda handles malformed input *)
lda::notamatrix="The input to lda is not a matrix.";
lda[dataset_/;Not[MatrixQ[dataset]],OptionsPattern[]]:=Message[lda::notamatrix]



(* The following function generates a bar chart of the weighted contributions of each sensor to the overall discrimination *)
(* The contributions are weighted by the weight of the factors themselves, represented by the eigenvalues from lda *)
(* This is because a sensor that contributes a lot to an unimportant factor is still unimportant in the overall discrimination *)
Remove[groupcontribs];

groupcontribs::usage="groupcontribs[eigensystem,numberofgroups,sensornames]\ngroupedcontributions[eigensystem,numberofgroups]\nThe function generates a bar chart of the contributions of each group of variables to the overall discrimination. Before summing, the contributions to each factor are weighted by the corresponding eigenvalue of the factor. This is needed so that a group that contributes a lot to an unimportant factor is still reported as unimportant in the overall discrimination.";
groupcontribs::numgroups="The number of variables in the eigensystem (`1`) is not an exact multiple of the number of groups provided (`2`).";

groupcontribs[eigensystem_,numberofgroups_,sensornames_:Null]:=Module[
{eigenvals=eigensystem[[1]],eigenvecs=eigensystem[[2]],variablespergroup,sqweightedeigenvecs,barvalues},

If[Mod[Length@eigenvecs,numberofgroups]==0,
variablespergroup=(Length@eigenvecs)/numberofgroups,
Message[groupcontribs::numgroups,Length@eigenvecs,numberofgroups];Abort[]
];

(* eigenvectors are weighted by the corresponding eigenvalues, then squared *)
sqweightedeigenvecs=(Normalize[eigenvals,Total]eigenvecs)^2;

barvalues=Round[100#,1]&@
Normalize[#,Total]&@
Table[
Chop@Total[sqweightedeigenvecs[[All,i;;variablespergroup-1+i]]^2,Infinity],
{i,1,Last@Dimensions@eigenvecs-variablespergroup+1,variablespergroup}];

BarChart[
barvalues,
ChartLabels->Map[Style[#,FontSize->20,Black]&,If[sensornames===Null,Array[group,(Last@Dimensions@eigenvecs)/variablespergroup],sensornames]],
BarSpacing->Large,
PlotLabel->Style["% contribution of each sensor (overall)",Black,FontFamily->"Arial",FontSize->20],
ImageSize->Scaled[0.25]
]
]



(* ::Section:: *)
(*Detection of outliers using PCA*)


outlierPCA[set_]:=Module[
{workingdata,eigenvectors,eigenvalues,PCs,contributions},
workingdata=If[
NumberQ[set[[ 1,2]]],
set[[All,2;;]],(* the first row contains data, and not variable labels: do not discard it! *)
set[[2;;,2;;]](* the first row contains variable labels: discard it *)
];
{eigenvalues,eigenvectors}=Eigensystem@Correlation@workingdata;
PCs=Standardize[workingdata] . Transpose[eigenvectors];
contributions=Round[100Normalize[eigenvalues,Total],0.1];
Show[{
ListPlot[
MapIndexed[Labeled[#1,First@#2]&,PCs[[All,1;;2]]],
AspectRatio->1,PlotStyle->PointSize[0.02],
Axes->False,Frame->True,FrameStyle->Directive[Black,FontSize->14],
FrameLabel->{"PC1 ("<>ToString[contributions[[1]]]<>"%)","PC2 ("<>ToString[contributions[[2]]]<>"%)"},
Epilog->Inset[Style[set[[2,1]],Red,Bold,FontSize->18],Scaled[{0.9,0.9}]]
],
Graphics@{Opacity[0],EdgeForm[{Gray,Dashed,Thick}],Ellipsoid[Mean@PCs[[All,1;;2]],6Covariance@PCs[[All,1;;2]]]}
},
PlotRange->All,PlotRangePadding->Scaled[0.1]
]
]


(* ::Section:: *)
(*Saving out function definitions*)


Save["C:\\Users\\Marco\\Documents\\Alabama\\Dissemination\\Papers\\2015 Alie Wallace coumarins\\Data\\currentLDAfunctions.m",{lda,groupcontribs,outlierPCA}]
