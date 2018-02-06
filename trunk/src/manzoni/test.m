(* ::Package:: *)

Save["data.inp",data];
plot=ListPlot[data[[All,3;;4]]];
Export["test.png",plot];



