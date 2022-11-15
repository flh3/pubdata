/* 
When sharing this program and data sets please include the citation:
  Stefanski, Leonard A. (2007). Residual (Sur)Realism. "The American 
      Statistician," 61, pp 163-177.
Online version location ...
http://www4.stat.ncsu.edu/~stefanski/NSF_Supported/Hidden_Images/Residual_Surrealism_TAS_2007.pdf
*/ run quit;

%let dataname=mizzo_1_data_yx1x5;
%let dx=5;
filename mypath url 
'http://www4.stat.ncsu.edu/~stefanski/NSF_Supported/Hidden_Images/MIZZOU_Tiger_Files/';

data &dataname;
infile mypath(&dataname..txt);
input y x1-x&dx.;
%let xdata=x1-x&dx.;
run;quit;

proc reg;
title1 "Correct Model";
model y = &xdata / collin collinoint;
output out=regoutcm r=rcm p=pcm;
run; quit;

proc gplot;
title2 "Residual vs Fitted Values";
plot rcm*pcm;
symbol v=dot height=.2;
run; quit;
