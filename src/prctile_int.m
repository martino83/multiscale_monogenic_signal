%------------------------------------------------------------------------
% Copyright or Â© or Copr. CREATIS laboratory, Lyon, France.
% 
% Contributor: Martino Alessandrini, Post Doctoral Fellow at the 
% Centre de Recherche en Acquisition et Traitement de l'Image pour la Santé
% CREATIS (CNRS 5220, INSERM U630, INSA, Claude Bernard Lyon 1 University) 
% in France (Lyon).
% 
% Date of creation: March 26th 2012
% 
% E-mail of the author: martino.alessandrini@creatis.insa-lyon.fr
% 
% This folder provides a MATLAB implementation of an Optical Flow
% estimation algorithm based on the monogenic phase. Given two input images
% the algorithm compute the displacement field between the two by assuming
% the conservation of the monogenic phase. This feature is much less
% sensitive to changes in the illumination conditions as compared to the
% traditional pixel intensity. To reduce dependency on the size of the
% windowing function, the computation is carried out at different scales in
% a coarse-to-fine fashion. The estimation is then refined iteratively in a
% pyramidal scheme.
% 
% The algorithm herein implemented is described in:
% M. Alessandrini, A. Basarab, H. Liebgott and O. Bernard, "Myocardial 
% Motion Estimation from Medical Images Using the Monogenic Signal, 
% accepted for buplication to IEEE Transactions on Image Processing
%------------------------------------------------------------------------

function y = prctile_int(x,p)
%PRCTILE gives the percentiles of the sample in X.
%	Y = PRCTILE(X,P) returns a value that is greater than P percent
%	of the values in X. For example, if P = 50  Y is the median of X. 
%
%	P may be either a scalar or a vector. For scalar P, Y is a row	
%	vector containing Pth percentile of each column of X. For vector P,
%	the ith row of Y is the P(i) percentile of each column of X.

%	Copyright (c) 1993 by The MathWorks, Inc.
%	$Revision: 1.1 $  $Date: 1993/05/24 18:56:10 $

[prows pcols] = size(p);
if prows ~= 1 & pcols ~= 1
    error('P must be a scalar or a vector.');
end

if any(p > 100) | any(p < 0)
    error('P must take values between 0 and 100');
end

xx = sort(x);
[m,n] = size(x);

if m==1 | n==1
    m = max(m,n);
    n = 1;
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx(:); max(x)];
else
    q = 100*(0.5:m - 0.5)./m;
    xx = [min(x); xx; max(x)];
end


q = [0 q 100];

y = interp1(q,xx,p);
