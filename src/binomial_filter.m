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

function [If bfilt] = binomial_filter(I,sigma)

if sigma == 0
    If = I;
    bfilt = [];
    return
end

if sigma < 0
    error('filter variance must be positive');
end

% filter with binomial filter of variance sigma
N = 4*sigma^2; % filter order
if floor(N) ~= N
    error('binomial filter error: 4*sigma^2 must be an integer');
end

% bild binomial filter
bfilt = zeros(N+1,1);
for kk = 0:N
    bfilt(kk+1) =  nchoosek(N,kk);
end
bfilt = bfilt/sqrt(sum(bfilt(:).^2));
If = I; % separable convolution
for kk = 1:size(I,3)
    If(:,:,kk) = conv2(bfilt,bfilt,If(:,:,kk),'same');
end