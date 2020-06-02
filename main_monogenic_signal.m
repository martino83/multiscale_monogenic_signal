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


%------------------------------------------------------------------------
% This main file shows how to compute the monogenic signal on a single 
% test image, representing a 2D chirp function. The image contains a full 
% 360 degrees rotation and the frequency varies linearly.
% Image Source: http://bigwww.epfl.ch/demo/monogenic/index.html?ni=10&nf=11
%------------------------------------------------------------------------

clc
close all
clear 
addpath src

% given an input image computes the monogenic signal and extracts its
% features
im = imread('data/images/chirp2D.jpg');
im = im2graydouble(im);


% parameters settings
filter_family = 'DoP';
ratio = .98;
lambda = 10;
filter_parameters = [lambda, ratio];

orient_mode = 'robust'; 
freq_mode = 'robust';

% computation
parms = struct('filter_family',filter_family,...
    'filter_parameters',filter_parameters,...
    'orient_mode',orient_mode,...
    'freq_mode',freq_mode,'sigma',.5);
[fM feats] = compute_monogenic_signal(im,parms);

% display
figure,
imagesc(im), colormap gray, 
axis image off, title('2D Chirp','fontsize',12),

figure
theta = atan2(feats.sin_theta,feats.cos_theta);
imagesc(imrotate(-theta,90)), axis image, 
title('MONOGENIC ORIENTATION','fontsize',12),
set(gca,'Xtick',[])
set(gca,'Ytick',[])

figure
imagesc(feats.freq), axis image, title('MONOGENIC FREQUENCY','fontsize',12), 
caxis([0 2.214])
set(gca,'Xtick',[])
set(gca,'Ytick',[])





