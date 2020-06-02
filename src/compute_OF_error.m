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
% M. Alessandrini, A. Basarab, H. Liebgott and O. Bernard, "Multiscale 
% Optical Flow Computation from the Monogenic Signal", submitted fot
% buplication to IEEE Transactions on Image Processing
%------------------------------------------------------------------------

function [AE EE ME] = compute_OF_error(ue,ve,uc,vc,mask)

me = sqrt(ue.^2 + ve.^2);
mc = sqrt(uc.^2 + vc.^2);

mask = mask & abs(me)>0 & abs(mc)>0;
idx_eval = find(mask==1);

%%% angluar error of Barron
ae = (uc.*ue + vc.*ve + 1)./sqrt((uc.^2+vc.^2+1).*(ue.^2+ve.^2+1));
ae = acos(ae);
ae = ae*180/pi;

AE.image = ae; % image for display
AE.image(mask == 0) = 0; % image for display

% ae = ae(mask==1); 
% ae = ae(:);% for statistics computation
AE.mean = mean(ae(idx_eval));
AE.std = std(ae(idx_eval));

%%% magnitude error
mag_e = sqrt((ue).^2 + (ve).^2);
mag_c = sqrt((uc).^2 + (vc).^2);
me = abs(mag_e - mag_c);

ME.image = me;
ME.image(mask == 0) = 0; % image for display

% me = me(mask == 1);
ME.mean = mean(me(idx_eval));
ME.std = std(me(idx_eval));

%%% endpoint error
ee = sqrt((ue-uc).^2 + (ve-vc).^2);
EE.image = ee;

% ee = ee(mask==1); ee = ee(:);
EE.mean = mean(ee(idx_eval));
EE.std = std(ee(idx_eval));

%%% robusteness
AE.R2p5 = robustness(ae,2.5);
AE.R5 = robustness(ae,5);
AE.R10 = robustness(ae,10);

EE.R0p5 = robustness(ee,0.5);
EE.R1 = robustness(ee,1);
EE.R2 = robustness(ee,2);

%%% accuracy
AE_Ax = accuracy(ae,[50 75 95]);
EE_Ax = accuracy(ee,[50 75 95]);
AE.A50 = AE_Ax(1);
AE.A75 = AE_Ax(2);
AE.A95 = AE_Ax(3);

EE.A50 = EE_Ax(1);
EE.A75 = EE_Ax(2);
EE.A95 = EE_Ax(3);

function RX = robustness(err,X)
N = numel(err);
NX = sum(err > X);
RX = NX/N;

function AX = accuracy(err,X)
AX = prctile_int(err,X); 
