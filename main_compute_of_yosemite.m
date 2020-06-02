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
% This main file runs the multiscale monogenic oprical flow algorithm on
% frame 9 of the yosemite sequence.
% Yosemite sequence was obtained at Prof. Michael Black home page:
% http://www.cs.brown.edu/~black/images.html
%------------------------------------------------------------------------
clc
close all
clear
addpath src


%-- 1/ load data
data_folder = 'data/sequences/yosemite';

% load benchmark
v_bench_name = [data_folder '/yos_img_09_bench.mat']; 

% load cloud mask (yosemite without clouds)
mask_cloud_name = [data_folder '/yos_img_09_mask_cloud.mat'];

% load three frames (center difference method)
im0 = imread([data_folder '/yos_img_08.tif']);
im1 = imread([data_folder '/yos_img_09.tif']);
im2 = imread([data_folder '/yos_img_10.tif']);

im0 = im2graydouble(im0);
im1 = im2graydouble(im1);
im2 = im2graydouble(im2);

I = cat(3,im0,im1,im2);
I = I - min(I(:));
I = 255 * I./max(I(:));

%-- 2/ set parameters for the computation of the monogenic signal:
% - SQF specification () 
filter_family = 'DoP';
ratio = .98;  % ratio parameter
lambda_0 = 3; % wavelength (1/f_0)
filter_parameters = [lambda_0, ratio];

% - monogenic orientation computation mode ('robust' or 'classical')
orient_mode = 'robust';
freq_mode = 'robust'; 

% - setup parameters structure
monogenic_parms = struct('filter_family',filter_family,...
    'filter_parameters',filter_parameters,...
    'orient_mode',orient_mode,...
    'freq_mode',freq_mode,'sigma',2);

%-- 3/ motion estimation parameters
% type 'help bspline_optical_flow_multires' for more info
options.monogenic_parms = monogenic_parms;
options.J_coarsest = 5;
options.J_finest = 3;
options.binomial_sigma = 1;
options.bspline = 'b5';
options.dt_threshold = .1; 
options.max_magnitude = 5;
options.padding = 0;
options.mode = 'spatial_affine';
nscales = 4;

% compute of
[Ve elapsedTime] = pyramidal_refinement_monogenic_of(I,'options',options,...
    'nscales',nscales,'mult',1.5);
ue = Ve(:,:,1);
ve = Ve(:,:,2);

% load benchmark and mask (for errors computation)
load(v_bench_name);

load(mask_cloud_name);
mask = mask_cloud;
mask(:,[1:21 end-20:end]) = 0; 
% to remove boundary effects a frame of 20 pixels is considered in the 
% horizontal direction. Not considering boundary artifacts is coherent with 
% the middelbury database evaluation criterion. See: 
% S. Baker, D. Scharstein, J. Lewis, S. Roth, M. Black,
% and R. Szeliski. "A database and evaluation methodology for
% optical flow." International Conference of Computer Vision, (ICCV 2007)


[AE EE ME] = compute_OF_error(ue,ve,u,v,mask);

figure,
im_ref = I(:,:,1);
[xI yI] = meshgrid(1:size(im_ref,2),1:size(im_ref,1));
dd = 5;
imagesc(im_ref), colormap gray, axis image, hold on,
xd = xI(1:dd:end,1:dd:end);
yd = yI(1:dd:end,1:dd:end);

ud = ue(1:dd:end,1:dd:end);
vd = ve(1:dd:end,1:dd:end);
quiver(xd,yd,ud,vd,0,'b');
title('RED = benchmark, BLUE = computed field')

% benchmark
u_bench = u(1:dd:end,1:dd:end);
v_bench = v(1:dd:end,1:dd:end);
quiver(xd,yd,u_bench,v_bench,0,'r');

disp('******************** Results ********************')                    
disp(['Angular Error = ' num2str(AE.mean) ' pm ' num2str(AE.std)]);
disp(['Magnitude Error = ' num2str(ME.mean) ' pm ' num2str(ME.std)]);
disp(['EndPoint Error = ' num2str(ME.mean) ' pm ' num2str(ME.std)]);
disp('*************************************************')
