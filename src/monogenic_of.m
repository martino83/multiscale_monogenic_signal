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

%------------------------------------------------------------------------
% Multiscale Monogenic Optical Flow Algorithm.
%
% [v elapsed_time] = monogenic_of(I,options)
%
% Input:
% I        - subsequent frames in the video sequence built likewise:
%            I = cat(3,im1,im2) if forward differencing is employed
%            I = cat(3,im1,im2,im3) if center differences
% options  - structure containing the parameters for the optical flow
%            computation. Field values:
%            .monogenic_parms: contains the parameters for the monogenic
%             signal computation. type 'help compute_monogenic_signal' for
%             further details
%            .binomial_sigma: variance of the binomial filter for image
%             pre_filtering (default is 1);
%            .J_coarsest: coarsest scale for moments computation
%             (default is 5)
%            .J_finest: finest scale for moments computation
%             (default is 2)
%            .bspline: bspline used, specified as a string 
%             (ex 'b3' or 'b5'). odd orders are expected.
%            .system_cond: the solution is accepted if the system matrix is
%             well conditioned. The control can be done on the condition
%             number or on the minimum eigenvalue. This field is specified
%             as a cell {check_tyoe,value}. Where check_type can be
%             'min_eig' or 'cond_number'.
%            .mode: can be 'lucas_kanade' (translation only) or
%             'spatial_affine' (default).
%            .mask: binary mask indicating the points where the optical 
%             flow should be computed. Default is all the image
%
% Output
% v_of     - estimated displacement:
%            u = v_of(:,:,1); horizontal component
%            v = v_of(:,:,1); vertical component
% elapsed_time - time elapsed for the optical flow computation
%
%------------------------------------------------------------------------

function [v elapsed_time] = monogenic_of(I,options)


if options.J_coarsest < options.J_finest
    error('condition J_coarsest >= J_finest violated');
end
        
% sequence pre-smoothing
if options.binomial_sigma > 0
    I = binomial_filter(I,options.binomial_sigma);    
end

% adjust mask size
if numel(options.mask) == 1 && options.mask == 1
   options.mask = ones([size(I,1) size(I,2)]); 
end

% build bspline pyramid, and filters for moments (I and II order) computaion
[b wfilt w2filt] = ...
    build_b_pyramid(options.J_coarsest,options.bspline);
N = str2double(options.bspline(2)); % bspline order

% filters for two scale filters computation
H = two_scale_filters_computation(options.J_coarsest,N,options.mode);

tic;
% tic;
[M_Ix2 M_Iy2 M_Ixy M_Ixt M_Iyt M_It] = ...
    compute_multiresolution_moments(I,H,options,b,wfilt,w2filt);

% compute velocity
Nscales = options.J_coarsest - options.J_finest + 1; 
scale_id = options.J_coarsest+1;

if numel(options.dt_threshold) == 4
    [x1 y1 x2 y2] = ...
        deal(options.dt_threshold(1),options.dt_threshold(2),options.dt_threshold(3),options.dt_threshold(4));
    [X Y] = meshgrid(1:size(I,2),1:size(I,1));
    idxNoiseRegion = X > x1 & X < x2 & Y > y1 & Y < y2;
    maskNoiseRegion = zeros(size(I,1),size(I,2));
    maskNoiseRegion(idxNoiseRegion) = 1;
else
    maskNoiseRegion = 0;
end



%%%% velocity computation at the coarsest scale
[vJcoarsest EJcoarsest] = velocity_computation_coarsest_scale...
    (M_Ix2{scale_id},M_Ixy{scale_id},M_Iy2{scale_id},M_Ixt{scale_id},...
    M_Iyt{scale_id},M_It{scale_id},options,b{scale_id},...
    options.mask(1:2^(scale_id-1):end,1:2^(scale_id-1):end),...
    maskNoiseRegion(1:2^(scale_id-1):end,1:2^(scale_id-1):end));

E{Nscales}.val = EJcoarsest;
E{Nscales}.scale = options.J_coarsest;
DD = 1;
if (options.J_coarsest == 0) || (options.J_coarsest == options.J_finest)
    [P] = size(vJcoarsest,3);
    [idx] = ( sqrt(sum(vJcoarsest(:,:,1:2).^2,3)) > options.max_magnitude);
    idx = repmat(idx,[1 1 P]);
    vJcoarsest(idx) = 0;
    v{Nscales}.val = vJcoarsest;
    v{Nscales}.scale = options.J_coarsest;
    
else
    v{Nscales}.val = vJcoarsest;
    v{Nscales}.scale = options.J_coarsest;

    % posteriori test
    vjj = vJcoarsest;
    Ejj = EJcoarsest;
    clear vJcoarsest EJcoarsest
    
    kk = 0;
    for scale_id = options.J_coarsest:-1:(options.J_finest+1)
        
       % transfer coarsest estimate to finer scale
       switch options.mode
           case {'spatial_affine','lucas_kanade'}
               final_size = size(M_Ix2{scale_id}.m_00);
           case 'spatio_temporal_affine'
               final_size = size(M_Ix2{scale_id}.m_000);
       end
       [v_transf E_transf] = interp_velocity_field(vjj,Ejj,final_size);
       
       % recompute velocity field if necessar
       [vjj Ejj] = velocity_computation_gen_scale...
            (M_Ix2{scale_id},M_Ixy{scale_id},M_Iy2{scale_id},...
            M_Ixt{scale_id},M_Iyt{scale_id},M_It{scale_id},...
            options,b{scale_id},...
            v_transf,E_transf,...
            options.mask(1:2^(scale_id-1):end,1:2^(scale_id-1):end),...
            maskNoiseRegion(1:2^(scale_id-1):end,1:2^(scale_id-1):end));
       debug()
       kk = kk+1;
       E{Nscales-kk}.val = Ejj;
       E{Nscales-kk}.scale = scale_id-1;
       v{Nscales-kk}.val = vjj;
       v{Nscales-kk}.scale = scale_id-1;
       
    end
    
    [P] = size(vjj,3);
    [idx] = ( sqrt(sum(vjj(:,:,1:2).^2,3)) > options.max_magnitude);
    idx = repmat(idx,[1 1 P]);
    vjj(idx) = 0;
    v{1}.val = vjj;
    v{1}.scale = options.J_finest;
end
elapsed_time = toc;
            


    
