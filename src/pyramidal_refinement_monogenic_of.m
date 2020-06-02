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
% [v_of elapsed_time] = pyramidal_refinement_monogenic_of...
%     (im1,im2,options,varargin)
%
% Input:
% I        - subsequent frames in the video sequence built likewise:
%            I = cat(3,im1,im2) if forward differencing is employed
%            I = cat(3,im1,im2,im3) if center differences
% varargin - variable input arguement are specified by the couple:
%            (parameter_name,value), where parameter_name can be:
%            'options': structure containing the parameters for the optical flow
%            computation at each refinement iteration.
%            'nscales': specifies the number of refinement iterations
%                       (defalt is 3)
%            'mult'   : multiplicative factor between successive sqf 
%                       wavelengths (default is 1.5)
%
% Output
% v_of     - estimated displacement:
%            u = v_of(:,:,1); horizontal component
%            v = v_of(:,:,1); vertical component
% elapsed_time - time elapsed for the optical flow computation
%
%
% example of usage
% v = pyramidal_refinement_monogenic_of(I,'options',options,'nscales',4,'mult',1.5)
%------------------------------------------------------------------------

function [v_of elapsed_time] = pyramidal_refinement_monogenic_of...
    (I,varargin)


% fix variable arguements
nscales = 3;
mult = 1.5;
options = 0;
for ii = 1:2:numel(varargin)
    switch varargin{ii}
        case 'nscales'
            nscales = varargin{ii+1};
        case 'mult'
            mult = varargin{ii+1};
        case 'options'
            options = varargin{ii+1};
        otherwise
            error('undefined parameter');
    end
end

% fix options structure
options = fix_options(options);


nframes = size(I,3);
wavelength = options.monogenic_parms.filter_parameters(1);
imJ = I;

% initialization
[xI yI] = meshgrid(1:size(imJ,2),1:size(imJ,1));
uJ = zeros(size(imJ,1),size(imJ,2)); % 
vJ = zeros(size(uJ));

freq_vec = 1/wavelength*mult.^(0:nscales-1);
wave_vec = 1./freq_vec;

optionsJ = options;
interp_type = 'linear';
elapsed_time = 0;

for J = 1:nscales % main cycle
    
    % wavelength used in teh current iteration
    wavelengthJ = wave_vec(J);
    optionsJ.monogenic_parms.filter_parameters(1) = wavelengthJ;
    disp(['---- computing scale ' num2str(J) ', wavelength = ' num2str(wavelengthJ)])
    
    % compute incrementsl of
    tic
    [dV] = monogenic_of(imJ,optionsJ);
    dd = 2^dV{1}.scale;
    
    % get dense field (interpolation)
    if dd > 1
        x = xI(1:dd:end,1:dd:end);
        y = yI(1:dd:end,1:dd:end);
  
        duJ = interp2(x,y,dV{1}.val(:,:,1),xI,yI,interp_type,0);
        dvJ = interp2(x,y,dV{1}.val(:,:,2),xI,yI,interp_type,0);        
    else         
        duJ = dV{1}.val(:,:,1);
        dvJ = dV{1}.val(:,:,2);
    end
    
    % final displacement
    uJ = uJ + duJ;
    vJ = vJ + dvJ;
    elapsed_time = elapsed_time + toc;
        
    % subtract displacement (warp)
    if J < nscales            
        switch nframes
            case 2
                imJ(:,:,2) = interp2(xI,yI,I(:,:,2),xI+uJ,yI+vJ,interp_type,0);
            case 3 
                imJ(:,:,3) = interp2(xI,yI,I(:,:,3),xI+uJ,yI+vJ,interp_type,0);
                imJ(:,:,1) = interp2(xI,yI,I(:,:,1),xI-uJ,yI-vJ,interp_type,0);
        end     
    end    
end

elapsed_time = elapsed_time./nscales; % time per iteration
v_of = cat(3,uJ,vJ); % solution

function optionsOut = fix_options(optionsIn)

if ~isstruct(optionsIn), 
    optionsOut.binomial_sigma = 1; 
    optionsOut.J_coarsest = 5;
    optionsOut.J_finest = 2; 
    optionsOut.bspline = 'b5'; 
    optionsOut.dt_threshold = 0.1; 
    optionsOut.system_cond = {'min_eig',18}; 
    optionsOut.max_magnitude = 2; 
    optionsOut.mode = 'spatial_affine'; 
    optionsOut.mask = 1; 
    optionsOut.padding = 0; 
    optionsOut.monogenic_parms.filter_parameters = [8 .98];
    optionsOut.monogenic_parms.filter_family = 'DoP';
else
    optionsOut = optionsIn;
    if ~isfield(optionsOut,'binomial_sigma'), optionsOut.binomial_sigma = 1; end
    if ~isfield(optionsOut,'J_coarsest'), optionsOut.J_coarsest = 5; end
    if ~isfield(optionsOut,'J_finest'), optionsOut.J_finest = 2; end
    if ~isfield(optionsOut,'bspline'), optionsOut.bspline = 'b5'; end
    if ~isfield(optionsOut,'dt_threshold'), optionsOut.dt_threshold = 0.1; end
    if ~isfield(optionsOut,'system_cond'), optionsOut.system_cond = {'min_eig',20}; end
    if ~isfield(optionsOut,'max_magnitude'), optionsOut.max_magnitude = 2; end
    if ~isfield(optionsOut,'mode'), optionsOut.mode = 'spatial_affine'; end
    if ~isfield(optionsOut,'mask'), optionsOut.mask = 1; end
    if ~isfield(optionsOut,'padding'), optionsOut.padding = 0; end
    if ~isfield(optionsOut,'filter_parameters'), optionsOut.filter_parameters = [8 .98]; end
    if ~isfield(optionsOut,'monogenic_parms'), 
        optionsOut.monogenic_parms.filter_parameters = [8 .98];
        optionsOut.monogenic_parms.filter_family = 'DoP';
    end
end

