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

function [v E idxSufficientDt] = velocity_computation_coarsest_scale...
    (M_Ix2,M_Ixy,M_Iy2,M_Ixt,M_Iyt,M_It,...
    options,wJ,mask,maskNoise)

% Compute velocity estimate at coarsest scale given image moments. 
% Input:
% M_Ix2,M_Ixy,M_Iy2,M_Ixt,M_Iyt,M_It2 - image moments
% mode - which velocity model is adopted:
%                       'spatio_temporal_affine' as for myocardium 
%                       'spatial_affine' as for natural images.
% eigthresh - threshold on the minimum eigenvalue of A
% idxEval - pixels where the fied has to be computed
% wJ - window used for moments computation
%
% Output
% v - velocity solution. 3D matrix containing velocity estimate for each
% pixel
% F - 2D matrix containing a flag value for each pixel :
%                   0 if matrix is well conditioned
%                   1 if ill conditioned (beta)
%                   2 if null velocity (tau)
% E - 2D matrix of error in the solution
%                                 
%
% Martino Alessandrini
% Post Doc CREATIS
% June 2011

[v E] = initialize_variables(M_Ix2,options.mode);
[N M P] = size(v); 

if maskNoise == 0
    dt_threshold = options.dt_threshold*median(M_It.m_00(:)); 
    idxSufficientDt = verify_dt_condition(M_It,options.mode,dt_threshold);
else
    idxSufficientDt = verify_dt_condition(M_It,options.mode,maskNoise);
end
    
idxEval = find(mask & idxSufficientDt);
l1norm_wJ = norm(wJ,1);

clear mask %idxSufficientDt

for ii = 1:numel(idxEval) % scan each pixel
    [A_ii b_ii] = compute_system_terms(...
        M_Ix2,M_Ixy,M_Iy2,M_Ixt,M_Iyt,...
        idxEval(ii),options.mode);
    [I,J] = ind2sub([N M],idxEval(ii));  

    testpassed = stability_test(A_ii,options);
    if ~testpassed % v remains zero
        E(I,J) = 100;
        
    else % solve system
        v_sol = pinv(A_ii)*b_ii;
%         v_sol = A_ii\b_ii;
        
        err = norm(A_ii*v_sol-b_ii,2)/l1norm_wJ;
        v(I,J,:) = v_sol;
        E(I,J) = err;
    end
    debug();
end
debug();    
    
function [v E] = initialize_variables(mom,mode)
switch mode
    case 'lucas_kanade'
        [N M] = size(mom.m_00);
        P = 2;
    case 'spatial_affine'
        [N M] = size(mom.m_00);
        P = 6;
    case 'spatio_temporal_affine'
        [N M] = size(mom.m_000);
        P = 8;
end
v = zeros(N,M,P);
E = zeros(N,M);
