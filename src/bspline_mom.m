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

function M = bspline_mom(I,mode,H,wfilt0,wfilt1,wfilt2,symmetry)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiscale moments computation using cubic b-spline interpolation. As
% described in the paper:
% IEEE Trans Image Process. 2004 Apr;13(4):484-95. "Multiresolution moment 
% filters: theory and applications." Sühling M, Arigovindan M, Hunziker P, 
% Unser M.
%
% INPUT DATA
% image_feature -  is the feature I want to compute the moments of. In
%                  order to compute time moments it must be a 3D matrix 
%                  containing 3 frames. Moments are computed for the
%                  central one
% 
%  window - window for local moments (cubic b-spline in this case)
%  J - coarsest scale
%
% OUTPUT DATA
% M - cell array of moments. M{j+1} contains the moments image at scale j.
%     Finer scale is j = 0,..,J. Note that each image is scaled by a factor
%     2 with respect to the finer scale one.
%
% Implementation by:
% Martino Alessandrini
% PostDoc CREATIS
% 14/6/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% filters for 0-order moments computation. 
% m_j_pqr denotes moment at scale j of order p, q and r in x, y and t
% direction respectively.


% symmetry input is [symmetry_x simmetry_y]
% it defines the symmetry of the input signal
M = cell(size(H));
% finest scale moments
M{1} = compute_zero_moments(I,mode,wfilt0,wfilt1,wfilt2,symmetry);
% padding rule for each moment

J = numel(H)-1;
if J > 0
    for jj = 1:J
        M{jj+1} = moments_pyramid(M{jj},H{jj},mode,symmetry);
        debug();
    end
end
debug()


    


