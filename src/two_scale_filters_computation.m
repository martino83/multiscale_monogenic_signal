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
% The function implements fir filtering the multiresolution moments
% computation. 
% 
% The implementation follows the paper:
% Sühling M, Arigovindan M, Hunziker P, Unser M., "Multiresolution moment 
% filters: theory and applications." EEE Trans Image Process. 2004 Apr;
% 13(4):484-95. 
%------------------------------------------------------------------------

function H = two_scale_filters_computation(J,N,mode)

l = -(N+1)/2:(N+1)/2;
[u2N, c] = u2N_FIR_coefs(N);
h = c*u2N; % h00

h = fliplr(h);
l = fliplr(l);

H = cell(1,J+1);

switch mode
    case 'lucas_kanade'
        H{1}.h_00 = h;        
        for kk = 2:J+1
            H{kk}.h_00 = H{1}.h_00;
        end
    case {'spatial_affine','spatio_temporal_affine'}
        H{1}.h_00 = h;
        H{1}.h_10 = l.*h;
        H{1}.h_11 = h;
        H{1}.h_20 = 1*l.^2.*h;
        H{1}.h_21 = 2*l.*h;
        H{1}.h_22 = h;
        for kk = 2:J+1
            H{kk}.h_00 = H{1}.h_00;
            H{kk}.h_10 = 2^(kk-1)*H{1}.h_10;
            H{kk}.h_11 = H{1}.h_11;
            H{kk}.h_20 = 2^(2*(kk-1))*H{1}.h_20;
            H{kk}.h_21 = 2^(kk-1)*H{1}.h_21;
            H{kk}.h_22 = H{1}.h_22;
        end
end
