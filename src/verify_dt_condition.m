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

function idx = verify_dt_condition(M_It,mode,noiseLevel)

coef = .3;
switch mode
    case {'lucas_kanade','spatial_affine'}
        if numel(noiseLevel) == 1
            idx = abs(M_It.m_00) > noiseLevel;
        else
            noiseLevel = abs(M_It.m_00(noiseLevel==1));
            noiseLevel = coef*mean(noiseLevel(:));
            idx = abs(M_It.m_00) > noiseLevel;
            debug()
        end
            
    case 'spatio_temporal_affine'
        if numel(noiseLevel) == 1
            idx = abs(M_It.m_000) > noiseLevel;
        else
            noiseLevel = abs(M_It.m_000(noiseLevel==1));
            noiseLevel = coef*mean(noiseLevel(:));
            idx = abs(M_It.m_000) > noiseLevel;
        end
end
        
