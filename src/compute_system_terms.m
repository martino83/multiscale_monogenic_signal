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

function [A b] = compute_system_terms(M_Ix2,M_Ixy,M_Iy2,M_Ixt,...
    M_Iyt,ii,mode)

switch mode
    case 'spatial_affine'
        A = [M_Ix2.m_00(ii) M_Ixy.m_00(ii) M_Ix2.m_10(ii) M_Ix2.m_01(ii) M_Ixy.m_10(ii) M_Ixy.m_01(ii); ...
               0            M_Iy2.m_00(ii) M_Ixy.m_10(ii) M_Ixy.m_01(ii) M_Iy2.m_10(ii) M_Iy2.m_01(ii); ...
               0            0                  M_Ix2.m_20(ii) M_Ix2.m_11(ii) M_Ixy.m_20(ii) M_Ixy.m_11(ii); ...
               0            0                         0           M_Ix2.m_02(ii) M_Ixy.m_11(ii) M_Ixy.m_02(ii); ...
               0            0                         0                  0           M_Iy2.m_20(ii) M_Iy2.m_11(ii); ...
               0            0                         0                  0                  0           M_Iy2.m_02(ii)];
        A = A + triu(A,1)';       
        b = -[M_Ixt.m_00(ii) M_Iyt.m_00(ii) M_Ixt.m_10(ii) M_Ixt.m_01(ii) M_Iyt.m_10(ii) M_Iyt.m_01(ii)]';
    
    case 'lucas_kanade'
        A = [M_Ix2.m_00(ii) M_Ixy.m_00(ii);M_Ixy.m_00(ii) M_Iy2.m_00(ii)];
        b = -[M_Ixt.m_00(ii) M_Iyt.m_00(ii)]';
        
     case 'spatio_temporal_affine'
        A = [M_Ix2.m_000(ii) M_Ixy.m_000(ii) M_Ix2.m_100(ii) M_Ix2.m_010(ii)  M_Ix2.m_001(ii) M_Ixy.m_100(ii) M_Ixy.m_010(ii) M_Ixy.m_001(ii); ...
               0            M_Iy2.m_000(ii)  M_Ixy.m_100(ii) M_Ixy.m_010(ii)  M_Ixy.m_001(ii) M_Iy2.m_100(ii) M_Iy2.m_010(ii) M_Iy2.m_001(ii); ...
               0            0                M_Ix2.m_200(ii) M_Ix2.m_110(ii)  M_Ix2.m_101(ii) M_Ixy.m_200(ii) M_Ixy.m_110(ii) M_Ixy.m_101(ii); ...
               0            0                         0      M_Ix2.m_020(ii)  M_Ix2.m_011(ii) M_Ixy.m_110(ii) M_Ixy.m_020(ii) M_Ixy.m_011(ii); ...
               0            0                         0                  0    M_Ix2.m_002(ii) M_Ixy.m_101(ii) M_Ixy.m_011(ii) M_Ixy.m_002(ii) 
               0            0                         0                  0            0       M_Iy2.m_200(ii) M_Iy2.m_110(ii) M_Iy2.m_101(ii); ...
               0            0                         0                  0            0             0         M_Iy2.m_020(ii) M_Iy2.m_011(ii);
               0            0                         0                  0            0             0                   0     M_Iy2.m_002(ii)];
        A = A + triu(A,1)';       
        b = -[M_Ixt.m_000(ii) M_Iyt.m_000(ii) M_Ixt.m_100(ii) M_Ixt.m_010(ii) M_Ixt.m_001(ii) M_Iyt.m_100(ii) M_Iyt.m_010(ii) M_Iyt.m_001(ii)]';
        
end
