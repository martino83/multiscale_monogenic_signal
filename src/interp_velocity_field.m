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

function [v_transf E_transf] = interp_velocity_field...
    (v_coars,E_coars,final_size)

v_transf = interp_3D_data(v_coars,final_size);
E_transf = interp_3D_data(E_coars,final_size);
E_transf(E_transf<0) = 0;

debug()

function v_transf = interp_3D_data(v_coars,final_size)

[N M P] = size(v_coars);
[nF mF] = deal(final_size(1),final_size(2));
v_transf = zeros(nF,mF,P);
[XI YI] = meshgrid(1:.5:M,1:.5:N); 

for pp = 1:P
     
    v_interp = interp2(v_coars(:,:,pp),XI,YI,'linear',0); 
    
    if mod(nF,2) == 0
%         v_interp = wextend('ar','asymw',v_interp,1,'d');
        v_interp = padarray(v_interp,[1 0],'replicate','post');
    end
    if mod(mF,2) == 0
%         v_interp = wextend('ac','asymw',v_interp,1,'r');
        v_interp = padarray(v_interp,[0 1],'replicate','post');
    end
    
    v_transf(:,:,pp) = v_interp;
        
end
