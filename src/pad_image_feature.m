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

function [Ipad convmode] = pad_image_feature(I,padvalue,symmetry)

[npad_x npad_y] = deal(padvalue(1),padvalue(2)); % pixels to add in x y direction
[sym_x sym_y] = deal(symmetry{1},symmetry{2}); % padding rule

% paddind x direction
Ipad = I;

if npad_x > 0
    switch sym_x
        case 'anti_mirr'
            Ipad = wextend('2D','asymw',I,[0 npad_x]);
            convmode = 'valid';
        case 'mirr'
            Ipad = wextend('2D','symw',I,[0 npad_x]);
            convmode = 'valid';
        case 'zero'
            convmode = 'same';
    end
end

if npad_y > 0
    switch sym_y
        case 'anti_mirr'
            Ipad = wextend('2D','asymw',Ipad,[npad_y 0]);
            convmode = 'valid';
        case 'mirr'
            Ipad = wextend('2D','symw',Ipad,[npad_y 0]);
            convmode = 'valid';
        case 'zero'
            convmode = 'same';
            
    end
end
