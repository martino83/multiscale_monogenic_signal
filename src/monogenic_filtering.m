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
% Computation of the monogenic signal.
% Input
% im                 - input grayvalued image
% filter_family      - sqf filter family. Several families are implemented:
%                      'loggabor'         from [1];
%                      'cauchy'           from [2];
%                      'poisson'          from [3];
%                      'dop'              from [4];
%                      'gaussian1storder' from [2];
% filter_parameters - [wavelength bandwidth_parameter] 
%
% Output
% w, q              - the monogenic signal is w + i*real(q) + j*imag(q);
% grad_q, grad_w    - gradient of filter responses, useful for time and
%                     spatial derivatives of the phase vector.
%
% references:
% [1] P. D. Kovesi, “Invariant measures of image features from phase 
% information,” Ph.D. dissertation, 1996.
% [2] D. Boukerroui, J. A. Noble, and M. Brady, “On the choice of band-pass
% quadrature filters,” Journal of Mathematical Imaging and Vision.
% [3] M. Felsberg and G. Sommer, “The monogenic scale-space: A unifying 
% approach to phase-based image processing in scale-space,” Journal of 
% Mathematical Imaging and Vision, vol. 21, pp. 5–26, 2004.
% [4] L. Wietzke, G. Sommer, and O. Fleischmann, “The geometry of 2d image 
% signals,” in Computer Vision and Pattern Recognition, 2009. CVPR 2009. 
% IEEE Conference on, june 2009, pp. 1690 –1697.
%------------------------------------------------------------------------
function [w q grad_w grad_q] = ...
    monogenic_filtering(im,filter_family,filter_parameters)


IM = fft2(im);
clear im

[H HR ux uy] = quadrature_filters(size(IM),...
    filter_family,filter_parameters);

[w q] = apply_filtering(IM,H,HR);
[w_x q_x] = apply_filtering(2j*pi*ux.*IM,H,HR);
[w_y q_y] = apply_filtering(2j*pi*uy.*IM,H,HR);
grad_w = cat(3,w_x,w_y);
grad_q = cat(3,q_x,q_y);

end

function [w q] = apply_filtering(IM,H,HR)

w = real(ifft2(H.*IM));
q = ifft2(HR.*H.*IM);

end




