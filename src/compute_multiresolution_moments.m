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
% Computation of the image moments at different scales. Ôutputs are cells
% with a number of entries equal to the number of scales. Each entry is an
% image of the computed moment. Image sizes decrease in diadic way.
%
% Monogenic signal computation and spatial and temporal derivatived of the
% phase vector follow the paper:
%
% /M. Felsberg, "Optical Flow Estimation From Monogenic Phase", First 
% International Workshop, IWCM 2004
% 
% Multiresolution moments computation follws the papers:
% 
% /Suhling, M.; Arigovindan, M.; Jansen, C.; Hunziker, P.; Unser, M.; 
% "Myocardial motion analysis from B-mode echocardiograms," Image 
% Processing, IEEE Transactions on , vol.14, no.4, pp.525-536, April 2005
% 
% /Suhling, M.; Arigovindan, M.; Hunziker, P.; Unser, M.; , 
% "Multiresolution moment filters: theory and applications," Image 
% Processing, IEEE Transactions on , vol.13, no.4, pp.484-495, April 2004
%------------------------------------------------------------------------


function [M_Ix2 M_Iy2 M_Ixy M_Ixt M_Iyt M_It] = compute_multiresolution_moments(I,H,...
    options,b,wfilt1,wfilt2)
                       
nframes = size(I,3);         
switch nframes

    case 2 % forward differencing

        [fM_ref feats] = compute_monogenic_signal(I(:,:,1),options.monogenic_parms);
        [fM_new] = compute_monogenic_signal(I(:,:,2),options.monogenic_parms);

        % spatial derivatives
        r1_x = feats.freq.*feats.cos_theta.^2;
        r2_y = feats.freq.*feats.sin_theta.^2;
        r1_y = feats.freq.*feats.sin_theta.*feats.cos_theta;
        r2_x = r1_y;

        % time derivative
        r_diff = phase_vector_derivative(fM_ref,fM_new);

        Ix2 = r1_x;
        Iy2 = r2_y;
        Ixy = r2_x;
        Ixt = real(r_diff);
        Iyt = imag(r_diff);
        It = abs(r_diff);

    case 3

        %%%%%%%%% pahse vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fM_old] = compute_monogenic_signal(I(:,:,1),options.monogenic_parms);
        [fM_ref feats] = compute_monogenic_signal(I(:,:,2),options.monogenic_parms);
        [fM_new] = compute_monogenic_signal(I(:,:,3),options.monogenic_parms);
        
        % spatial derivatives
        r1_x = feats.freq.*feats.cos_theta.^2;
        r2_y = feats.freq.*feats.sin_theta.^2;
        r1_y = feats.freq.*feats.sin_theta.*feats.cos_theta;
        r2_x = r1_y;
        
        % time derivative
        r_diff_forward = phase_vector_derivative(fM_ref,fM_new);
        r_diff_backward = phase_vector_derivative(fM_old,fM_ref);
        r_diff = .5*(r_diff_forward+r_diff_backward);


        Ix2 = r1_x;
        Iy2 = r2_y;
        Ixy = r2_x;
        Ixt = real(r_diff);
        Iyt = imag(r_diff);
        It = abs(r_diff);
        debug();

end        
    


% boundary conditions:
% Assume I is padded with mirror reflection in x and y [mirror_x mirror_y]
% then it is:
% Ix2 -> [mirror_x mirror_y];
% Iy2 -> [mirror_x mirror_y];
% Ixy -> [anti_mirror_x anti_mirror_y];
% Ixt -> [anti_mirror_x mirror_y];
% Iyt -> [mirror_x anti_mirror_y];

% these properties must be respected by the moments:
% even moments -> same symmetry as input
% odd moments -> opposite symmetry as input

w0 = b{1}; % windows for moments computation
w1 = wfilt1{1};
w2 = wfilt2{1};

% it is:
% w0 -> even
% w1 -> odd
% w2 -> even

% {'mirr','mirr'} defines the kind of padding that is applied to the
% original signal

if options.padding
    pad_Ix2 = {'mirr','mirr'};
    pad_Iy2 = {'mirr','mirr'};
    pad_Ixy = {'anti_mirr','anti_mirr'};
    pad_Ixt = {'anti_mirr','anti_mirr'};
    pad_Iyt = {'mirr','anti_mirr'};
    pad_It = {'mirr','mirr'};
    
else
    pad_Ix2 = {'zero','zero'};
    pad_Iy2 = {'zero','zero'};
    pad_Ixy = {'zero','zero'};
    pad_Ixt = {'zero','zero'};
    pad_Iyt = {'zero','zero'};
    pad_It = {'zero','zero'};
end
        
switch options.mode
    case 'lucas_kanade'
        M_Ix2 = bspline_mom(Ix2,'lucas_kanade',H,w0,w1,w2,pad_Ix2);
        M_Iy2 = bspline_mom(Iy2,'lucas_kanade',H,w0,w1,w2,pad_Iy2);
        M_Ixy = bspline_mom(Ixy,'lucas_kanade',H,w0,w1,w2,pad_Ixy);
        M_Ixt = bspline_mom(Ixt,'lucas_kanade',H,w0,w1,w2,pad_Ixt);
        M_Iyt = bspline_mom(Iyt,'lucas_kanade',H,w0,w1,w2,pad_Iyt);
        M_It = bspline_mom(It.^2,'time_average',H,w0,w1,w2,pad_It);
    case 'spatial_affine'
        M_Ix2 = bspline_mom(Ix2,'spatial_affine_full',H,w0,w1,w2,pad_Ix2);
        M_Iy2 = bspline_mom(Iy2,'spatial_affine_full',H,w0,w1,w2,pad_Iy2);
        M_Ixy = bspline_mom(Ixy,'spatial_affine_full',H,w0,w1,w2,pad_Ixy);
        M_Ixt = bspline_mom(Ixt,'spatial_affine_reduced',H,w0,w1,w2,pad_Ixt);
        M_Iyt = bspline_mom(Iyt,'spatial_affine_reduced',H,w0,w1,w2,pad_Iyt); 
        M_It = bspline_mom(It.^2,'time_average',H,w0,w1,w2,pad_It);
end


function r_diff = phase_vector_derivative(fM_t,fM_tp1)

% time derivative
p_new = fM_tp1.p;
q_new = fM_tp1.q;
p_ref = fM_t.p;
q_ref = fM_t.q;

p_diff = p_ref.*p_new + real( q_ref.*conj(q_new) );
q_diff = p_ref.*q_new - q_ref.*p_new;
r_diff = q_diff./(abs(q_diff)+eps).*atan2(abs(q_diff),p_diff);

    