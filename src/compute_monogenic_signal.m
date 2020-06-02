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
% Computation of the monogenic signal and associated image features. 
% [fM feats] = compute_monogenic_signal(im,parms)
% Input:
% im    - input image.
% parms - parameters for monogenic signal computation. It's a structure
%         with the following fields:
%         .filter_family: defines the sqf used for the monogenic signal
%          computation. Default is 'DoP'. type 'help
%          monogenic_filtering' for further details on available families
%         .filter_parameters: contains the parameters (typically center
%          wavelength and bandwidth) of the sqf. type 'help
%          monogenic_filtering' for further details.
%         .orient_mode: can be 'robust' or 'classical' depending on wheter
%          robust or pointwise computation is used.
%         .freq_mode: can be 'robust' or 'classical' depending on wheter
%          robust or pointwise computation is used.
%         .sigma: is the standard deviation of the gaussian kernel used in
%          the robust orientation computation.
% Output:
% fM    - its fiels .p and .q (complex) represent the monogenic signal. In
%         the quaternion formulation it would be: 
%         fM.p + i*real(fM.q) + 1j*imag(fM.q).
% feats - monogenic features:
%         .A: amplitude;
%         .theta: orientation;
%         .cos_theta, .sin_theta: cosinus and sinus of theta;
%         .phi: monogenic phase;
%         .nu: frequency;
%
% An example of usage is given in the 'main_monogenic_signal.m'
%
% The code takes inspiration from other implementation available on the
% internet, namely the one by Michael Felsberg:
% http://www.cvl.isy.liu.se/research/ima/monogenic-signal-and-scale-space
% and the one by Peter Kovesi:
% http://www.csse.uwa.edu.au/~pk/research/matlabfns/
%
% Monogenic features are computed as decribed in:
% /M. Felsberg, "Optical Flow Estimation From Monogenic Phase", First 
% International Workshop, IWCM 2004
%
% The robust orientation and frequency computation follows:
% /Unser, M.; Sage, D.; Van De Ville, D.; , "Multiresolution Monogenic 
% Signal Analysis Using the Riesz–Laplace Wavelet Transform," Image 
% Processing, IEEE Transactions on , vol.18, no.11, pp.2402-2418, Nov. 2009
%------------------------------------------------------------------------

function [fM feats] = compute_monogenic_signal(im,varargin)

if nargin < 2 || isempty(varargin{1})
    parms = 0;
else
    parms = varargin{1};
end

parms = fix_parms(parms);

im = double(im);

[p q grad_p grad_q] = monogenic_filtering(im,...
            parms.filter_family,parms.filter_parameters);
                
fM.p = p; % monogenic signal
fM.q = q;

[A theta cos_theta sin_theta phi nu] = ...
    monogenic_features_comp(p,q,grad_p,grad_q,parms);

feats = struct('A',A,'theta',theta,...
    'phi',phi,...
    'freq',nu,'sin_theta',sin_theta,...
    'cos_theta',cos_theta);


function parmsOut = fix_parms(parmsIn)

if ~isstruct(parmsIn) 
   parmsOut.filter_family = 'DoP';
   parmsOut.filter_parameters = [3 .98];
   parmsOut.orient_mode = 'robust';
   parmsOut.freq_mode = 'robust';
   parmsOut.sigma = .5;
else
    parmsOut = parmsIn;
    if ~isfield(parmsOut,'filter_family'), parmsOut.filter_family = 'DoP'; end
    if ~isfield(parmsOut,'filter_parameters'), parmsOut.filter_parameters = [3 .98]; end
    if ~isfield(parmsOut,'orient_mode'), parmsOut.orient_mode = 'robust'; end   
    if ~isfield(parmsOut,'freq_mode'), parmsOut.freq_mode = 'robust'; end
    if ~isfield(parmsOut,'sigma'), parmsOut.sigma = .5; end
end