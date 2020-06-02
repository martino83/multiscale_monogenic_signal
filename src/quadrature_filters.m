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
% computation of sqf for monogenic signal computation
% family - if 'LogGabor' then: parameters = [wavelength, sigmaOnf];
%          if 'Cauchy' then: parameters = [wavelength, a]
%          if 'Poisson' the: parameters = [wavelength]
%          if 'gaussian1storder' then: parameters = [wavelength]
%------------------------------------------------------------------------

function [H HR u1 u2] = quadrature_filters(siz,family,parameters)



%%%% frequency plane
sx = siz(2);
sy = siz(1);
        
if sx == 1
    sx = sy;
    sy = 1;
end
        
if mod(sx,2)
    xrange = (-(sx-1)/2:(sx-1)/2)'/(sx-1);
else
    xrange = (-sx/2:(sx/2-1))'/sx;     
end


if sy == 1
    
    u1 = xrange;
    u2 = 1;
    radius = abs(ifftshift(u1));   % Quadrant shift to put 0 frequency at the corners    
    radius(1,1) = 1;
    HR = []; 
    lp = 1;
    
else


    if mod(sy,2)
        yrange = (-(sy-1)/2:(sy-1)/2)/(sy-1);
    else
        yrange = (-sy/2:(sy/2-1))/sy;     
    end
    [u1,u2] = meshgrid(xrange, yrange);
    u1 = ifftshift(u1);   % Quadrant shift to put 0 frequency at the corners
    u2 = ifftshift(u2);
    radius = sqrt(u1.^2 + u2.^2);    
    radius(1,1) = 1;
    HR = (- 1i*u1 + u2)./radius; 
    lp = lowpassfilter(siz,.4,10);   % Radius .4, 'sharpness' 10

end

switch lower(family)
    case 'loggabor' % kavesi
%       parameters meaning
%       sigmaOnf = 0.75 - 1 octave
%                = 0.55 - 2 octaves
%                = 0.41 - 3 octaves
        
        wavelength = parameters(1);
        sigmaOnf = parameters(2);
        fo = 1.0/wavelength;                  % Centre frequency of filter.
        logGabor = exp((-(log(radius/fo)).^2) / (2*log(sigmaOnf)^2));  
        logGabor = logGabor.*lp;              % Apply low-pass filter
        logGabor(1,1) = 0;                    
        H = logGabor;
        
    case 'cauchy'
%       parameters meaning        
%       a = 1.8597 - 2.5 octaves
        wavelength = parameters(1);
        a = parameters(2);
        fo = 1.0/wavelength;
        sigma = a./(fo);
        Cauchy = (radius).^a.*exp(-sigma*radius).*lp;
        Cauchy = Cauchy./max(Cauchy(:));
        Cauchy(1,1) = 0;
        H = Cauchy;
        
    case 'poisson' % ferlsberg
%       parameters meaning
%       s = ( 1+exp(2*pi*f0) ) / ( exp(2*pi*f0) - 1 )
%       circa 1.5 octaves (fixed)
        wavelength = parameters(1);
        fo = 1./wavelength;
        scale = ( 1+exp(2*pi*fo) ) / ( exp(2*pi*fo) - 1 );
        Poisson = exp(-2*pi*radius*(scale-1)) - 2*exp(-2*pi*radius*(scale)) + ...
            exp(-2*pi*radius*(scale+1)).*lp;
        Poisson = Poisson./max(Poisson(:));
        Poisson(1,1) = 0;
        H = Poisson;
        
    case 'dop'
        lambda = parameters(1);
        ratio = parameters(2);
        if ratio >= 1
            error('ratio must be smaller than 1 !!!')
        end
        s2 = lambda/(2*pi*(ratio-1))*log(ratio);
        s1 = ratio*s2;
        % f0 = 1/(2*pi*(s1-s2)*log(s2/s1))
        DoP = exp(-2*pi*radius*s1) - exp(-2*pi*radius*s2);
        DoP = DoP.*lp;
        DoP(1,1) = 0;
        H = DoP;
        
    case 'gaussian1storder'
%       parameters meaning        
%       2.5 octaves (fixed)
        wavelength = parameters(1);
        sigma = wavelength;
        GaussDer = (radius).*exp(-sigma^2*radius.^2).*lp;
        GaussDer = GaussDer./max(GaussDer(:));
        GaussDer(1,1) = 0;
        H = GaussDer;
        
    otherwise
        error('undefined family');
        
end
