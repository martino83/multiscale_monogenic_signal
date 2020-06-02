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


function [A theta cosinus sinus phi nu] = ...
    monogenic_features_comp(p,q,grad_p,grad_q,parms)

%%% features computation
A = sqrt(p.^2 + abs(q).^2); % enveloppe

%%% orientation computation
% window for computation (ignored in the point case)
if strcmp(parms.orient_mode,'robust')
    s = parms.sigma;
    x = -floor(4*s):floor(4*s);
    w = exp(-x.^2/(s/2)); w = w./sum(w(:));
else
    w = 0;
end

[theta cosinus sinus] = orientation_estimate(q,w,parms.orient_mode);

%%% (wrapped) phase computation
r = real(q).*cosinus + imag(q).*sinus;
phi = angle(p + 1i*sign(real(q)).*abs(r));

%%% frequency computation
nu = compute_frequency(p,q,grad_p,grad_q,...
    cosinus,sinus,parms.freq_mode);
end

function nu = compute_frequency(p,q,grad_p,grad_q,cosinus,sinus,mode)

switch mode
    case 'robust'
        
        eps = 1e-7;
        r1 = real(q);
        r2 = imag(q);
        w = grad_p(:,:,1) + 1i*grad_p(:,:,2);
        r = real(grad_q(:,:,1)) + 1i*imag(grad_q(:,:,2));
        amp = (p.^2 + abs(q).^2);
        
        q = r2.*sinus + r1.*cosinus;
        nu = - q.*(real(w).*cosinus + imag(w).*sinus) + ...
            p.*(real(r)+imag(r));
        nu = nu./(amp+eps);
        nu = max(nu,0);
                
    otherwise
        eps = 1e-7;
        r1 = real(q);
        r2 = imag(q);
        w = grad_p(:,:,1) + 1i*grad_p(:,:,2);
        r = real(grad_q(:,:,1)) + 1i*imag(grad_q(:,:,2));
        amp = (p.^2 + abs(q).^2);
                
        nu = - (r1.*real(w) + r2.*imag(w)) + ...
            p .* (real(r)+imag(r));
        nu = nu./(amp+eps);
        
        nu = max(nu,0);
        
      
end
end



function [theta cosinus sinus] = orientation_estimate(q,w,mode)
w = w(:);

switch mode
    case 'robust'
        
        f1 = real(q);
        f2 = imag(q);
        J11 = f1.^2;
        J12 = f1.*f2;
        J22 = f2.^2;
        clear f1 f2
        
        N = floor(numel(w)/2);        
        J11 = conv2(w,w',padarray(J11,[N N],'symmetric'),'valid');
        J12 = conv2(w,w',padarray(J12,[N N],'symmetric'),'valid');
        J22 = conv2(w,w',padarray(J22,[N N],'symmetric'),'valid');
        
        theta = 1/2*atan2(2*J12,J22-J11);
      
        sinus = cos(theta);
        cosinus = sin(theta);
        
    otherwise % trivial pointwise estimate
        
%         theta=-angle(q.^2)/2;
        theta= - angle(q);
        theta = wrapToPi(2*theta + pi)/2;
        sinus = cos(theta);
        cosinus = sin(theta);
        
end
end 
