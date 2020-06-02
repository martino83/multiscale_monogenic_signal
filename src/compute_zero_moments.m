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

function M = compute_zero_moments(I,mode,w0,w1,w2,symmetry)

npad = floor(numel(w0)/2); % padded pixels


switch mode
    case {'lucas_kanade','time_average'}
        [Ipad convmode] = pad_image_feature(I,[npad npad],symmetry);
        I = Ipad;
        clear Ipad
        
        m_0_00 = conv2(w0,w0,I,convmode);
        M = struct('m_00',m_0_00);
        
    case 'spatio_temporal_time_average'
        for kk = 1:size(I,3)
            [Ipad(:,:,kk) convmode] = pad_image_feature(I(:,:,kk),[npad npad],symmetry);
        end
        I = Ipad;
        clear Ipad
        m_0_xx0 = zeros(size(I,1),size(I,2));
        
        for kk = 1:size(I,3)
            m_0_xx0 = m_0_xx0 + w0(kk+1)*I(:,:,end-kk+1);
        end
        m_0_000 = conv2(w0,w0,m_0_xx0,convmode);
        M = struct('m_000',m_0_000);
        clear m_0_000
         
    case 'spatial_affine_full'
        [Ipad convmode] = pad_image_feature(I,[npad npad],symmetry);
        % compute moments at scale 0 of order p, q
        I = Ipad;
        clear Ipad
        % y - filtering
        m_0_x0 = conv2(w0,1,I,convmode);
        m_0_x1 = conv2(w1,1,I,convmode);
        m_0_02 = conv2(w2,w0,I,convmode);
        clear I
        
        % x - filtering
        m_0_00 = conv2(1,w0,m_0_x0,convmode);
        m_0_10 = conv2(1,w1,m_0_x0,convmode);
        m_0_20 = conv2(1,w2,m_0_x0,convmode);

        m_0_01 = conv2(1,w0,m_0_x1,convmode);
        m_0_11 = conv2(1,w1,m_0_x1,convmode);

        clear m_0_x0 m_0_x1 w0 w1 w2

        M = struct('m_00',m_0_00,'m_01',m_0_01,...
            'm_10',m_0_10,'m_11',m_0_11,'m_20',m_0_20,...
            'm_02',m_0_02);

        clear m_0_00 m_0_01 m_0_10 m_0_11 m_0_02 m_0_20

    case 'spatial_affine_reduced'
        [Ipad convmode] = pad_image_feature(I,[npad npad],symmetry);
        % compute moments at scale 0 of order p, q and r
        I = Ipad;
        clear Ipad
        % y - filtering
        m_0_x0 = conv2(w0,1,I,convmode);
        m_0_x1 = conv2(w1,1,I,convmode);

        clear I
        % x - filtering
        m_0_01 = conv2(1,w0,m_0_x1,convmode);
        m_0_00 = conv2(1,w0,m_0_x0,convmode);
        m_0_10 = conv2(1,w1,m_0_x0,convmode);

        clear m_0_x0 m_0_x1

        M = struct('m_00',m_0_00,'m_01',m_0_01,...
            'm_10',m_0_10);
        
        clear m_0_000 m_0_001 m_0_010 m_0_100     

%     case 'time_average'
%         m_0_00 = conv2(w0,w0,I,'valid');
%         M = struct('m_00',m_0_00);
%         clear I m_0_00
        
    case 'spatio_temporal_affine_full'
        % compute moments at scale 0 of order p, q
        for kk = 1:size(I,3)
            [Ipad(:,:,kk) convmode] = pad_image_feature(I(:,:,kk),[npad npad],symmetry);
        end
        I = Ipad;
        clear Ipad
        % time integration
        m_0_xx0 = zeros(size(I,1),size(I,2));
        m_0_xx1 = zeros(size(I,1),size(I,2));
        m_0_xx2 = zeros(size(I,1),size(I,2));
        
        for kk = 1:size(I,3)
            m_0_xx0 = m_0_xx0 + w0(kk+1)*I(:,:,end-kk+1);
            m_0_xx1 = m_0_xx1 + w1(kk+1)*I(:,:,end-kk+1);
            m_0_xx2 = m_0_xx2 + w2(kk+1)*I(:,:,end-kk+1);
        end
            
        
        % y - filtering
        m_0_x00 = conv2(w0,1,m_0_xx0,convmode);
        m_0_x10 = conv2(w1,1,m_0_xx0,convmode);
        
        m_0_x01 = conv2(w0,1,m_0_xx1,convmode);
        m_0_x11 = conv2(w1,1,m_0_xx1,convmode);
        
        m_0_020 = conv2(w2,w0,m_0_xx0,convmode);
        m_0_002 = conv2(w0,w0,m_0_xx2,convmode);
        
        clear I m_0_xx0 m_0_xx1 m_0_xx2
        
        % x - filtering
        m_0_000 = conv2(1,w0,m_0_x00,convmode);
        m_0_100 = conv2(1,w1,m_0_x00,convmode);
        m_0_200 = conv2(1,w2,m_0_x00,convmode);

        m_0_001 = conv2(1,w0,m_0_x01,convmode);
        m_0_101 = conv2(1,w1,m_0_x01,convmode);
        
        m_0_010 = conv2(1,w0,m_0_x10,convmode);
        m_0_110 = conv2(1,w1,m_0_x10,convmode);
        m_0_011 = conv2(1,w0,m_0_x11,convmode);
        
        clear m_0_x00 m_0_x10 m_0_x01 m_0_x11 w0 w1 w2

        M = struct('m_000',m_0_000,'m_001',m_0_001,'m_010',m_0_010,...
                    'm_011',m_0_011,'m_100',m_0_100,'m_101',m_0_101,'m_110',m_0_110,...
                    'm_002',m_0_002,'m_020',m_0_020,'m_200',m_0_200);

        clear m_0_000 m_0_001 m_0_010 m_0_011 m_0_100 m_0_101 m_0_110 ...
                    m_0_002 m_0_020 m_0_200
                
    case 'spatio_temporal_affine_reduced'
        for kk = 1:size(I,3)
            [Ipad(:,:,kk) convmode] = pad_image_feature(I(:,:,kk),[npad npad],symmetry);
        end
        I = Ipad;
        clear Ipad
        % compute moments at scale 0 of order p, q
        
        % time integration
        m_0_xx0 = zeros(size(I,1),size(I,2));
        m_0_xx1 = zeros(size(I,1),size(I,2));
        
        for kk = 1:size(I,3)
            m_0_xx0 = m_0_xx0 + w0(kk+1)*I(:,:,end-kk+1);
            m_0_xx1 = m_0_xx1 + w1(kk+1)*I(:,:,end-kk+1);
        end
            
        
        % y - filtering
        m_0_x00 = conv2(w0,1,m_0_xx0,convmode);
        m_0_x10 = conv2(w1,1,m_0_xx0,convmode);
        
        m_0_x01 = conv2(w0,1,m_0_xx1,convmode);
        
        clear I m_0_xx0 m_0_xx1
        
        % x - filtering
        m_0_000 = conv2(1,w0,m_0_x00,convmode);
        m_0_100 = conv2(1,w1,m_0_x00,convmode);
        
        m_0_001 = conv2(1,w0,m_0_x01,convmode);
        m_0_010 = conv2(1,w0,m_0_x10,convmode);
        
        
        clear m_0_x00 m_0_x10 m_0_x01 

        M = struct('m_000',m_0_000,'m_001',m_0_001,'m_010',m_0_010,...
            'm_100',m_0_100);
        clear m_0_000 m_0_001 m_0_010 m_0_100
        
    otherwise
        error('wrong mode option')
end

        
% function w = compute_moment_filter(n,window)
% w = [1 0 -1].^n; w = w.*window;
