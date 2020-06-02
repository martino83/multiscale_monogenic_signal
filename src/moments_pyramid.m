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

function [M_coarser] = moments_pyramid(M_finer,H_finer,mode,symmetry)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the moments an finer scale compute the ones at the coarser scale
% using the a trous algorithm in:
% IEEE Trans Image Process. 2004 Apr;13(4):484-95. "Multiresolution moment 
% filters: theory and applications." Sühling M, Arigovindan M, Hunziker P, 
% Unser M.
%
% INPUT DATA
% M_finer - image moments matrix at coarser scale
% H_finer - filters for passing from one scale to the clarser one
% J - considered scale
%
% OUTPUT DATA
% M_coearser - moments at next scale
% H_coarser - filters for next scale
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch mode
    case 'spatial_affine_full'
        
        pad_00_x = padding_rule(0,symmetry{1});
        pad_01_x = padding_rule(0,symmetry{1});
        pad_02_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        pad_01_y = padding_rule(1,symmetry{2});
        pad_02_y = padding_rule(2,symmetry{2});
        pad_10_y = padding_rule(0,symmetry{2});
        pad_20_y = padding_rule(0,symmetry{2});
        pad_11_x = padding_rule(1,symmetry{1});
        pad_11_y = padding_rule(1,symmetry{2});
        pad_10_x = padding_rule(1,symmetry{1});
        pad_20_x = padding_rule(2,symmetry{1});
        
        m_00_x = filter_and_downsample_x(M_finer.m_00,H_finer.h_00,pad_00_x);
        m_01_x = filter_and_downsample_x(M_finer.m_01,H_finer.h_00,pad_01_x);
        m_02_x = filter_and_downsample_x(M_finer.m_02,H_finer.h_00,pad_02_x);

        % implementing filter bank
        M_coarser.m_00 = filter_and_downsample_y(m_00_x,H_finer.h_00,pad_00_y);
        M_coarser.m_01 = filter_and_downsample_y(m_00_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_01_x,H_finer.h_11,pad_01_y);

        M_coarser.m_02 = filter_and_downsample_y(m_00_x,H_finer.h_20,pad_00_y) + ...
            filter_and_downsample_y(m_01_x,H_finer.h_21,pad_01_y) + ...
            filter_and_downsample_y(m_02_x,H_finer.h_22,pad_02_y);

        clear m_00_x m_01_x m_02_x    
        
        %
        m_00_y = filter_and_downsample_y(M_finer.m_00,H_finer.h_00,pad_00_y);
        m_10_y = filter_and_downsample_y(M_finer.m_10,H_finer.h_00,pad_10_y);
        m_20_y = filter_and_downsample_y(M_finer.m_20,H_finer.h_00,pad_20_y);

        % implementing filter bank
        
        M_coarser.m_10 = filter_and_downsample_x(m_00_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_10_y,H_finer.h_11,pad_10_x);

        M_coarser.m_20 = filter_and_downsample_x(m_00_y,H_finer.h_20,pad_00_x) + ...
            filter_and_downsample_x(m_10_y,H_finer.h_21,pad_10_x) + ...
            filter_and_downsample_x(m_20_y,H_finer.h_22,pad_20_x);

        clear m_00_y m_10_y m_20_y
        
        M_coarser.m_11 = ...
            filter_and_downsample_xy(M_finer.m_00,H_finer.h_10,H_finer.h_10,pad_00_x,pad_00_y) + ...
            filter_and_downsample_xy(M_finer.m_10,H_finer.h_11,H_finer.h_10,pad_10_x,pad_10_y) + ...
            filter_and_downsample_xy(M_finer.m_01,H_finer.h_10,H_finer.h_11,pad_01_x,pad_01_y) + ...
            filter_and_downsample_xy(M_finer.m_11,H_finer.h_11,H_finer.h_11,pad_11_x,pad_11_y);

    case 'spatial_affine_reduced'
        
        pad_00_x = padding_rule(0,symmetry{1});
        pad_01_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        pad_01_y = padding_rule(1,symmetry{2});
        pad_10_y = padding_rule(0,symmetry{2});
        pad_10_x = padding_rule(1,symmetry{1});
        
        m_00_x = filter_and_downsample_x(M_finer.m_00,H_finer.h_00,pad_00_x);
        m_01_x = filter_and_downsample_x(M_finer.m_01,H_finer.h_00,pad_01_x);

        % implementing filter bank
        
        M_coarser.m_00 = filter_and_downsample_y(m_00_x,H_finer.h_00,pad_00_y);
        M_coarser.m_01 = filter_and_downsample_y(m_00_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_01_x,H_finer.h_11,pad_01_y);
        clear m_00_x m_01_x
        %
        
        m_00_y = filter_and_downsample_y(M_finer.m_00,H_finer.h_00,pad_00_y);
        m_10_y = filter_and_downsample_y(M_finer.m_10,H_finer.h_00,pad_10_y);
        
        
        % implementing filter bank
        M_coarser.m_10 = filter_and_downsample_x(m_00_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_10_y,H_finer.h_11,pad_10_x);
        %
        clear m_00_y m_10_y
        
    case 'spatio_temporal_affine_full'
        
        pad_00_x = padding_rule(0,symmetry{1});
        pad_01_x = padding_rule(0,symmetry{1});
        pad_02_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        pad_01_y = padding_rule(1,symmetry{2});
        pad_02_y = padding_rule(2,symmetry{2});
        pad_10_y = padding_rule(0,symmetry{2});
        pad_20_y = padding_rule(0,symmetry{2});
        pad_11_x = padding_rule(1,symmetry{1});
        pad_11_y = padding_rule(1,symmetry{2});
        pad_10_x = padding_rule(1,symmetry{1});
        pad_20_x = padding_rule(2,symmetry{1});
        
        % t=0
        m_000_x = filter_and_downsample_x(M_finer.m_000,H_finer.h_00,pad_00_x);
        m_010_x = filter_and_downsample_x(M_finer.m_010,H_finer.h_00,pad_01_x);
        m_020_x = filter_and_downsample_x(M_finer.m_020,H_finer.h_00,pad_02_x);

        % implementing filter bank
        M_coarser.m_000 = filter_and_downsample_y(m_000_x,H_finer.h_00,pad_00_y);
        M_coarser.m_010 = filter_and_downsample_y(m_000_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_010_x,H_finer.h_11,pad_01_y);

        M_coarser.m_020 = filter_and_downsample_y(m_000_x,H_finer.h_20,pad_00_y) + ...
            filter_and_downsample_y(m_010_x,H_finer.h_21,pad_01_y) + ...
            filter_and_downsample_y(m_020_x,H_finer.h_22,pad_02_y);

        clear m_000_x m_010_x m_020_x    
        
        %
        m_000_y = filter_and_downsample_y(M_finer.m_000,H_finer.h_00,pad_00_y);
        m_100_y = filter_and_downsample_y(M_finer.m_100,H_finer.h_00,pad_10_y);
        m_200_y = filter_and_downsample_y(M_finer.m_200,H_finer.h_00,pad_20_y);

        % implementing filter bank
        
        M_coarser.m_100 = filter_and_downsample_x(m_000_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_100_y,H_finer.h_11,pad_10_x);

        M_coarser.m_200 = filter_and_downsample_x(m_000_y,H_finer.h_20,pad_00_x) + ...
            filter_and_downsample_x(m_100_y,H_finer.h_21,pad_10_x) + ...
            filter_and_downsample_x(m_200_y,H_finer.h_22,pad_20_x);

        clear m_000_y m_100_y m_200_y
        
        M_coarser.m_110 = ...
            filter_and_downsample_xy(M_finer.m_000,H_finer.h_10,H_finer.h_10,pad_00_x,pad_00_y) + ...
            filter_and_downsample_xy(M_finer.m_100,H_finer.h_11,H_finer.h_10,pad_10_x,pad_10_y) + ...
            filter_and_downsample_xy(M_finer.m_010,H_finer.h_10,H_finer.h_11,pad_01_x,pad_01_y) + ...
            filter_and_downsample_xy(M_finer.m_110,H_finer.h_11,H_finer.h_11,pad_11_x,pad_11_y);

        % t=1
        m_001_x = filter_and_downsample_x(M_finer.m_001,H_finer.h_00,pad_00_x);
        m_011_x = filter_and_downsample_x(M_finer.m_011,H_finer.h_00,pad_01_x);
        
        % implementing filter bank
        M_coarser.m_001 = filter_and_downsample_y(m_001_x,H_finer.h_00,pad_00_y);
        M_coarser.m_011 = filter_and_downsample_y(m_001_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_011_x,H_finer.h_11,pad_01_y);

        clear m_001_x m_011_x    
        
        %
        m_001_y = filter_and_downsample_y(M_finer.m_001,H_finer.h_00,pad_00_y);
        m_101_y = filter_and_downsample_y(M_finer.m_101,H_finer.h_00,pad_10_y);
        
        % implementing filter bank
        
        M_coarser.m_101 = filter_and_downsample_x(m_001_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_101_y,H_finer.h_11,pad_10_x);

        clear m_001_y m_101_y
    
        pad_00_x = padding_rule(0,symmetry{1});
        pad_01_x = padding_rule(0,symmetry{1});
        pad_02_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        pad_01_y = padding_rule(1,symmetry{2});
        pad_02_y = padding_rule(2,symmetry{2});
        pad_10_y = padding_rule(0,symmetry{2});
        pad_20_y = padding_rule(0,symmetry{2});
        pad_11_x = padding_rule(1,symmetry{1});
        pad_11_y = padding_rule(1,symmetry{2});
        pad_10_x = padding_rule(1,symmetry{1});
        pad_20_x = padding_rule(2,symmetry{1});
        
        % t=0
        m_000_x = filter_and_downsample_x(M_finer.m_000,H_finer.h_00,pad_00_x);
        m_010_x = filter_and_downsample_x(M_finer.m_010,H_finer.h_00,pad_01_x);
        m_020_x = filter_and_downsample_x(M_finer.m_020,H_finer.h_00,pad_02_x);

        % implementing filter bank
        M_coarser.m_000 = filter_and_downsample_y(m_000_x,H_finer.h_00,pad_00_y);
        M_coarser.m_010 = filter_and_downsample_y(m_000_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_010_x,H_finer.h_11,pad_01_y);

        M_coarser.m_020 = filter_and_downsample_y(m_000_x,H_finer.h_20,pad_00_y) + ...
            filter_and_downsample_y(m_010_x,H_finer.h_21,pad_01_y) + ...
            filter_and_downsample_y(m_020_x,H_finer.h_22,pad_02_y);

        clear m_000_x m_010_x m_020_x    
        
        %
        m_000_y = filter_and_downsample_y(M_finer.m_000,H_finer.h_00,pad_00_y);
        m_100_y = filter_and_downsample_y(M_finer.m_100,H_finer.h_00,pad_10_y);
        m_200_y = filter_and_downsample_y(M_finer.m_200,H_finer.h_00,pad_20_y);

        % implementing filter bank
        
        M_coarser.m_100 = filter_and_downsample_x(m_000_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_100_y,H_finer.h_11,pad_10_x);

        M_coarser.m_200 = filter_and_downsample_x(m_000_y,H_finer.h_20,pad_00_x) + ...
            filter_and_downsample_x(m_100_y,H_finer.h_21,pad_10_x) + ...
            filter_and_downsample_x(m_200_y,H_finer.h_22,pad_20_x);

        clear m_000_y m_100_y m_200_y
        
        M_coarser.m_110 = ...
            filter_and_downsample_xy(M_finer.m_000,H_finer.h_10,H_finer.h_10,pad_00_x,pad_00_y) + ...
            filter_and_downsample_xy(M_finer.m_100,H_finer.h_11,H_finer.h_10,pad_10_x,pad_10_y) + ...
            filter_and_downsample_xy(M_finer.m_010,H_finer.h_10,H_finer.h_11,pad_01_x,pad_01_y) + ...
            filter_and_downsample_xy(M_finer.m_110,H_finer.h_11,H_finer.h_11,pad_11_x,pad_11_y);

%         t = 2
        M_coarser.m_002 = filter_and_downsample_xy(...
            M_finer.m_002,H_finer.h_00,H_finer.h_00,pad_00_x,pad_00_y);
        
        % t=1
        m_001_x = filter_and_downsample_x(M_finer.m_001,H_finer.h_00,pad_00_x);
        m_011_x = filter_and_downsample_x(M_finer.m_011,H_finer.h_00,pad_01_x);
        
        % implementing filter bank
        M_coarser.m_001 = filter_and_downsample_y(m_001_x,H_finer.h_00,pad_00_y);
        M_coarser.m_011 = filter_and_downsample_y(m_001_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_011_x,H_finer.h_11,pad_01_y);

        clear m_001_x m_011_x    
        
        %
        m_001_y = filter_and_downsample_y(M_finer.m_001,H_finer.h_00,pad_00_y);
        m_101_y = filter_and_downsample_y(M_finer.m_101,H_finer.h_00,pad_10_y);
        
        % implementing filter bank
        
        M_coarser.m_101 = filter_and_downsample_x(m_001_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_101_y,H_finer.h_11,pad_10_x);
        clear m_001_y m_101_y
        
         

    case 'spatio_temporal_affine_reduced'
        
        pad_00_x = padding_rule(0,symmetry{1});
        pad_01_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        pad_01_y = padding_rule(1,symmetry{2});
        pad_10_y = padding_rule(0,symmetry{2});
        pad_10_x = padding_rule(1,symmetry{1});
        
        % t=0
        m_000_x = filter_and_downsample_x(M_finer.m_000,H_finer.h_00,pad_00_x);
        m_010_x = filter_and_downsample_x(M_finer.m_010,H_finer.h_00,pad_01_x);
        
        % implementing filter bank
        M_coarser.m_000 = filter_and_downsample_y(m_000_x,H_finer.h_00,pad_00_y);
        M_coarser.m_010 = filter_and_downsample_y(m_000_x,H_finer.h_10,pad_00_y) + ...
            filter_and_downsample_y(m_010_x,H_finer.h_11,pad_01_y);

        clear m_000_x m_010_x
        %
        m_000_y = filter_and_downsample_y(M_finer.m_000,H_finer.h_00,pad_00_y);
        m_100_y = filter_and_downsample_y(M_finer.m_100,H_finer.h_00,pad_10_y);

        % implementing filter bank
        
        M_coarser.m_100 = filter_and_downsample_x(m_000_y,H_finer.h_10,pad_00_x) + ...
            filter_and_downsample_x(m_100_y,H_finer.h_11,pad_10_x);

        clear m_000_y m_100_y
        
        % t=1
        m_001_x = filter_and_downsample_x(M_finer.m_001,H_finer.h_00,pad_00_x);
        
        % implementing filter bank
        M_coarser.m_001 = filter_and_downsample_y(m_001_x,H_finer.h_00,pad_00_y);
        
        clear m_001_x    
        
        
    case {'time_average','lucas_kanade'}
        pad_00_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        
        M_coarser.m_00 = filter_and_downsample_xy(...
            M_finer.m_00,H_finer.h_00,H_finer.h_00,pad_00_x,pad_00_y);

    case 'spatio_temporal_time_average'
        pad_00_x = padding_rule(0,symmetry{1});
        pad_00_y = padding_rule(0,symmetry{2});
        
        M_coarser.m_000 = filter_and_downsample_xy(...
            M_finer.m_000,H_finer.h_00,H_finer.h_00,pad_00_x,pad_00_y);
        
end
%**************************************************************************

function Ifilt = filter_and_downsample_x(I,v,pad_x)

npad = floor(numel(v)/2);
        
switch pad_x
    case 'mirr'
        I = pad_image_feature(I,[npad 0],{'mirr','mirr'});
        Ifilt = conv2(1,v,I,'valid');
        Ifilt = Ifilt(:,1:2:end);
        
    case 'anti_mirr'
        I = pad_image_feature(I,[npad 0],{'anti_mirr','mirr'});
        Ifilt = conv2(1,v,I,'valid');
        Ifilt = Ifilt(:,1:2:end);
        
    case 'zero'
        Ifilt = conv2(1,v,I,'same');
        Ifilt = Ifilt(:,1:2:end);
        
end


function Ifilt = filter_and_downsample_y(I,v,pad_y)

npad = floor(numel(v)/2);
        
switch pad_y
    case 'mirr'
        I = pad_image_feature(I,[0 npad],{'mirr','mirr'});
        Ifilt = conv2(v,1,I,'valid');
        Ifilt = Ifilt(1:2:end,:);
        
    case 'anti_mirr'
        I = pad_image_feature(I,[0 npad],{'mirr','anti_mirr'});
        Ifilt = conv2(v,1,I,'valid');
        Ifilt = Ifilt(1:2:end,:);
        
    case 'zero'
        Ifilt = conv2(v,1,I,'same');
        Ifilt = Ifilt(1:2:end,:);
        
end


function Ifilt = filter_and_downsample_xy(I,h,v,pad_x,pad_y)
Ifilt = filter_and_downsample_x(I,h,pad_x);
Ifilt = filter_and_downsample_y(Ifilt,v,pad_y);

function moment_pad = padding_rule(moment_order,input_symmetry)
% padding rule single direction

switch input_symmetry
    case 'mirr'
        if mod(moment_order,2) == 0
            moment_pad = 'mirr';
        else
            moment_pad = 'anti_mirr';
        end
        
    case 'anti_mirr'
        if mod(moment_order,2) == 0
            moment_pad = 'anti_mirr';
        else
            moment_pad = 'mirr';
        end
        
    case 'zero'
        moment_pad = 'zero';
                
end
        

