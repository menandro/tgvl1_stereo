%% TGV-L1 stereo for rectified images
% Solves the disparity u between two images (assumes rectified images)

%   Author: Menandro Roxas
%   Contact: menandro.roxas@gmail.com
%   $Date: 2020-11-16 (Monday, 16 Nov 2020) $

% Copyright 2020, Menandro Roxas.
%
%                         All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
%
% THE AUTHOR DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 

close all
clear all

im1 = imread('input1.png');
im2 = imread('input2.png');
[height_max, width_max, ch] = size(im1);

%% PARAMETERS - CHANGE IF YOU WANT
alpha0 = 5.0;
alpha1 = 1.1;
timestep_lambda = 1;
tensor_ab = [4 0.2];
l = timestep_lambda;

eta_p = 3;
eta_q = 2;
lambda = 3.0;
innerIter = 5; % Number of iterations of TGV-L1 per warping
warpIter = 20; % Number of warping per pyramid level

nPyr = int8(60); % Max number of pyramids
fScale = 0.9; % Scaling of preceding levels %0.9
minWidth = 16; % Minimum width of the pyramid

flowscale = 100; % For visualization only

%% Construct image pyramid
if ch==3
    gray1 = rgb2gray(im1);
    gray2 = rgb2gray(im2);
else
    gray1 = im1;
    gray2 = im2;
end
Isub1 = double(gray1)/256.0;
Isub2 = double(gray2)/256.0;
[I1, nLevel, sizes] = utils_pyramid(Isub1, nPyr, fScale, minWidth);
[I2, ~, ~]          = utils_pyramid(Isub2, nPyr, fScale, minWidth);

%% Iteration
fprintf('Iterating');
for L=nLevel:-1:1
    [M, N] = size(I1{L});
    fprintf('.');
    curSize = sizes{L};

    % Initialize / Upscale variables
    if L==nLevel
        u = deal(zeros(curSize));
        u_ = u;
    else
        % Propagate to the next level
        u =(1/fScale)*imresize(u,curSize);
        u_ =(1/fScale)*imresize(u_,curSize);
    end

    % Calculate diffusion tensor
    [tensor, ~, ~] = utils_calc_tensor(I1{L}, tensor_ab, 2);
    a = tensor{1};
    b = tensor{2};
    c = tensor{3};
    eta_u = (a.^2 + b.^2 + 2*c.^2 + (a+c).^2 + (b+c).^2).*(alpha1.^2);% + 0*w.^2;
    eta_v = zeros(M,N,2);
    eta_v(:,:,1) = (alpha1.^2).*(b.^2 + c.^2) + 4*alpha0.^2;
    eta_v(:,:,2) = (alpha1.^2).*(a.^2 + c.^2) + 4*alpha0.^2;


    for warp = 1:warpIter
        % Calculate warping vectors
        warpX = u_;
        warpY = zeros(curSize);

        % Warp I2
        Iwarp = utils_warping(I2{L}, warpX, warpY, 'linear');
        [Iu, Iz] = utils_compute_derivative(I1{L}, Iwarp, false);

        %%%%%%%%%%%%%%%%%
        %TGVL1 iteration
        %%%%%%%%%%%%%%%%%
        % Initialize iteration variables
        p = zeros(M,N,2);      
        q = zeros(M,N,4); 
        du = deal(zeros(curSize));
        
        v = zeros(M,N,2);
        v_ = v;
        grad_v = zeros(M,N,4);

        qc = zeros(M, N, 4);
        div_q = zeros(M, N, 2);
        dq_tensor = zeros(M, N, 2);

        tau = 1;
        sigma = 1/tau;

        for k = 1:innerIter
            % update timesteps
            if(sigma < 1000)
                mu = 1/sqrt(1+0.7*tau*l);
            else
                mu = 1;
            end

            % --------- Update dual variables tgv ---------------
            u_x = utils_dxp(u_) - v_(:,:,1);
            u_y = utils_dyp(u_) - v_(:,:,2);
            du_tensor_x = a.*u_x + c.*u_y; 
            du_tensor_y = c.*u_x + b.*u_y;
            p(:,:,1) = p(:,:,1) + (alpha1*sigma/eta_p).*du_tensor_x;
            p(:,:,2) = p(:,:,2) + (alpha1*sigma/eta_p).*du_tensor_y;
            % projection
            reprojection = max(1.0, sqrt(p(:,:,1).^2 + p(:,:,2).^2));
            p(:,:,1) = p(:,:,1)./reprojection;
            p(:,:,2) = p(:,:,2)./reprojection;

            grad_v(:,:,1) = utils_dxp(v_(:,:,1));
            grad_v(:,:,2) = utils_dyp(v_(:,:,2));
            grad_v(:,:,3) = utils_dyp(v_(:,:,1));
            grad_v(:,:,4) = utils_dxp(v_(:,:,2));
            q = q + alpha0*sigma/eta_q.*grad_v;
            % projection
            reproject = max(1.0, sqrt(q(:,:,1).^2 + q(:,:,2).^2 + q(:,:,3).^2 + q(:,:,4).^2));
            q(:,:,1) = q(:,:,1)./reproject;
            q(:,:,2) = q(:,:,2)./reproject;
            q(:,:,3) = q(:,:,3)./reproject;
            q(:,:,4) = q(:,:,4)./reproject;


            % --------- Update primal variables L1 using thresholding ---------------
            Tp1 = a.*p(:,:, 1) + c.*p(:,:,2);
            Tp2 = c.*p(:,:, 1) + b.*p(:,:,2);
            div_p = utils_dxm(Tp1) + utils_dym(Tp2);
            tau_eta_u = tau./eta_u;
            uhat = u_ + tau_eta_u.*div_p;

            dun = uhat - u;
            rho = Iu.*dun + Iz;
            upper = lambda*tau_eta_u.*(Iu.*Iu);
            lower = -lambda*tau_eta_u.*(Iu.*Iu);
            case1 = (rho <= upper) & (rho >= lower);
            case2 = (rho < lower);
            case3 = (rho > upper);
            case4 = (Iu == 0);
            du(case1) = dun(case1) - rho(case1)./Iu(case1);
            du(case2) = dun(case2) + (lambda*tau_eta_u(case2)).*Iu(case2);
            du(case3) = dun(case3) - (lambda*tau_eta_u(case3)).*Iu(case3);
            du(case4) = dun(case4);

            u = u + du;

            qc(:,:,1) = [q(:,1:end-1, 1), zeros(M,1)];
            qc(:,:,2) = [q(1:end-1,:, 2); zeros(1,N)];
            qc(:,:,3) = [q(1:end-1,:, 3); zeros(1,N)];
            qc(:,:,4) = [q(:,1:end-1, 4), zeros(M,1)];

            qw_x = [zeros(M,1,1), q(:,1:end-1,1)];
            qw_w = [zeros(M,1,1), q(:,1:end-1,4)];
            qn_y = [zeros(1,N,1); q(1:end-1,:,2)];
            qn_z = [zeros(1,N,1); q(1:end-1,:,3)];

            div_q(:,:,1) = (qc(:,:,1) - qw_x) + (qc(:,:,3) - qn_z);
            div_q(:,:,2) = (qc(:,:,4) - qw_w) + (qc(:,:,2) - qn_y);

            dq_tensor(:,:,1) = a.*p(:,:,1) + c.*p(:,:,2); 
            dq_tensor(:,:,2) = c.*p(:,:,1) + b.*p(:,:,2);

            v = v_ + tau./eta_v.*(alpha1.*dq_tensor + alpha0.*div_q);

            u_ = u + mu.*(u - u_);
            v_ = v + mu.*(v - v_);

            sigma = sigma/mu;
            tau = tau*mu;
        end
        
        imshow(abs(u)/(flowscale/L), [0 1], 'ColorMap', jet);
    end
    %%%%%%%%%%%%%%%%%
    %TGVL1 end
    %%%%%%%%%%%%%%%%%

end
