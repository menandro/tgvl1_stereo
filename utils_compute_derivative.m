%% Compute 1D derivative for stereo

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

function [Iu, Iz] = utils_compute_derivative(I0, I1, isFisheye)
    mask = [1 -8 0 8 -1]/12.0;
    if isa(I0,'double')==false
        I0 = double(I0);
    end
    if isa(I1,'double')==false
        I1 = double(I1);
    end
    if (~isFisheye)
        Iu = (imfilter(I0, mask, 'replicate') + imfilter(I1, mask, 'replicate'))*0.5; %x derivative
        Iz = (I1-I0);
    else
        Iu = (imfilter(I0, mask, 'replicate') + imfilter(I1, mask, 'replicate'))*0.5; %x derivative
        Iz = (I1-I0);
    end
end