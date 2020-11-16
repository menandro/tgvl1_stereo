% Solves the warping of the image using the optical flow

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

function Iw = utils_warping(I1, u, v, mode) %input - I,u,v %output - Iw
    [height,width, ~] = size(I1);
    idx = double(repmat((1:width), height,1));
    idy = double(repmat((1:height)', 1,width));
    Idxx = idx + u;
    Idyy = idy + v;
    Idxx(Idxx>width) = width;
    Idxx(Idxx<1) = 1;
    Idyy(Idyy>height) = height;
    Idyy(Idyy<1) = 1;

    if isequal(mode,'linear')
        Iw = interp2(I1, Idxx, Idyy, 'linear');
    elseif isequal(mode,'nearest')
        Iw = interp2(double(I1),Idxx,Idyy,'linear');
        Iw = im2bw(Iw,0.9);
    elseif isequal(mode,'prionearest')
        Iw = interp2(double(I1),Idxx,Idyy,'nearest');
        Iw = round(Iw);
    elseif isequal(mode,'truenearest');
        Iw = interp2(double(I1),Idxx,Idyy,'nearest');
    elseif isequal(mode,'cubic')
        Iw = interp2(I1,Idxx,Idyy,'cubic');
    elseif isequal(mode,'spline')
        Iw = interp2(I1,Idxx,Idyy,'spline');
    elseif isequal(mode,'zeronearest');
        Iw = interp2(double(I1),Idxx,Idyy,'nearest',1);
    end


