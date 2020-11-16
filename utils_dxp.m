function [dx] = utils_dxp(u)
[M N] = size(u);
dx = [u(:,2:end) u(:,end)] - u;