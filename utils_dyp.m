function [dy] = utils_dyp(u)
[M N] = size(u);
dy = [u(2:end,:); u(end,:)] - u;