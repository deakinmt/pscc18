function [ lz,lr,lx,R,X ] = lambdas( Z )
%LAMBDAS Summary of this function goes here
%   Detailed explanation goes here
R = real(Z);
X = imag(Z);

lz = R/X;
lr = R/abs(Z);
lx = X/abs(Z);

end

