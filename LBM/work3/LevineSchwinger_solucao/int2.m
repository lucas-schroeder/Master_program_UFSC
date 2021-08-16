function stuff = int2(x,z)
% int2.m
%
% This Matlab script defines one of the integral equations
% from Eq. (VI.4) of Levine and Schwinger and is used in
% levine.m to calculate the end correction vs. z=ka for
% sound wave reflection at the open end of a cylindrical
% pipe.
%
% By Gary P. Scavone, CCRMA, Stanford University, 1999.

stuff = log(1./(2*besseli(1,x).*besselk(1,x)))...
	./(x.*(x.^2 + z^2).^(0.5));
