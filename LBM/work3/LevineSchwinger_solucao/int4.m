function stuff = int1(x,z)
% int4.m
%
% This Matlab script defines one of the integral equations
% from Eq. (VI.4) of Levine and Schwinger and is used in
% levine.m to calculate the end correction vs. z=ka for
% sound wave reflection at the open end of a cylindrical
% pipe.
%
% By Gary P. Scavone, CCRMA, Stanford University, 1999.

stuff = log(pi*abs(besselj(1,x)).*(besselj(1,x).^2 + ...
	bessely(1,x).^2).^(0.5))./(x.*(z^2 - x.^2).^(1/2));
