function stuff = int3(x,z)
% int3.m
%
% This Matlab script defines the integral equation from 
% Eq. (VI.5) of Levine and Schwinger and is used in
% levine.m to calculate the reflection coefficient magnitude
% vs. z=ka for sound wave reflection at the open end of a
% cylindrical pipe.
%
% By Gary P. Scavone, CCRMA, Stanford University, 1999.

stuff = atan(besselk(1,x)./(pi*besseli(1,x)))...
	.*(1 - z./((z^2 + x.^2).^(1/2))).*(1./x);
