function [Bz,Br] = B_ring(z,b,Io)
%B_RING Compute analytic value of magnetic field due to a ring of current.
%
%   Bz = B_RING(z,b,I) = 2*pi*b*b*I./(z.*z + b*b)^1.5
%
%   Bz is the magnetic field at a point z on the z-axis due to a ring of
%   radius b in the x-y plane centered on the origin.  A current of
%   magnitude I flows counter-clockwise (as viewed from a point on the +z
%   axis). The array z contains the measurement points.
%
%   [Bz,Br] = B_RING([z,r],b,I) Observation position is at distance r
%   perpendicular to the z-axis. (z and r can be single-column arrays).
%
%   To put B in MKS units, multiply by mu_o/(4*pi).

if (size(z,2) == 2)
  r = z(:,2);
  z = z(:,1);
end

if (~exist('r')) | (size(z,2) == 2)
  Bz = 2*pi*b*b*Io./(z.*z + b*b)^1.5;
  return;
end

% See http://www.netdenizen.com/emagnet/offaxis/iloopoffaxis.htm
if (exist('r'))
  alpha = r/b;
  beta  = z/b;
  gamma = z./r;
  Q     = (1+alpha).^2 + beta.^2;
  k     = sqrt(4*alpha./Q);
  [K,E] = ellipke(k*k);

  Bz = 2*(Io/b)*(1./sqrt(Q)).*( K+E.*(1-alpha.^2-beta.^2)./(Q-4*alpha) );
  Br = 2*(Io/b)*(gamma./sqrt(Q)).*( -K+E.*(1+alpha.^2+beta.^2)./(Q-4*alpha) );
end