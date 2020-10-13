function [Bo,L,x] = biot_savart(xo, x, I)
%BIOT_SAVART Evaluate Biot-Savart integral
%
%   Bo = biot_savart(xo, x, I) computes, for each row i in xo,
%
%   Bo(i,:) = I * sum_k^N ( L(k,:) x (xo-x(k,:)) ) / |xo - x(k,:)|^3
%
%   where L = diff(x) and N = size(x,1).

L = diff(x);
x = x(1:end-1,:) + L/2; % Centers of segments.

Bo = zeros(size(xo));  
for i = 1:size(xo,1)
    Xo      = repmat(xo(i,:),size(L,1),1);
    X       = Xo - x;
    Rd      = (sum( X.^2 , 2)).^1.5;
    Bo(i,1) = I*sum( ( X(:,3).*L(:,2)-X(:,2).*L(:,3) )./Rd );
    Bo(i,2) = I*sum( ( X(:,1).*L(:,3)-X(:,3).*L(:,1) )./Rd );
    Bo(i,3) = I*sum( ( X(:,2).*L(:,1)-X(:,1).*L(:,2) )./Rd );
end