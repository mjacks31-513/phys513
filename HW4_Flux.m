
Nsl      = 256; % Number of loop segments.
xyz      = zeros(Nsl+1,3);
th       = 2*pi*([0:Nsl]')/Nsl;
xyz(:,1) = cos(th);
xyz(:,2) = sin(th);
xyz(:,3) = 0;

% Compute Bz at z = 2, r = [0, dr, 2*dr, ..., 1].
% Flux is approximated as Bz*2*pi*dr;

Io = 1;
b   = 1;

Ns = 4; % Number of annular area segments.
r = b*[1:Ns]'/Ns;
xo = zeros(size(r));

xo(:,1) = r;
xo(:,2) = 0;
xo(:,3) = 2;

for i = 1:Ns
    Bd(i,:) = HW4_BiotSavart(xo(i,:),xyz,Io);
end

dr = r(2)-r(1);
F = sum(Bd(:,3)*2*pi*dr);

fprintf('HW4.3 answer: Flux = %.2f * mu_o/(4*pi)\n', F);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check: When zo is large and number of segments large, 
% flux should approach Bzo*pi*r^2, where Bzo is on-axis field of ring.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns = 256; % Number of annular area segments.
r = b*[1:Ns]'/Ns;
xo = zeros(size(r));

zo = 1000;

xo(:,1) = r;
xo(:,2) = 0;
xo(:,3) = zo;

for i = 1:Ns
    Bd(i,:) = HW4_BiotSavart(xo(i,:),xyz,Io);
end

dr = r(2)-r(1);
F1 = sum(Bd(:,3)*2*pi.*r*dr);

Bc = HW4_B_ring(zo,b,Io); 

F2 = Bc*pi*b^2;

fprintf('%% difference between computed and limiting value = %.2f\n',100*(1-F1/F2));
