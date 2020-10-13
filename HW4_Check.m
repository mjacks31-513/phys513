%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wire loop of radius 1 in the x-y plane with center at the origin.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b  = 1; % Wire radius
Io = 1; % Current in wire
zo = 1; % Observation point

xo = [0 0 zo];  

fprintf('Circular current loop; on-axis observation point.\n');

for k = 1:7
    Ns = 2*2^k;

    xyz      = zeros(Ns+1,3);
    th       = 2*pi*([0:Ns]')/Ns;
    xyz(:,1) = cos(th);
    xyz(:,2) = sin(th);
    xyz(:,3) = 0;

    [Bd,L,xc] = HW4_BiotSavart(xo,xyz,Io);	% Discrete solution

    Bc = HW4_B_ring(xo(3),b,Io);     % Continuous solution

    E(k) = 1-Bd(3)/Bc;

    fprintf(['Bz percent error from exact solution, '...
         'Ns = %3d: %8.3f%%\n'],Ns,100*E(k));

    if k == 4
        figure(1);clf;
        plot(xyz(:,1),xyz(:,2),'k-','Marker','o');
        hold on;grid on;
        plot(xc(:,1),xc(:,2),'k.','MarkerSize',20);
        legend('Segments and end points','Computed centers','Location','Best');
        axis square
    end

end


fprintf('Circular current loop; off-axis observation point.\n');

xo = [0 1 zo];  
r = sqrt(xo(1)^2 + xo(2)^2);

for k = 1:7
    Ns = 2*2^k;

    xyz      = zeros(Ns+1,3);
    th       = 2*pi*([0:Ns]')/Ns;
    xyz(:,1) = cos(th);
    xyz(:,2) = sin(th);
    xyz(:,3) = 0;

    [Bd,L,xc] = HW4_BiotSavart(xo,xyz,Io);	% Discrete solution

    [Bzc,Brc] = HW4_B_ring([xo(3),r],b,Io);     % Continuous solution

    Bmc = sqrt(Bzc^2 + Brc^2);
    Bmd = sqrt(sum(Bd.^2));
    
    E(k) = 100*(1-Bmd/Bmc);

    fprintf(['|B| percent error from exact solution, '...
         'Ns = %3d: %8.3f%%\n'],Ns,100*E(k));

end