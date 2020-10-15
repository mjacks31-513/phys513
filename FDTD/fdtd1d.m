
% The following is a modified version of fdtd_original.m obtained from
% https://www.mathworks.com/matlabcentral/fileexchange/7459-fdtd1d-m
% The additions include the ability to extract the field at certain
% locations and to define different runs.
% To speed up the run, comment out the plotting in the loop.

% The basics of the algorithm when sigma = 0 is described in
% https://my.ece.utah.edu/~ece6340/LECTURES/lecture%2014/FDTD.pdf
% See also
% http://www.astrosen.unam.mx/~aceves/Fisica_Computacional/ebooks/sullivan_emsimulation_fdtd.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scott Hudson, WSU Tri-Cities
%1D electromagnetic finite-difference time-domain (FDTD) program.
%Assumes Ey and Hz field components propagating in the x direction.
%Fields, permittivity, permeability, and conductivity
%are functions of x. Try changing the value of "profile".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','docked');
%close all;
clear all;

L = 5;      % Domain length in meters
N = 500;    % # spatial samples in domain
x = linspace(0,L,N); % x coordinate of spatial samples

eps0 = 8.854e-12; %permittivity of free space
mu0 = pi*4e-7; %permeability of free space
animate = 0; % Live
run = 2;
if run == 0
    profile = 0; % eps = eps_o, mu = mu_o, sigma = 0.
    source = 2;  % Gaussian pulse at left boundary
    xg = L;
    Niter = 600; % # of iterations to perform
end
if run == 1
    profile = 5; % Step in conductivity
    source = 2;  % Gaussian pulse at left boundary
    xg = L/4;
    Niter = 300; % # of iterations to perform
end
if run == 2
    profile = 5; % Step in conductivity
    source = 1;  % sinusoid
    sigmac = 0.5;
    Niter = 700; % # of iterations to perform
    fs = 300e6; % Source frequency in Hz
    xg = L/2; % Position of change in conductivity
    xp = L/2; % Position of probe
    Ip = find(x >= xp,1)-1; % Index associated with xp
end
fprintf('\n')
Ig = find(x >= xg,1)-1;

ds = L/N; %spatial step in meters
dt = ds/300e6; % "magic time step"
showWKB=0; % if =1 then show WKB appromination at end

%scale factors for E and H
ae = ones(N,1)*dt/(ds*eps0);
am = ones(N,1)*dt/(ds*mu0);
as = ones(N,1);
epsr = ones(N,1);
mur= ones(N,1);
sigma = zeros(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here we specify the epsilon, sigma, and mu profiles. I've
%predefined some interesting examples. Try profile = 1,2,3,4,5,6 in sequence.
%You can define epsr(i), mur(i) (relative permittivity and permeability)
%and sigma(i) (conductivity) to be anything you want.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N
   epsr(i) = 1;
   mur(i) = 1;
   w1 = 0.5;
   w2 = 1.5;
   if (profile == 0)
       % eps = eps_o, mu = mu_o, sigma = 0.
   end
   if (profile==1) %dielectric window
      if (abs(x(i)-L/2)<0.5) epsr(i)=4; end
   end
   if (profile==2)%dielectric window with smooth transition
   	if (abs(x(i)-L/2)<1.5) epsr(i)=1+3*(1+cos(pi*(abs(x(i)-L/2)-w1)/(w2-w1)))/2; end
      if (abs(x(i)-L/2)<0.5) epsr(i)=4; end
   end
   if (profile==3) %dielectric discontinuity
      if (x(i)>L/2) epsr(i) = 9; end
   end
  	if (profile==4) %dielectric disontinuity with 1/4-wave matching layer
      if (x(i)>(L/2-0.1443)) epsr(i) = 3; end 
      if (x(i)>L/2) epsr(i) = 9; end
   end
   if (profile==5) %conducting half space
      if (x(i)>=xg) sigma(i) = sigmac; end
      %if (x(i)>3*L/4) sigma(i) = sigmac; end
   end
   if (profile==6) %sinusoidal dielectric
      epsr(i) = 1+sin(2*pi*x(i)/L)^2;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ae = ae./epsr;
am = am./mur;
ae = ae./(1+dt*(sigma./epsr)/(2*eps0));
as = (1-dt*(sigma./epsr)/(2*eps0))./(1+dt*(sigma./epsr)/(2*eps0));

%plot the permittivity, permeability, and conductivity profiles
figure(1)
    subplot(3,1,1);
    plot(x,epsr);
    grid on;
    axis([3*ds L min(epsr)*0.9 max(epsr)*1.1]);
    title('relative permittivity');
    subplot(3,1,2);
    plot(x,mur);
    grid on;
    axis([3*ds L min(mur)*0.9 max(mur)*1.1]);
    title('relative permeabiliity');
    subplot(3,1,3);
    plot(x,sigma);
    grid on;
    axis([3*ds L min(sigma)*0.9-0.001 max(sigma)*1.1+0.001]);
    title('conductivity');
    drawnow
%initialize fields to zero
Hz = zeros(N,1);
Ey = zeros(N,1);

figure(2);clf
    set(gcf,'doublebuffer','on'); % For smoother graphics
    grid on;
    plot(x,Ey,'b');
    hold on;
    axis([3*ds L -2 2]);
    plot(x,377*Hz,'r');
    %legend('E_y','377H_z')

for iter=1:Niter
    % Source
    if source == 1
        % "smooth turn on" sinusoidal source 
        Ey(100) = Ey(100)+2*(1-exp(-((iter-1)/50)^2))*sin(2*pi*fs*dt*iter);
        %Ey(3) = sin(2*pi*fs*dt*iter); 
        %Ey(251) = sin(2*pi*fs*dt*iter);
        %Ey(100) = sin(2*pi*fs*dt*iter);
    end
    if source == 2
        % Gaussian pulse
        Ey(3)=exp(-((iter-10)/2)^2);
    end
    
    % The next 10 or so lines of code are where we actually integrate Maxwell's
    % equations. All the rest of the program is basically bookkeeping and plotting.
    Hz(1) = Hz(2); % Absorbing boundary conditions for left-propagating waves
    for i=2:N-1 % Update H field
      Hz(i) = Hz(i)-am(i)*(Ey(i+1)-Ey(i));
    end
    Ey(N) = Ey(N-1); % Absorbing boundary conditions for right-propagating waves
    for i=2:N-1 % Update E field
      Ey(i) = as(i)*Ey(i)-ae(i)*(Hz(i)-Hz(i-1));
    end
    
    Hz_Ip(iter) = Hz(Ip);
    Ey_Ip(iter) = Ey(Ip);

    if (animate || iter == Niter)
        figure(2);hold off;
           plot(x,Ey,'b');
           hold on;
           axis([3*ds L -2 2]);
           plot(x,377*Hz,'r');
           title('E_y (blue), 377H_z (red)');
           %legend('E_y','377H_z') % Slows down rendering
           xlabel('x [m]');
           p = patch([xg 5.0 5.0 xg],[-2.0 -2.0 2.0 2.0],[-1,-1,-1,-1],[0.5,0.5,0.5]);
           set(p,'FaceAlpha',0.5,'EdgeColor','none');
           drawnow;
           %pause(0);
    end

    if iter > 1
        fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    end
    fprintf('Step = %04d',iter);

end
fprintf('\n');

figure(3);clf;hold on;box on;grid on;
    plot(Ey_Ip,'b');
    plot(377*Hz_Ip,'r');
    xlabel('step');
    legend('E_y','377H_z','Location','NorthWest')
    title(sprintf('At x = %.4f',x(Ip)));
    drawnow
    
%WKB prediction for the fields
if (showWKB==1)
	phase = cumsum((epsr).^0.5)*ds;
	beta0 = 2*pi*fs/(300e6);
   theory = sin(2*pi*fs*(Niter+4)*dt-beta0*phase)./(epsr.^0.25);
   input('press enter to show WKB theory');
   plot(x,theory,'b.');
   theory = sin(2*pi*fs*(Niter+4)*dt-beta0*phase).*(epsr.^0.25);
   plot(x,theory,'r.');
   title('E (blue), 377*H (red), WKB theory (points)');
end
