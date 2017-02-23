clear all
close all

% Get ready to generate movie file
vidObj = VideoWriter('ls_swiral.avi');
open (vidObj);

N = 150;            % number of grid points in one direction
R1=.15;             % initial radius of circles 
h = 1/(N-1);        % grid spacing
dt = 0.2*h;         % time step
tfin = 2.6;         % total simulation time
nit = tfin/dt;      % number of time steps
vn = 0.0;           % normal speed this must be positive
x = -0:h:1;
y = x;
[X,Y] = meshgrid(x);

% Create a velocity field
u = -(sin(pi*X).^2).*cos(pi*Y).*sin(pi*Y);
v =  (sin(pi*Y).^2).*cos(pi*X).*sin(pi*X);

%   Initialize the level set function
a1 = 0.5;
b1 = 0.25;
phi = ((X-a1).*(X-a1)+(Y-b1).*(Y-b1)).^.5-R1;


%   arrays for the von Neumann boundary conditions
for i = 1:N
    ip(i) = i + 1;
    im(i) = i - 1;
end
im(1) = 1;
ip(N) = N;

%   begin simulation loop!
for iter = 1:nit
	for i = 1:N
        for j = 1:N
            dmx = (phi(i, j) - phi(im(i), j))/h;                     % x backward difference
            dpx = (phi(ip(i), j) - phi(i, j))/h;                     % x forward difference
            dmy = (phi(i, j) - phi(i, im(j)))/h;                     % y backward difference
            dpy = (phi(i, ip(j)) - phi(i, j))/h;                     % y forward difference
            dcx = (phi(ip(i), j) - phi(im(i), j))/2*h;
            dcy = (phi(i, ip(j)) - phi(i, im(j)))/2*h;
            fluxx = max(abs(max(dmx, 0)), abs(min(dpx, 0)));        % Godunov Flux x direction
            fluxy = max(abs(max(dmy, 0)), abs(min(dpy, 0)));        % Godunov Flux y direction
            phin(i,j) = phi(i,j) - vn*(fluxx^2 + fluxy^2)^(1/2)*dt - v(i,j)*((dpx + dmx)/2)*dt - u(i,j)*((dpy + dmy)/2)*dt;         % advance by dt            
        end
    end 
    phi = phin;
    
    %  Reinitialization
    for t = 1:5
        for i = 1:N
            for j = 1:N
                dmx = (phi(i, j) - phi(im(i), j))/h;                     % x backward difference
                dpx = (phi(ip(i), j) - phi(i, j))/h;                     % x forward difference
                dmy = (phi(i, j) - phi(i, im(j)))/h;                     % y backward difference
                dpy = (phi(i, ip(j)) - phi(i, j))/h;                     % y forward difference
                if (phin(i, j) > 0)
                    fluxx = max(abs(max(dmx, 0)), abs(min(dpx, 0)));        % Godunov Flux x direction
                    fluxy = max(abs(max(dmy, 0)), abs(min(dpy, 0)));        % Godunov Flux y direction
                else
                    fluxx = max(abs(max(dpx, 0)), abs(min(dmx, 0)));        % Godunov Flux x direction
                    fluxy = max(abs(max(dpy, 0)), abs(min(dmy, 0)));        % Godunov Flux y direction
                end
                sgn = phin(i, j)/(phin(i, j)^2 + h^2)^(1/2);
                grp = (fluxx^2 + fluxy^2)^(1/2);         
                phir(i, j) = phi(i, j) - (grp - 1)*sgn*dt;                   % advance by dt
            end
        end 
        phi = phir;   
    end
    phi = phir;

    %   Plotting
    contour(X,Y,phi,[0,0],'r');
    axis([0 1 0 1])
    axis('square')
    hold on
    streamslice(X,Y,u,v);
    hold off
    pause(.001)
    
    %   Calculating the area
    A(iter) = sum(trapz((1-heaviside(phi))))*h*h;
    F = getframe;
    writeVideo(vidObj, F);
end
close(vidObj);