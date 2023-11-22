% Title: Vorticity-Stream Function at Square Corner
% Author: Jonathan Blisko
% Date: 05/13/20

clear all
close all
clc

%% Initialize Variables
% X-Direction
L=1;
deltax=0.01;
nx=L/deltax+1;
x=0:deltax:L;

% Y-Direction
H=L;
deltay=deltax;
ny=H/deltay+1;
y=0:deltay:H;

% Constants
U=1;
p=(cos(pi/(nx-1))+(deltax/deltay)^2*cos(pi/(ny-1)))/(1+(deltax/deltay)^2);
w_opt=2/(1+sqrt(1-p^2));
hf=floor((nx-1)/2);
Re=200;
deltat=deltax^2/(4*U);
c=deltat/(2*deltax);
r=deltat/deltax^2;
SORConvergence=1e-8;
SteadyConvergence=1e-6;
k=1;
max_it=50000;
diff=1;

%% Initial Conditions
w=zeros(nx,ny);
psi=zeros(nx,ny);
u=zeros(nx,ny);
v=zeros(nx,ny);

% Inlet and Oulet Stream Function
v(2:hf,ny)=U;
psi(2:hf+1,ny)=v(2:hf+1,ny).*linspace(L/2,0,hf)';
%psi(2:hf+1,ny)=deltax*v(2:hf+1,ny)+psi(2:hf+1,ny-1);

% Exclude Upper Right Quadrant
w(hf+2:nx,hf+2:ny)=NaN;
psi(hf+2:nx,hf+2:ny)=NaN;
u(hf+2:nx,hf+2:ny)=NaN;
v(hf+2:nx,hf+2:ny)=NaN;

%% Solving Equations

while (diff>=SteadyConvergence && k<max_it)
    w_old=w;
    
    % Transport Equation: Left Column
    for ii=2:hf
        for j=2:ny-1
            w(ii,j)=w(ii,j)-c*u(ii,j)*(w(ii+1,j)-w(ii-1,j))-c*v(ii,j)*(w(ii,j+1)-w(ii,j-1))+(r/Re)*(w(ii+1,j)+w(ii-1,j)+w(ii,j+1)+w(ii,j-1)-4*w(ii,j));
        end
    end
    
    % Transport Equation: Right Row
    for ii=hf+1:nx-1
        for j=2:hf
            w(ii,j)=w(ii,j)-c*u(ii,j)*(w(ii+1,j)-w(ii-1,j))-c*v(ii,j)*(w(ii,j+1)-w(ii,j-1))+(r/Re)*(w(ii+1,j)+w(ii-1,j)+w(ii,j+1)+w(ii,j-1)-4*w(ii,j));
        end
    end
    
    SOR_error=1;
    while (SOR_error>=SORConvergence)
        psiold=psi;
        % SOR Method: Left Column
        for j=ny-1:-1:2
            for ii=2:hf
                psi(ii,j)=(w_opt/4)*(psi(ii+1,j)+psi(ii-1,j)+psi(ii,j+1)+psi(ii,j-1)+deltax^2*w(ii,j))+(1-w_opt)*psi(ii,j);
            end
        end
    
        % SOR Method: Right Row
        for ii=hf+1:nx-1
            for j=2:hf
                psi(ii,j)=(w_opt/4)*(psi(ii+1,j)+psi(ii-1,j)+psi(ii,j+1)+psi(ii,j-1)+deltax^2*w(ii,j))+(1-w_opt)*psi(ii,j);
            end
        end
        
        psi(nx,2:ny-1)=psi(nx-1,2:ny-1);
        
        SOR_error=max(max(abs(psi-psiold)));
    end
    
    % Velocity Update: Left Column
    for ii=2:hf
        for j=2:ny-1
            u(ii,j)=(psi(ii,j+1)-psi(ii,j-1))/(2*deltax);
            v(ii,j)=-(psi(ii+1,j)-psi(ii-1,j))/(2*deltax);
        end
    end
    
    % Velocity Update: Right Row
    for ii=hf+1:nx-1
        for j=2:hf
            u(ii,j)=(psi(ii,j+1)-psi(ii,j-1))/(2*deltax);
            v(ii,j)=-(psi(ii+1,j)-psi(ii-1,j))/(2*deltax);
        end
    end
    
    u(nx,2:ny-1)=u(nx-1,2:ny-1);    % Outlet Velocity
    v(nx,2:ny-1)=v(nx-1,2:ny-1);
    
    % Vorticity BC's: Left Column
    for ii=2:hf
        w(ii,1)=2*(psi(ii,1)-psi(ii,2))/(deltax^2);    % Bottom
        w(ii,ny)=2*(psi(ii,ny)-psi(ii,ny-1))/(deltax^2);    % Inlet
        w(1,ii)=2*(psi(1,ii)-psi(2,ii))/(deltax^2);    % L.H.S.
    end
    
    % Vorticity BC's: Right Row
    for ii=hf+1:nx-1
        w(ii,1)=2*(psi(ii,1)-psi(ii,2))/(deltax^2);    %Bottom
        w(ii,hf+1)=2*(psi(ii,hf+1)-psi(ii,hf))/(deltax^2);    % Top R.H.S.
        w(1,ii)=2*(psi(1,ii)-psi(2,ii))/(deltax^2);    % L.H.S.
        w(hf+1,ii)=2*(psi(hf+1)-psi(hf,ii))/(deltax^2);    % R.H.S.
    end
    w(nx,:)=w(nx-1,:);  % Outlet
    
    diff=max(max(abs(w-w_old)));
    
    k=k+1;
end

figure
[X,Y]=meshgrid(x,y);
vv=[-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1,0,1,5e-6,1e-4,1e-3,1e-2,1e-1,-0.8:0.1:-0.2,-.04,.005,.02,-.25,-.15,-.075,-.02];
%vv=[0:0.05:0.5,-.1:.01:.1];
contour(Y,X,psi,vv);
axis equal
colorbar
colormap('parula')
hold on
t = (1/8:1/4:1)'*2*pi;
x = cos(t)+sqrt(2)/2+hf*deltax;
y = sin(t)+sqrt(2)/2+hf*deltax;
fill(y,x,[220 220 220]/255)
xlim([0 1])
ylim([0,1])
xlabel('X-Direction','Interpreter','Latex')
ylabel('Y-Direction','Interpreter','Latex')
title('Stream Function Contour Plot','Interpreter','Latex','Fontsize',14)
set(gcf,'RendererMode','manual')

figure
ww=[-250:50:-100,-10:2:10,10:5:30,-1e-2,1e-2,1e-4,-100:10:-20,-32.4,-24.5];
contour(Y,X,w,ww)
axis equal
colorbar
colormap('parula')
hold on
fill(y,x,[220 220 220]/255)
xlim([0 1])
ylim([0,1])
xlabel('X-Direction','Interpreter','Latex')
ylabel('Y-Direction','Interpreter','Latex')
title('Vorticity Contour Plot','Interpreter','Latex','Fontsize',14)
set(gcf,'RendererMode','manual')

figure
zz=[1e-5,1e-4,1e-3,1e-2,1e-1,1:1:5,-3,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1,-5e-3,1.8];
contour(Y,X,v,zz)
axis equal
colorbar
colormap('parula')
hold on
fill(y,x,[220 220 220]/255)
xlim([0 1])
ylim([0,1])
xlabel('X-Direction','Interpreter','Latex')
ylabel('Y-Direction','Interpreter','Latex')
title('Vertical Velocity Contour Plot','Interpreter','Latex','Fontsize',14)
set(gcf,'RendererMode','manual')

figure
kk=[1e-5,1e-4,1e-3,1e-2,1e-1,1,5,10,-1e-5,-1e-4,-1e-3,-1e-2,-1e-1,-1:-1:-4,0.05];
contour(Y,X,u,kk)
axis equal
colorbar
colormap('parula')
hold on
fill(y,x,[220 220 220]/255)
xlim([0 1])
ylim([0,1])
xlabel('X-Direction','Interpreter','Latex')
ylabel('Y-Direction','Interpreter','Latex')
title('Horizontal Velocity Contour Plot','Interpreter','Latex','Fontsize',14)
set(gcf,'RendererMode','manual')
