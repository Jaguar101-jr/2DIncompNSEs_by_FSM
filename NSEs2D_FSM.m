%% A MATLAB PROGRAM TO SOLVE THE 2D Incompressible NSEs with Vortocity by FOURIER-SPECTRAL METHODS.
%% A Matlab code written by and developed by HUSSEIN A. H. Muhammed Oct. 2022.
%% B.Sc.H AND M.Sc. (Honuors).
%% 

%Simulation Property Setting
GridSize=128;
Visc= 0.005;

% Prepare the movie file.

    vidObj = VideoWriter('2DIncompNSEs_by_FSM.mp4');
    open(vidObj);

%Space Setting
h=2*pi/GridSize;
axis=h*[1:1:GridSize];
[x,y]=meshgrid(axis, axis);

%Time Setting
FinTime=70;
dt=0.1;
t=0;

%Movie File Data Allocation Setup
FrameRate=10;
Mov(10)=struct('cdata',[],'colormap',[]);

k=0;
j=1;

%Defining Initial Vorticity Distribution
H=exp (-((x-pi+pi/5).^2 + (y-pi+pi/5).^2)/(0.3)) - exp(-((x-pi-pi/5).^2 ...
+(y-pi+pi/5).^2)/(0.2)) + exp (-((x-pi-pi/5).^2 + (y-pi-pi/5).^2)/(0.4));

%Adding Random Noise to Initial Vorticity
epsilon = 0.3;
Noise=random('unif',-1,1,GridSize, GridSize);

%Note that for Low Viscosities Adding Noise to Non?Trivial Vorticity
%Distribution results in blow up , so either do pure noise or smooth data

w = H + epsilon * Noise;
w_hat = fft2(w);

%%%%%%%%% Method Begins Here %%%%%%%%%%

kx=1i*ones(1,GridSize)'*(mod((1:GridSize)-ceil(GridSize/2+1),GridSize)-floor(GridSize/2));
ky=1i*(mod((1:GridSize)'-ceil(GridSize/2+1) , GridSize) - floor(GridSize/2))*ones(1,GridSize);
AliasCor=kx<2/3*GridSize&ky<2/3*GridSize;

Lap_hat=kx.^2+ky.^2;

ksqr=Lap_hat; ksqr(1,1)=1;

while t<FinTime
    
psi_hat = -w_hat./ksqr;

u=real(ifft2(ky.*psi_hat));
v=real(ifft2(-kx.*psi_hat));
wx=real(ifft2(kx.*w_hat));
wy=real(ifft2(ky.*w_hat));
VgradW=u.*wx+v.*wy;
VgradW_hat=fft2(VgradW);
VgradW_hat=AliasCor.*VgradW_hat;

%Crank?Nicholson Update Method
w_hat_update=1./(1/dt-0.5*Visc*Lap_hat).*((1/dt+0.5*Visc*Lap_hat).*w_hat - VgradW_hat);
if (k==FrameRate)

w=real(ifft2(w_hat_update));
%Vel=sqrt(u.^2+v.^2); %This is for plotting velocity

contourf(x,y,w,80);
colorbar;

shading flat; colormap ('jet');
drawnow
Mov(j)=getframe;
k=0;
j=j+1
end
w_hat=w_hat_update;
t=t+dt;
k=k+1;
title(sprintf('NSEs in Fourier-Domain:propagation time(s) = %.2f' , t));
%Write each frame to the file
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
end
%Close the file.
close(vidObj);

