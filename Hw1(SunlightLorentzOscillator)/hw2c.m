
lambda=533; %wavelength in nm
c=300; % speed of light in nm/fs
k=2*pi/lambda;
w=k*c;
r=1
t=0:.01:15; %time that goes from 0 to 15 fs
x=cos(w*t-k*r);
y=4*cos(w*t-k*r+pi/4);
plot3(t,x,y)