clear

K=zeros(9);
K(1,1)=1600;
K(1,4)=-1600;
K(4,4)=3200;
K(4,1)=-1600;
K(4,7)=-1600;
K(7,7)=1600;
K(7,4)=-1600;

amu=1.66053892*10^(-27); %#kg in 1amu
C=12.0107*amu;
O=15.9994*amu;

M=zeros(9);
M(1,1)=1/sqrt(O);
M(2,2)=1/sqrt(O);
M(3,3)=1/sqrt(O);
M(4,4)=1/sqrt(C);
M(5,5)=1/sqrt(C);
M(6,6)=1/sqrt(C);
M(7,7)=1/sqrt(O);
M(8,8)=1/sqrt(O);
M(9,9)=1/sqrt(O);

H=M*K*M;

[V,D]=eig(H);

w=sqrt(eig(H));

c=2.99*10^10; % speed of light in cm/s

% lambda=c./w;
% 
% k=2*pi./lambda

k=w/(c*2*pi)

