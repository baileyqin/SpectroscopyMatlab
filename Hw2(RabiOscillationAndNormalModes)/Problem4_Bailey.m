clear
X = load('CisDichloroethyleneForceMatrix.txt');
M=zeros(18);

au=1.66053892*10^(-27); %#kg in 1amu
C=12*au;
H=1.008*au;
Cl=35.4527*au;

for j=1:3
    M(j,j)= 1/sqrt(C);
end

for j=4:6
    M(j,j)= 1/sqrt(C);
end

for j=7:9
    M(j,j)= 1/sqrt(Cl);
end

for j=10:12
    M(j,j)= 1/sqrt(H);
end

for j=13:15
    M(j,j)= 1/sqrt(Cl);
end

for j=16:18
    M(j,j)= 1/sqrt(H);
end

bohr = 5.29e-11; % in m 
eV = 1.602176e-19; % in J 
hartree = 27.21138505*eV; % in J 

K=X*hartree/bohr/bohr;

H=M*K*M;

[V,D]=eig(H);

V

Z=eig(H);

Z

w=sqrt(Z);

c=2.99*10^10; % speed of light in cm/s

% lambda=c./w;

k=w/(c*2*pi);

k


geomEq = load('CisDichloroethyleneGeometry.txt'); 

close all 

for j = 1:18 
modeindex = j; 
geomt = zeros(size(geomEq)); 
displacement = 3e-14*M*V(:,modeindex); 
figure(j) 
plotlim = [-3,3]; 
xlim(plotlim); ylim(plotlim); zlim(plotlim); 

plotgeom(geomEq, j, 'g') 
plotgeom(geomEq + displacement,j,'r'); 
plotgeom(geomEq - displacement,j,'r'); 

end

