clear all

% Number1
% part b
hbar=1.055*10^-34;% in Js
mu=1.6605*10^-27;% kg
angstrom=10^10; %angstrom to m
cmtoJ=1.986*10^-23; %inverse cm to J

dx=0.005;
R=0:dx:5;

w=214.5*cmtoJ/hbar; %in 1/s
m=126.904*mu/2; %in kg
De=12440;% in inverse cm
R0=2.666;% in angstrom
B = sqrt(m*w^2/(2*De*cmtoJ))/angstrom; % in inverse angstrom

figure(1)

PESg=De*(1-exp(-B*(R-R0))).^2; %in inverse cm
plot(R,PESg,'red')
axis([2 5 0 20*10^3])
hold on

w2=125.65*cmtoJ/hbar; % in inverse seconds
De2=4349;%in inverse cm
R02=3.024;%in angstroms
B2=sqrt(m*w2^2/(2*De2*cmtoJ))/angstrom; %in inverse angstrom

PESe=De2*(1-exp(-B2*(R-R02))).^2+15769;%in inverse cm
plot(R,PESe,'red')
axis([2 5 0 20000])
hold on

%number 2
%part b
DeltaX = dx/angstrom; %in m
alpha = hbar^2/(2*m*DeltaX^2); %Kinetic term in J
alpha = alpha/cmtoJ; % kinetic term in inverse cm

H = zeros(length(R));
for i=1:length(R)
    H(i,i) = 0.5*(PESg(i)+2*alpha);
end
for i=1:(length(R)-1)
    H(i,i+1)=-alpha;
end

H=H+transpose(H); %hamiltonian for the ground state in inverse cm

[V, E] = eig(H); % wavefunctions and energy levels of ground state

for i=1:100
    plot(R,V(:,i)*500 + E(i,i),'blue')
end

H2 = zeros(length(R));
for i=1:length(R)
    H2(i,i) = 0.5*(PESe(i) + 2*alpha);
end

for i=1:(length(R)-1)
    H2(i,i+1)=-alpha;
end

H2=H2+transpose(H2); %hamiltonian for the excited state in inverse cm

[V2, E2] = eig(H2); %wavefunctions and energy levels of excited state

for i=1:100
    plot(R,V2(:,i)*500 + E2(i,i),'blue')
end

%Number 3
%part a

for i=1:length(R)
    deltaE(i)=E2(i,i)-E(1,1);
    FCF(i)=dot(V(:,1),V2(:,i));
end
FCF=FCF.^2;
figure(2)
plot(deltaE, FCF)
hold on;
axis([15000 25000 0 0.2])

%part b

kb=1.38*10^-23;%boltzman constant in J/K

for j=1:4
for i=1:length(R)
    deltaE(i)=1/(E2(i,i)-E(j,j))*10^7;
    FCF(i)=dot(V(:,j)*exp(-E(j,j)*cmtoJ/kb/300),V2(:,i));
end
FCF=FCF.^2;
figure(j+2)
bar(deltaE, FCF);
hold on;
axis([500 700 0 .1])
end

%part c

for i=1:length(R)
    deltaE(i)=-1/(E(i,i)-E2(i,i))*10^7;
    FCF(i)=dot(V(:,i),V2(:,1));
end
FCF=FCF.^2;
figure(7)
bar(deltaE, FCF)
hold on;
axis([600 750 0 1])

%number 4

De1 = 7000; %in inverse cm
R1 = 2.9; %in angstroms
wex=125.69*cmtoJ/hbar; %in inverse s
Bex=sqrt((wex^2*m)/(2*De1*cmtoJ))/angstrom; %in inverse angstrom

PESg = 15769+De1*(1-exp(-Bex*(R-R1))).^2;
H3=zeros(length(R));

for i=1:length(R)
        H3(i,i) = 0.5*(PESg(i) + 2*alpha);
end

for i=1:(length(R)-1)
        H3(i,i+1) = -alpha;
end

H3=H3+transpose(H3);

[V3,E3]=eig(H3);

Psi = zeros(length(E3),1);

figure(8)
for t =0:2000
    for i=1:length(R)
        FCF(i)=dot(V(:,1),V3(:,i));
        Psi(:) = Psi(:)+FCF(i)*exp(-1i/hbar*E3(i,i)*cmtoJ*(t)*10^(-15))*(V3(:,i));
    end
    plot(R, abs(Psi*500+E3(:,1))+E(1,1), R,(PESg-(15769+De1*(1-exp(-Bex*(R0-R1))).^2))+100)
    axis([2 5 100 300])
    Psi = zeros(length(E3),1);
    title(num2str(t));
    drawnow
end

