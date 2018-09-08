clear;
A=[cos(2*pi/3) -sin(2*pi/3) 0; sin(2*pi/3) cos(2*pi/3) 0; 0 0 1];
x=[-126.3139598;0;-190.73011914];
Z=[0 0 224 ; 0 0 137; 0 0 0 ;0 0 -146; x'; (A*x)'; (A*A*x)'].*10^-12;
M=[1.008;12.01;12.01;12.01;18.998403 ;18.998403 ;18.998403 ];
M=M*1.66053892*10^-27;
Ixx=0;
for i=1:7
    c=(Z(i,2)^2+Z(i,3)^2)*M(i);
    Ixx=Ixx+c;
end
Iyy=0;
for i=1:7
    c=(Z(i,1)^2+Z(i,3)^2)*M(i);
    Iyy=Iyy+c;
end
Izz=0;
for i=1:7
    c=(Z(i,1)^2+Z(i,2)^2)*M(i);
    Izz=Izz+c;
end
Ixy=0; %stuff is wrong
for i=1:7
    c=-Z(i,1)*Z(i,2)*M(i);
    Ixy=Ixy+c;
end

Ixz=0; %stuff is wrong
for i=1:7
    c=-Z(i,1)*Z(i,3)*M(i);
    Ixz=Ixz+c;
end

Iyz=0; %stuff is wrong
for i=1:7
    c=-Z(i,2)*Z(i,3)*M(i);
    Iyz=Iyz+c;
end

I=[Ixx,Ixy,Ixz;Ixy,Iyy,Iyz;Ixz,Iyz,Izz];

Id=diag(I)

hbar=1.05457173*10^-34;

figure(3)
hold on

for J=0:35
    for K=-J:J % the corresponding range of K 
        Energy = ((J*(J+1)*hbar^2)/(2*Id(1))+(hbar^2)*(K^2)*(1/(2*Id(3))-1/(2*Id(1)))); % Your expression for energy as a function of J and K 
        Energy = Energy/(hbar*2*pi)/10^9; % to convert to GHz 
        plot(J,Energy,'o','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k')
    end
end