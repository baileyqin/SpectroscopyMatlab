clear

x = load('solarspectrumEquallySpacedinFreq.txt');

c=300; %in nm/fs
ev=1.602176*10^-19; %J/ev
epsilon=8.854*10^-12*ev*10^-9; %vacuum permittivity in C^2/eVnm
hbar=6.582119*10^-16*10^15;%in eVfs
w0=c/(500)

i=x(:,1)/(ev*10^-3)/ev*10^-15*(10^9)^2; %intensity in eV/nm^2
w=x(:,2);
w=w./hbar;%in 1/fs

t=0:.1:100; %in fseconds

for j=1:length(i)
    z(j)=rand(1)*2*pi;
end

a=c/(495) - c/(505);
g=zeros(length(i));
for j=1:length(i)
    g(j) = exp(-((w(j) - w0)^2)/(2*a^2))/(a*sqrt(2*pi));
end


for j=1:length(i)
    E(j)=sqrt((i(j)*2)/(c*epsilon))*g(j);
end

En=zeros(1,length(t));

for j=1:length(t)
    for k=1:length(i)
        En(j)=En(j)+E(k)*sin(w(k));
    end
end
plot(t,En)


%for j=1:length()