clear

c = 3*10^8;  %in m/s
lambda = 530; %in nm
lambda = lambda * 10^-9 %m
w0 = c/lambda;  %in inverse seconds
k = 2.1*w0^2; %the constant Ne^2/epsilon*m
y = 5*10^12; %in inverse seconds
x=9 %chi because 1+chi=einfinity =10

syms w;

real(w) = 1 +x+ k*(w0^2-w^2)/((w0^2-w^2)^2 + w^2*y^2);
imaginary(w) = k*w*y/((w0^2-w^2)^2+w^2*y^2);

ni(w) = sqrt(-real + sqrt(real^2 + imaginary^2)) / sqrt(2);
a(w) = ni(w)/lambda;

h = ezplot(a(w),[10^10, 10^16]);

hold on
axis([0 1*10^15 0 10^7.4])

figure
h2 = ezplot(ni(w),[10^10, 10^16]);
axis([0 1*10^15 0 10])