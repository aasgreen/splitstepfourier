clear all;
distance = input('Enter fiber lengh (in units of L_D)');
beta2 = input('dispersion:1 for normal, -1 for anomalous');
N = input('nonlinear parameter N=');
mshape = input('m=0 for sech, m> 0 for super gaussian=');
chirp0=0;

nt = 1024; Tmax = 32;
step_num = round(20*distance*N^2);
deltaz = distance/step_num;
dtau = (2*Tmax)/nt;

tau = (-nt/2:nt/2-1)*dtau;
omega = (pi/Tmax)*[0:nt/2-1,-nt/2:-1];

if mshape ==0
    uu = sech(tau).*exp(-0.5i*chirp0*tau.^2);
else
    uu = exp(-0.5*(1+1i*chirp0).*tau.^(2*mshape));
end

temp = fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi);
figure; subplot(2,1,1)
    plot(tau, abs(uu).^2,'--k'); hold on
    axis([-5 5 0 inf]);
    subplot(2,1,2)
    plot(fftshift(omega)/(2*pi), abs(temp).^2, '--k'); hold on;
    axis([-.5 .5 0 inf]);

dispersion = exp(i*0.5*beta2*omega.^2*deltaz);
hhz = 1i*N^2*deltaz;

temp = uu.*exp(abs(uu).^2.*hhz/2);
for n=1:step_num
    f_temp = ifft(temp).*dispersion;
    uu = fft(f_temp);
    temp = uu.*exp(abs(uu).^2.*hhz);
end

uu = temp.*exp(-abs(uu).^2.*hhz/2.)
temp =fftshift(ifft(uu)).*(nt*dtau)/sqrt(2*pi)


subplot(2,1,1)
plot(tau,abs(uu).*2, '-k')
subplot(2,1,2)
plot(fftshift(omega/(2*pi)), abs(temp).^2,'-k')
