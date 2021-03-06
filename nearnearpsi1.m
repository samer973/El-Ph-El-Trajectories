function val = nearnearpsi1(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3,phi) 
% define some functions that are needed in the computation
% radical = @(s, t, sigma) sqrt(t.^2 - (s - sigma).^2);
% radicalp = @(s, t, sigma) (sqrt(t + (s - sigma)))/(sqrt(t - (s - sigma)));
% mybessel1 = @(z,s,t,omega) besselj(1,omega*radical(s, t, z));
% mybessel0= @(z,s,t,omega) besselj(0,omega*radical(s, t, z));
% specify the initial data for KG eq, coming from the initial data for
% Dirac's equation for the electron with mass omega that corresponds to the initial
% probability density function being normal distribution with mean zero and
% variance alpha, and the mixing angle theta

% add e^(ik\cdot x) factor for photon
psi0_gauss = @(s1,s2,s3, alpha1,alpha2,alpha3) ((1/(2*pi.*alpha1)).^(1/4)).*exp(-((s1+0).^2)/(4.*alpha1)) ...
            .* ((1/(2*pi.*alpha2)).^(1/4))*exp(-((s2+1.5).^2)/(4.*alpha2))  ...
            .* ((1/(2*pi.*alpha3)).^(1/4)).*exp(-((s3-.5).^2)/(4.*alpha3));

psi0_ppp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*sin(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_ppm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*sin(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_pmp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*cos(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_pmm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*cos(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mpp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*sin(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mpm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*sin(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mmp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*cos(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mmm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*cos(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);

if t1==0 && t2 == 0 && t3 ==0
     val = [psi0_mpp(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3),psi0_mpm(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3), ...
            psi0_mmp(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3),psi0_mmm(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3)];
return
end
% %compute the second integral in the KG solution
% y_2 = @(z,s,t,theta,alpha,omega) mybessel0(z,s,t,omega).*g(z,theta, alpha);
% q = 1i*(0.5)*omega*integral(@(z) y_2(z,s,t,theta,alpha,omega),s-t,s+t);
% % compute the first integral in the KG solution
% y_1 = @(z,s,t,theta,alpha,omega) mybessel1(z,s,t,omega).*f(z, theta, alpha).*radicalp(s, t, z);
% w = omega*0.5*integral(@(z) y_1(z,s,t,theta,alpha,omega),s-t,s+t);
% %Almost there! Now just need first part of KGE solution
% o = f(s-t, theta, alpha);
% %FINALLY put all three pieces together
% val = o-w-q;

%Transport:
transp_mpp = exp(1i .* phi)*psi0_pmp(s2+t2,s1-t1,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_mpm = exp(1i .* phi)*psi0_pmm(s2+t2,s1-t1,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);

transp_mmp = psi0_mmp(s1-t1,s2-t2,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_mmm = psi0_mmm(s1-t1,s2-t2,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);

%Integral 1
int1_mpp = -1*1i*omega/2*integral(@(tau) exp(1i .* phi).*psi0_pmm(s2+t2,s1-t1,tau,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
int1_mpm = -1*1i*omega/2*integral(@(tau) exp(1i .* phi).*psi0_pmp(s2+t2,s1-t1,tau,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);

int1_mmp = -1*1i*omega/2*integral(@(sigma) psi0_mpp(s1-t1,(sigma),s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s1-t1);
int1_mmm = -1*1i*omega/2*integral(@(sigma) psi0_mpm(s1-t1,(sigma),s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s1-t1);

%Integral 2
int2_mpp = -1*1i*omega/2*integral(@(sigma) exp(1i * phi)*psi0_ppp(s2+t2,sigma,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s1-t1,s2+t2);
int2_mpm = -1*1i*omega/2*integral(@(sigma) exp(1i * phi)*psi0_ppm(s2+t2,sigma,t3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s1-t1,s2+t2);

int2_mmp = -1*1i*omega*integral(@(sigma) exp(1i * phi)*psi0_pmp(s2-t2+2*sigma,s1-t1,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),(1/2)*(t2-s2+s1-t1),t2);
int2_mmm = -1*1i*omega*integral(@(sigma) exp(1i * phi)*psi0_pmm(s2-t2+2*sigma,s1-t1,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),(1/2)*(t2-s2+s1-t1),t2);

%Int 3
int3_mpp = -1*1i*omega/2*integral(@(sigma) psi0_mmp(s1-t1,sigma,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s1-t1);
int3_mpm = -1*1i*omega/2*integral(@(sigma) psi0_mmm(s1-t1,sigma,t3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s1-t1);

int3_mmp = -1*1i*omega/2*integral(@(tau) psi0_mmm(s1-t1,s2-t2,tau,theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
int3_mmm = -1*1i*omega/2*integral(@(tau) psi0_mmp(s1-t1,s2-t2,tau,theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);

%Final values
val_mpp = transp_mpp + int1_mpp + int2_mpp + int3_mpp;
val_mpm = transp_mpm + int1_mpm + int2_mpm + int3_mpm;

val_mmp = transp_mmp + int1_mmp + int2_mmp + int3_mmp;
val_mmm = transp_mmm + int1_mmm + int2_mmm + int3_mmm;

%Output
val = [val_mpp,val_mpm,val_mmp,val_mmm];
end