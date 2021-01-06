function val = farfarpsi2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3) 

psi0_gauss = @(s1,s2,s3, alpha1,alpha2,alpha3) ((1/(2*pi.*alpha1)).^(1/4)).*exp(-((s1+0).^2)/(4.*alpha1)) ...
            .* ((1/(2*pi.*alpha2)).^(1/4))*exp(-((s2+1.5).^2)/(4.*alpha2)) ...
            .* ((1/(2*pi.*alpha3)).^(1/4)).*exp(-((s3-.5).^2)/(4.*alpha3));
%recenter gaussians
psi0_ppp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*sin(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_ppm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*sin(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_pmp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*cos(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_pmm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) sin(theta1).*cos(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);

if t1==0 && t2 == 0 && t3 ==0
     val = [psi0_ppp(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3),psi0_ppm(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3), ...
            psi0_pmp(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3),psi0_pmm(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3)];
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
transp_ppp = psi0_ppp(s1+t1,s2+t2,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_ppm = psi0_ppm(s1+t1,s2+t2,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_pmp = psi0_pmp(s1+t1,s2-t2,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_pmm = psi0_pmm(s1+t1,s2-t2,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);

%Integral 1
int1_ppp = -1*1i*omega/2*integral(@(sigma) psi0_pmp(s1+t1,(sigma),s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
int1_ppm = -1*1i*omega/2*integral(@(sigma) psi0_pmm(s1+t1,(sigma),s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
int1_pmp = -1*1i*omega/2*integral(@(sigma) psi0_ppp(s1+t1,(sigma),s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
int1_pmm = -1*1i*omega/2*integral(@(sigma) psi0_ppm(s1+t1,(sigma),s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);

%Integral 2
int2_ppp = -1*1i*omega/2*integral(@(tau) psi0_ppm(s1+t1,s2+t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
int2_ppm = -1*1i*omega/2*integral(@(tau) psi0_ppp(s1+t1,s2+t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
int2_pmp = -1*1i*omega/2*integral(@(tau) psi0_pmm(s1+t1,s2-t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
int2_pmm = -1*1i*omega/2*integral(@(tau) psi0_pmp(s1+t1,s2-t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);

%Final values
val_ppp = transp_ppp + int1_ppp + int2_ppp;
val_ppm = transp_ppm + int1_ppm + int2_ppm;
val_pmp = transp_pmp + int1_pmp + int2_pmp;
val_pmm = transp_pmm + int1_pmm + int2_pmm;

%Output
val = [val_ppp,val_ppm,val_pmp,val_pmm];
end