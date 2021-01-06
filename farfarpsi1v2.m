function val = farfarpsi1v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3) 
mup=0;
mue1=-1.5;
mue2=0.5;

%psi0_gauss = @(s1,s2,s3, alpha1,alpha2,alpha3) ((1/(2*pi.*alpha1)).^(1/4)).*exp(-((s1+0).^2)/(4.*alpha1)) ...
%            .* ((1/(2*pi.*alpha2)).^(1/4))*exp(-((s2+1.5).^2)/(4.*alpha2)) ...
%            .* ((1/(2*pi.*alpha3)).^(1/4)).*exp(-((s3-.5).^2)/(4.*alpha3));
psi0_gauss=@(s1,s2,s3, alpha1,alpha2,alpha3) (gaussian(s1,alpha1,mup)).*(gaussian(s2,alpha2,mue1)).*(gaussian(s3,alpha3,mue2));
%recenter gaussians
psi0_mpp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*sin(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mpm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*sin(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mmp = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*cos(theta2).*sin(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);
psi0_mmm = @(s1,s2,s3,theta1,theta2,theta3,alpha1,alpha2,alpha3) cos(theta1).*cos(theta2).*cos(theta3).* psi0_gauss(s1,s2,s3,alpha1,alpha2,alpha3);

mppterm=cos(theta1).*sin(theta2).*sin(theta3);
mpmterm=cos(theta1).*sin(theta2).*cos(theta3);
mmpterm=cos(theta1).*cos(theta2).*sin(theta3);
mmmterm=cos(theta1).*cos(theta2).*cos(theta3);
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
transp_mpp = psi0_mpp(s1-t1,s2+t2,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_mpm = psi0_mpm(s1-t1,s2+t2,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_mmp = psi0_mmp(s1-t1,s2-t2,s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);
transp_mmm = psi0_mmm(s1-t1,s2-t2,s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3);

%Integral 1, electron 1 integral
%int1_mpp = -1*1i*omega/2*integral(@(sigma) psi0_mmp(s1-t1,(sigma),s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
%int1_mpm = -1*1i*omega/2*integral(@(sigma) psi0_mmm(s1-t1,(sigma),s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
%int1_mmp = -1*1i*omega/2*integral(@(sigma) psi0_mpp(s1-t1,(sigma),s3+t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);
%int1_mmm = -1*1i*omega/2*integral(@(sigma) psi0_mpm(s1-t1,(sigma),s3-t3,theta1,theta2,theta3,alpha1,alpha2,alpha3),s2-t2,s2+t2);

int1_mpp = -1*1i*omega/2*(mmpterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussint(s2-t2,s2+t2,alpha2,mue1)).*(gaussian(s3+t3,alpha3,mue2)));
int1_mpm = -1*1i*omega/2*(mmmterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussint(s2-t2,s2+t2,alpha2,mue1)).*(gaussian(s3-t3,alpha3,mue2)));
int1_mmp = -1*1i*omega/2*(mppterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussint(s2-t2,s2+t2,alpha2,mue1)).*(gaussian(s3+t3,alpha3,mue2)));
int1_mmm = -1*1i*omega/2*(mpmterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussint(s2-t2,s2+t2,alpha2,mue1)).*(gaussian(s3-t3,alpha3,mue2)));

%Integral 2, electron 2 integral
%int2_mpp = -1*1i*omega/2*integral(@(tau) psi0_mpm(s1-t1,s2+t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
%int2_mpm = -1*1i*omega/2*integral(@(tau) psi0_mpp(s1-t1,s2+t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
%int2_mmp = -1*1i*omega/2*integral(@(tau) psi0_mmm(s1-t1,s2-t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);
%int2_mmm = -1*1i*omega/2*integral(@(tau) psi0_mmp(s1-t1,s2-t2,(tau),theta1,theta2,theta3,alpha1,alpha2,alpha3),s3-t3,s3+t3);

int2_mpp = -1*1i*omega/2*(mpmterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussian(s2+t2,alpha2,mue1)).*(gaussint(s3-t3,s3+t3,alpha3,mue2)));
int2_mpm = -1*1i*omega/2*(mppterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussian(s2+t2,alpha2,mue1)).*(gaussint(s3-t3,s3+t3,alpha3,mue2)));
int2_mmp = -1*1i*omega/2*(mmmterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussian(s2-t2,alpha2,mue1)).*(gaussint(s3-t3,s3+t3,alpha3,mue2)));
int2_mmm = -1*1i*omega/2*(mmpterm.*(gaussian(s1-t1,alpha1,mup)).*(gaussian(s2-t2,alpha2,mue1)).*(gaussint(s3-t3,s3+t3,alpha3,mue2)));
%Final values
val_mpp = transp_mpp + int1_mpp + int2_mpp;
val_mpm = transp_mpm + int1_mpm + int2_mpm;
val_mmp = transp_mmp + int1_mmp + int2_mmp;
val_mmm = transp_mmm + int1_mmm + int2_mmm;

%Output
val = [val_mpp,val_mpm,val_mmp,val_mmm];
end