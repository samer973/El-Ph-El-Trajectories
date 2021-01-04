function val = psiArbiternon(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3,phi) 
% define some functions that are needed in the computation
% radical = @(s, t, sigma) sqrt(t.^2 - (s - sigma).^2);
% radicalp = @(s, t, sigma) (sqrt(t + (s - sigma)))/(sqrt(t - (s - sigma)));
% mybessel1 = @(z,s,t,omega) besselj(1,omega*radical(s, t, z));
% mybessel0= @(z,s,t,omega) besselj(0,omega*radical(s, t, z));
% specify the initial data for KG eq, coming from the initial data for
% Dirac's equation for the electron with mass omega that corresponds to the initial
% probability density function being normal distribution with mean zero and
% variance alpha, and the mixing angle theta
% val_ppp,val_ppm,val_pmp,val_pmm,val_mpp,val_mpm,val_mmp,val_mmm;
% if(s2+t2 > s3-t3)
%     val_mpp = 0;
%     val_mpm = 0;
%     val_mmp = 0;
%     val_mmm = 0;
%     val_ppm = 0;
%     val_pmm = 0;
%     val_ppp = 0;
%     val_pmp = 0;
% else
%     if(s1-t1 < s2+t2)
%         val_mpp = subsref(nearnearpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{1}}));
%         val_mpm = subsref(nearnearpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{2}}));
%         val_mmp = subsref(nearnearpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{3}}));
%         val_mmm = subsref(nearnearpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{4}}));
%     else
%         val_mpp = subsref(farfarpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{1}}));
%         val_mpm = subsref(farfarpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{2}}));
%         val_mmp = subsref(farfarpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{3}}));
%         val_mmm = subsref(farfarpsi1(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{4}}));
%     end
%     
%     if(s1+t1 > s3-t3)
%         val_ppm = subsref(nearnearpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{2}}));
%         val_pmm = subsref(nearnearpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{4}}));
%         val_ppp = subsref(nearnearpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{1}}));
%         val_pmp = subsref(nearnearpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3,phi), struct('type', '()', 'subs', {{3}}));
%     else
%         val_ppm = subsref(farfarpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{2}}));
%         val_pmm = subsref(farfarpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{4}}));
%         val_ppp = subsref(farfarpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{1}}));
%         val_pmp = subsref(farfarpsi2(s1,s2,s3,t1,t2,t3,alpha,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{3}}));
%     end
% end
        val_mpp = subsref(farfarpsi1v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{1}}));
        val_mpm = subsref(farfarpsi1v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{2}}));
        val_mmp = subsref(farfarpsi1v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{3}}));
        val_mmm = subsref(farfarpsi1v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{4}}));
        val_ppm = subsref(farfarpsi2v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{2}}));
        val_pmm = subsref(farfarpsi2v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{4}}));
        val_ppp = subsref(farfarpsi2v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{1}}));
        val_pmp = subsref(farfarpsi2v2(s1,s2,s3,t1,t2,t3,alpha1,alpha2,alpha3,omega,theta1,theta2,theta3), struct('type', '()', 'subs', {{3}}));
    
%Output
val = [val_mmm,val_mmp,val_mpm,val_mpp,val_pmm,val_pmp,val_ppm,val_ppp];
end