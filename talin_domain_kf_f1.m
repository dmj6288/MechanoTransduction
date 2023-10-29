function forcerates=talin_domain_kf_f1(f,k0,L0,L2)

n_tmp=10000;
f_0=1000;

Kt=4.1;
A=0.8;
force=500;
fx2=@(x) (Kt./A).*(x./L2+1./(4.*(1-x./L2).^2)-1/4);
xf0=@(f) (L0.*coth(f.*L0./Kt)-Kt./f).*(1+f./f_0);
xfm0=zeros(n_tmp+1,2);
xfm1=zeros(n_tmp+1,2);
for i=1:n_tmp+1
    x=0+0.999999*(i-1)*L2/n_tmp;
    fxm2(i,:)=[fx2(x) x];%force at a given extension %extension
    f0=0+(i-1)*force/n_tmp;
    if f0==0
        xfm0(i,:)=[0 f0];
        xfm1(i,:)=[0 f0];
    else
    xfm0(i,:)=[xf0(f0) f0];
    end
end
extension0=max(xfm0(:,1));
xf2_spline=spline(fxm2(:,1),fxm2(:,2));
fx0_spline=spline(xfm0(:,1),xfm0(:,2));
%reverse function of force-extension
xf2=@(f) ppval(xf2_spline,f);
fx0=@(x) ppval(fx0_spline,x);
% force dependent folding rate
kf_f=@(k0,f) k0.*exp(-(integral(xf2,0,f,'ArrayValued',true)-integral(xf0,0,f,'ArrayValued',true))./Kt); 
forcerates= kf_f(k0,f);
