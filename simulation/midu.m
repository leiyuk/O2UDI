function rou = midu(x,y,z,t)
%不考虑r<6378e3即碎片运动到地表下的情况，此时运动会停止
r=sqrt(x^2+y^2+z^2);
if(r<200e3+6378e3)
    g=9.8;
    m=29*1.661e-27;
    k=1.381e-23;
    T=310.7;
    rou=1.293*exp(-m*g*(r-6378e3)/k/T);
else
    rou=3.6e-10*exp(-(r-200e3-6378e3)./(37400+0.1/2*(r-200e3-6378e3)));
end


F=min(0.6,0.6/600e3*(r-6378e3));
e_r1=[x;y;z]/norm(r);
%不妨假设初始时刻太阳的方向矢量为(1,0,0)，即初始时刻太阳位于春分点上
eps0=23.5/180*pi;
reh=[cos(2*pi/365/86400*t);sin(2*pi/365/86400*t);0];
e_r2=[1 0 0;0 cos(eps0) -sin(eps0);0 sin(eps0) cos(eps0)]*reh;
nm=30/180*pi;
e_rm=[cos(nm) -sin(nm) 0;sin(nm) cos(nm) 0;0 0 1]*e_r2;
rou=rou*(1+F*dot(e_r1,e_rm));
end

