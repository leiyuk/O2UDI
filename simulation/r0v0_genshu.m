function [a,e,w] = r0v0_genshu(r0,v0,miu)
%r0,v0均为列阵，miu为引力常数
zita=norm(v0)^2/2-miu/norm(r0);
a=-miu/2/zita;
h=cross(r0,v0);
e=sqrt(1-norm(h)^2/miu/a);
cosi=h(3)/norm(h);
sinoumu=h(1)/norm(h)/sqrt(1-cosi^2);
cosoumu=h(2)/norm(h)/(-sqrt(1-cosi^2));
e_s=1/miu*cross(v0,h)-r0/norm(r0);
ON=[cosoumu;sinoumu;0];
cosw=dot(ON,e_s)/norm(e_s)/norm(ON);
if(dot(cross(ON,e_s),h)<0)
    w=2*pi-acos(cosw);
else
    w=acos(cosw);
end
end

