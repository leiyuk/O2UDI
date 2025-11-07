function [yspace,xspace,P_shang,P_xia] = spaceDistribution(lo,hi,dR,a,e)
%lo和hi分别为半径的上下限，dR为间隔，a、e为1维行向量，no为向量维度
%yspace大致代表着空间碎片在半径方向的一个分布
xspace=linspace(lo,hi-dR,(hi-lo)/dR)+dR/2;
yspace=zeros(1,(hi-lo)/dR);
[~,lie_a]=size(a);
[~,lie_e]=size(e);
if(lie_a~=lie_e)
    fprintf("error!spaceDistribution");
    pause;
end
P_shang=0;
P_xia=0;
for i=1:lie_a
    for j=1:(hi-lo)/dR
        r_down=lo+(j-1)*dR;
        r_up=lo+j*dR;
        if(a(i)*(1+e(i))>r_down&&a(i)*(1-e(i))<r_up)
            jifen_down=asin(max((r_down-a(i))/a(i)/e(i),-1));
            jifen_shang=asin(min((r_up-a(i))/a(i)/e(i),1));
            yspace(j)=yspace(j)+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
        end
    end
    if(a(i)*(1+e(i))>hi)
        jifen_down=asin(max((hi-a(i))/a(i)/e(i),-1));
        jifen_shang=asin(1);
        P_shang=P_shang+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
    end
    if(a(i)*(1-e(i))<lo)
        jifen_down=asin(-1);
        jifen_shang=asin(min((lo-a(i))/a(i)/e(i),1));
        P_xia=P_xia+((jifen_shang-jifen_down)-e(i)*(cos(jifen_shang)-cos(jifen_down)))/pi;
    end
end
xspace=xspace/1e3;%单位变为1km
end