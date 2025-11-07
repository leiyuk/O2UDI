function [lie_a_choose,lie_e_choose] = sousuo(aeda_sm,a,e,w)
%aeda_sm矩阵中，第一行是a，第二行是e，第三行是w
%sousuo函数返回值为最接近当前小碎片演化情况的aeda_sm矩阵的2*6行向量，第一行对应delta_a/s_m，第二行delta_e/s_m
%重点即在于这一算法，之后需要改进
%这里取为近地点差距的绝对值和远地点差距的绝对值的和最小的点
[~,max_lie]=size(aeda_sm);
chushai=ceil(0.02*max_lie);
cs=zeros(4,chushai);
mm=1;
for i=1:max_lie
    aeda_sm(4,i)=i;%标记一个列序号
end
for i=1:chushai
    p=abs(aeda_sm(1,i)*(1+aeda_sm(2,i))-a*(1+e))+abs(aeda_sm(1,i)*(1-aeda_sm(2,i))-a*(1-e));
    lie_min=i;
    for j=i+1:max_lie
        if(abs(aeda_sm(1,j)*(1+aeda_sm(2,j))-a*(1+e))+abs(aeda_sm(1,j)*(1-aeda_sm(2,j))-a*(1-e))<p)
            p=abs(aeda_sm(1,j)*(1+aeda_sm(2,j))-a*(1+e))+abs(aeda_sm(1,j)*(1-aeda_sm(2,j))-a*(1-e));
            lie_min=j;
        end
    end
    cs(:,mm)=aeda_sm(:,lie_min);
    aeda_sm(:,lie_min)=aeda_sm(:,i);
    aeda_sm(:,i)=cs(:,mm);
    mm=mm+1;
end
p=min(abs(cs(3,1)-w),2*pi-abs(cs(3,1)-w));
lie_min=1;
for i=2:chushai
    if(min(abs(cs(3,i)-w),2*pi-abs(cs(3,i)-w))<p)
        p=min(abs(cs(3,i)-w),2*pi-abs(cs(3,i)-w));
        lie_min=i;
    end
end
if(e>3e-3)
    lie_a_choose=cs(4,lie_min);
    lie_e_choose=cs(4,lie_min);
else
    lie_a_choose=cs(4,1);
    lie_e_choose=cs(4,lie_min);
end
end
