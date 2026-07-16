function [day,num] = search(real,a,e,w)
%输出位置最相近的前半碎片所在天数和序号
%real矩阵中，第一个维度存储a/e/w，第二个维度存储所在天数，第三个维度存储碎片序号
%sousuo函数返回值为找到碎片的天数和序号
%取为近地点差距的绝对值和远地点差距的绝对值的和最小的点
[~,max_day,max_num]=size(real);
N=max_day*max_num; %数据的总个数
% chushai=ceil(0.02*N);  %0.05可以调节
% disp(chushai)
chushai=20;

cs=zeros(6,chushai);
mm=1;
data=zeros(6,N);
for i=1:max_day
    for j=1:max_num
        m=j*(i-1)+j;
        data(1:3,m)=real(1:3,i,j);
        data(4,m)=i;
        data(5,m)=j;
        data(6,m)=m;
    end
end
%所有数存到二维数组data里并排序

for i=1:chushai
    p=abs(data(1,i)*(1+data(2,i))-a*(1+e))+abs(data(1,i)*(1-data(2,i))-a*(1-e));
    num_min=i;
    for j=i+1:N
        if(abs(data(1,j)*(1+data(2,j))-a*(1+e))+abs(data(1,j)*(1-data(2,j))-a*(1-e))<p)
            p=abs(data(1,j)*(1+data(2,j))-a*(1+e))+abs(data(1,j)*(1-data(2,j))-a*(1-e));
            num_min=j;
        end
    end
    cs(:,mm)=data(:,num_min);
    data(:,num_min)=data(:,i);
    data(:,i)=cs(:,mm);
    mm=mm+1;%轨道最接近的几个放到前面来,再比较w
end

p=min(abs(cs(3,1)-w),2*pi-abs(cs(3,1)-w));
num_min=1;
for i=2:chushai
    if(min(abs(cs(3,i)-w),2*pi-abs(cs(3,i)-w))<p)
        p=min(abs(cs(3,i)-w),2*pi-abs(cs(3,i)-w));
        num_min=i;
    end
end

if(e>3e-3)
    day=cs(4,num_min);
    num=cs(5,num_min);
else
    day=cs(4,1);
    num=cs(5,1);
end

end