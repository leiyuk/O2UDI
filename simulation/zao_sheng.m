function y = zao_sheng(sigma)%生成3*1的噪声矢量
while 1
    y=randn(3,1)*sigma;
    if(abs(y(1))<=3*sigma&&abs(y(2))<=3*sigma&&abs(y(3))<=3*sigma)
        return;
    end
end
end

