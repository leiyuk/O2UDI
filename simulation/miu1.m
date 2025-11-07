function y=miu1(x)
if(x<=-1.1)
    y=-0.6;
else
    if(x<0)
        y=-0.6-0.318*(x+1.1);
    else
        y=-0.95;
    end
end
end

