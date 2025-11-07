function y=miu2(x)
if(x<=-0.7)
    y=-1.2;
else
    if(x<-0.1)
        y=-1.2-1.333*(x+0.7);
    else
        y=-2;
    end
end
end

