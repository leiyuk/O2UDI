function y=miu3(x)
if(x<=-1.75)
    y=-0.3;
else
    if(x<-1.25)
        y=-0.3-1.4*(x+1.75);
    else
        y=-1;
    end
end
end

