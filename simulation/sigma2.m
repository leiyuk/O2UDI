function y=sigma2(x)
if(x<=-0.5)
    y=0.5;
else
    if(x<-0.3)
        y=0.5-(x+0.5);
    else
        y=0.3;
    end
end
end

