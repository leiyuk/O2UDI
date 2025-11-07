function y=sigma1(x)
if(x<=-1.3)
    y=0.1;
else
    if(x<-0.3)
        y=0.1+0.2*(x+1.3);
    else
        y=0.3;
    end
end
end

