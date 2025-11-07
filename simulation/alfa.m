function y=alfa(x)
if(x<=-1.95)
    y=0;
else
    if(x<0.55)
        y=0.3+0.4*(x+1.2);
    else
        y=1;
    end
end
end

