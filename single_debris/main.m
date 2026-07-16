clear;
% clf;

% 设定主文件夹路径
mainFolder = 'cosmos_1408';

% 获取文件夹信息
d = dir(mainFolder);

subFolders={};
for i=1:length(d)-2
    subFolders{i} = d(i+2).name;        % 提取名称
end

% 结果数组
% disp(subFolders);

R2_record={};
all_record={};
use_record={};
pointxy_record={};
linexy_record={};


num=653;
for i=1:num

    [R2,all,use,pointxy,linexy]=SingleDebris('cosmos_1408/',subFolders{i},...
    30,i*3,0);
    R2_record{end+1,1}      = R2;
    all_record{end+1,1}     = all;
    use_record{end+1,1}     = use;
    pointxy_record{end+1,1} = pointxy;
    linexy_record{end+1,1}  = linexy;

end



x=0;
for i=1:num
x=x+R2_record{i};
end
disp(['R2 average: ',num2str(x/num)])


x1=0;
for i=1:num
    x1=x1+abs(linexy_record{i}(2,1))/max(pointxy_record{i}(2,:));
end
disp(['y-axis relative error average: ',num2str(x1/num)])

x2=0;
for i=1:num
    x2=x2+use_record{i};
end
disp(['num of data average: ',num2str(x2/num)])

cal=0;
for i=1:num
    if linexy_record{i}(2,end)>linexy_record{i}(2,1)
        cal=cal+1;
    end
end
disp(['ratio of positive slope: ',num2str(cal/num)])

plot_NC;

