function isLeap = isLeapYear(year)
    % 判断是否为闰年
    % 输入：year - 要判断的年份
    % 输出：isLeap - 布尔值，如果是闰年则为true，否则为false

    if mod(year, 4) == 0
        if mod(year, 100) == 0
            if mod(year, 400) == 0
                isLeap = 1; % 能被400整除是闰年
            else
                isLeap = 0; % 能被100整除但不能被400整除不是闰年
            end
        else
            isLeap = 1; % 能被4整除但不能被100整除是闰年
        end
    else
        isLeap = 0; % 不能被4整除不是闰年
    end
end
