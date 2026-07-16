function date_str = get_target_date(x)
    % GET_TARGET_DATE 计算从2009年4月20日起经过 x 天后的日期
    % 输入: x (0 到 3650 之间的整数)
    % 输出: date_str (字符串，形如 '2009-4-20')
    
    % 1. 输入合法性校验 (可选，确保输入符合你的 1-3650 范围要求)
    if ~isnumeric(x) || floor(x) ~= x || x < 0 || x > 3650
        error('输入必须是 0 到 3650 之间的整数。');
    end
    
    % 2. 定义第 0 天的基准日期
    base_date = datetime(2009, 5, 10);
    
    % 3. 加上指定的天数 x
    % 使用 days() 函数将整数转换为时间间隔，直接与 datetime 相加
    target_date = base_date + days(x);
    
    % 4. 格式化输出为字符串
    % 'yyyy-M-d' 表示年份4位，月份和日期不补零（如 2009-4-20）
    % 如果你想让月份和日期强行补齐两位数（如 2009-04-20），可以改为 'yyyy-MM-dd'
    target_date.Format = 'yyyy-MM-dd'; 
    
    % 将 datetime 类型转换为标准的字符数组 (char) 输出
    date_str = char(target_date);
end