

% 作者：LT
% 本程序作用：
	% 异动点检测（LT 对 YSH方法 进行优化）

% 输入：
    % SP：特征矩阵(原始波形)。读取 SP 两个方法，参见 TXT
	% K1,K2, 异动点检测的2个阈值。这里 随机初始化。
		% K1,K2 越小，越敏感，检测出的异动点越多。
	% N：异动点检测的第3个阈值。
		% N 表示 保留 异动点 的比例。初始化为1，表示 保留 全部 异动点。
	
% 输出：
    % wave_filtered：滤波后波形(前24个为0)【包含数据(第一列)、时间戳(第二列)】
    % optimal_abnorm：异动点【包含数据(第一列)、时间戳(第二列)】	
	% 三个阈值 optimal_K1，optimal_K2，optimal_N
	

% 主要思想
	% 步骤：
		% A 预处理：滤波
		% B 第一次检测（baseline，YSH方法）
		% C 第二次过滤（LT方法）
	
	% A 预处理：滤波
		% 计算 原始波形SP的一阶差分（相邻两个元素差值 绝对值），
		% 计算 阈值 Δ： 比例因子*（最大值-最小值）
		% 如果 SP 中，波动量 小于 Δ 的 点很多，则 加大滤波强度（扩大 中值滤波 窗宽）

		
	% B 第一次检测：
	% K1和K2 初始化尽量小【0.001，0.005】，N初始化为1
		% 便于 baseline方法 敏感地更多的检测出 更多的 异动点。
		% 保留 这些 异动点，方便 第二次过滤。
    
	% C 第二次过滤		
		% 拿着B检测出来的异动点（此时K1和K2已固定），减少 异动点 个数，同时 保持 波形形状。
		% 设计损失，梯度下降 策略 优化 N。
			% 考虑因素有二：（思想由YSH提供）
				% 异动点比例N , 
				% 曲线（滤波后曲线，异动点描线）下面积之差。
			% 基于此，一共设计了 4 种损失函数。
				% 具体形式，请在本文件里搜索 关键词 '% 损失 todo'
			
		% 4 种损失 性能对比
			% 损失（1）和（4）目前是比较好用的。
			% 损失（2）和（3）不稳定，多数情况下，要么跑一次就结束了，要么N就一直停在初始值1。
		
		% 对 经过A预处理之后的 滤波波形，计算 二阶差分，并 升序排序。
			% 取出 二阶差分里面 后N %（取 大的） 的时间戳 index。
			% 在 滤波波形里 找到 index 对应的点，
			% 这些点，记作 当前 最优的 异动点序列。
			% 迭代结束，这些 异动点序列，就是 输出结果 optimal_abnorm。

% 程序开始
    clc;
    %    rng(1);
    rng('shuffle')        % 根据当前时间初始化生成器
    
    %     读取 MAT
    % load 'D:\2-代码工程\Z-Dataset\01-Internship\BaiDuYun\每一列特征-MAT\C_latest_parameter.mat'
    % todo 4
    load   'D:\2-代码工程\Z-Dataset\01-Internship\样本数据集1\6-parameter.mat';
    total_num = size(SP,1);
    ratio_1 = 0.0;
    ratio_2 = 0.2;
    Start = 1 + ratio_1 *total_num;
    End = round( ratio_2 *total_num );
    len = End - Start +1;
    SP = SP(Start: End, :);
    %     读取 MAT
    figure(1);
    plot(SP) ; % 查看 是否 掉点
	title('原始波形');

        
    % 最终得到的异动点，阈值
    % 返回值：
    optimal_K1=[];
    optimal_K2=[];
    optimal_N=[];
    optimal_abnorm = [];
    % 返回值
    
    % 阈值 初始化
    % todo1
    K1 = 0.001 +0.004*rand(1,1);
    K2 = 0.001 +0.004*rand(1,1);
    
    N = 1.0; 
    iteration = 100;
    
    %     K1=0.005; K2=0.05;    
    %     K1 = rand(1,1);  K2 = rand(1,1);    
    % 阈值
    
    % Adam参数
    beta1 = 0.9;
    beta2 = 0.999;
    eps = 1e-8;
    
    % todo2
    a = 0.01;
%     a = 0.001;
    % Adam参数
    
    % 损失
    loss = intmax;
    % 损失
    
    % (filter_num+1)/2 = 25前输出为0
    % 滤波窗宽
    filter_num=49;
    % 滤波窗宽

% 1预处理滤波
        % 	一阶滤波 得到 SP_1
        SP_1 = movmedian(SP, filter_num);
		% 对于matlab库函数（比如 movmedian），官网可以获取对应的C代码.

	ratio = 0.005;
	delta = ratio *(max(SP_1) - min(SP_1));
	SP_2 = SP_1;
	ratio_2 = 0.001;
	%
    diff_1 = abs( diff(SP_2) );
    index = find(diff_1<delta);     % 找出 diff_1 中 小于 delta 的 时间戳index			

	 if  size(index,1)  >= ratio_2 * len
            winSize = round(size(index,1)/80);
			SP_2 = movmedian(SP_2,  winSize );
     end
% 1预处理滤波

        % 	得到 滤波后 曲线 SP_2
        SP = SP_2;
        figure(2);
        plot(SP);       % 查看 预处理滤波 效果
		title('预处理滤波');

		clear func_AdamGradientDescent;         %  清空persistent 静态变量，以免 影响 下一轮 迭代。


        % 2 随机初始化K。
        % 		如果 异动点数量 小于 0.001*总数，则 重新 初始化K和N。
        step = 100;
        for  i = 1 : step
                [~, abnorms] = func_AbnormDetect(SP(:,1),K1,K2,filter_num);    
                if size(abnorms,1) >= 0.001*len
                    break;
                elseif  i==iteration
                    fprintf('初始化失败，请重新启动');
                    % break;
                    exit(1);
                else
                    K1 = 0.001 +0.099*rand(1,1);
                    K2 = 0.001 +0.099*rand(1,1);
                               N = 1.0;
                end
        end

% % 主函数
       % 检测
       % 滤波+异动点检测
       [wave_filtered, abnorms] = func_AbnormDetect(SP(:,1),K1,K2,filter_num);
       wave_filtered = wave_filtered(all(~isnan(wave_filtered),2),:);       % 删除含有NaN的行      
       % 对数据(第一列)做 二阶差分，保留 时间戳(第二列)
       second_order_dif = ( ( abs( diff(  abs(diff( abnorms(:, 1) ) ) ) ) ) )    ;       
       % 拼接 时间戳 （二阶差分输出 减少2个）
       second_order_dif = [second_order_dif, abnorms(3:end, 2) ];       
       % 对 二阶差分(第一列)排序，小→大
       sorted_dif = sortrows( second_order_dif, 1 );

       for  k = 1 : iteration
           % 取 后 N% 的 点
           cut_point = 1+round( (1-N)*size(sorted_dif, 1) );
           firstN_timestamp = sorted_dif( cut_point: end , 2);
           %  把 前一个点 也 上报
           firstN_timestamp = [firstN_timestamp;firstN_timestamp-1];
           % 去重
           firstN_timestamp=unique(firstN_timestamp,'rows');
           % 在 wave_filtered 里 取出 时间戳 为 firstN_timestamp 的 点
           firstN_abnorm = wave_filtered( firstN_timestamp+(filter_num+1)/2 , :);
           firstN_abnorm = sortrows(firstN_abnorm, 2);
           
		   % auc 是曲线下面积
           auc = trapz( wave_filtered(:,2) ,wave_filtered(:,1) );
           auc_hat = trapz( firstN_abnorm(:,2) ,firstN_abnorm(:,1) );
           rectangle =  (eps+(max(wave_filtered(:,1))-min(wave_filtered( (filter_num+1)/2:end ,1) ))*len);
           y        = auc /rectangle;
           y_hat = auc_hat /rectangle;
           
           % 损失 todo
           % （1）MSE  
			% f = @(N)    N + 1e+04*(  y - y_hat  )^2;
           
         % （2） 交叉熵  有点用  [有一些波形，动不了N，一直是1.0]
           CE = - y*log(y_hat)  - (1-y)*log(1-y_hat);
           f = @(N)    N * CE  ;
           
           % （3）两个log     [有一些波形，动不了N，一直是1.0]
%               f = @(N) log( eps  + N )/log(0.05) + ( log(  auc / auc_hat  ) )^2  ; 
           % （4）一个log
%              f = @(N)  N * ( log(  auc / auc_hat  ) )^2  ;  
     
           % 检测
           
           % GD
           fprintf("第%d次：N=%f, K1=%f, K2=%f, f(N)=%f\n", k, N,K1,K2, f(N));
           
           x= N;
           x_grad = (f(x+a) - f(x-a)) / (2*a);
           [x_fin, xn] = func_AdamGradientDescent(f, x,   x_grad,  a,  eps, beta1, beta2);                      
           N  = (xn);
           
           if N<=0 || N>1
               N = 1;
           end                      
           
           % GD
           
           %             if f(N,K1,K2) < loss
           loss = f(N);
           
           optimal_N=N;
           optimal_abnorm =  firstN_abnorm;
          
           my_plot( wave_filtered, optimal_abnorm);
           str = sprintf('N=%f, K1=%f, K2=%f.fig', N,K1,K2);
           path = strcat('用户六-损失1-中间结果','\', str);
           savefig('final.fig');
           %                 savefig(path);
           clf(figure(3));
           %             end
            
           if  size( optimal_abnorm,1 )/len   <  1000 / ( 24*60*60 )
               break;
           end           
           
           if x_fin
               fprintf('梯度太平滑，下降不动了。请重新启动。或者更换数据。');
               break;
           end
           
       end

delete(figure(1));
delete(figure(2));
delete(figure(3));

% % 程序结束
