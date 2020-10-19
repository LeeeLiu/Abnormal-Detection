

% func_AbnormDetect：异动点检测的 baseline
% 原作者：YSH
% 第二作者：LT
	% 在YSH版本上做 封装+删减：仅分析一列 特征
	
function [wave_filtered, abnorms] = func_AbnormDetect(SP,K1,K2, filter_num)

% 输入：
    % SP：特征矩阵(原始波形)。读取 SP 两个方法，参见 TXT
	% K1,K2, 异动点检测的关键阈值，范围是 0~1
	% filter_num：滤波窗宽，取 49
	
% 输出：
    % wave_filtered：滤波后波形(前24个点特征量为0)
    % abnorms：异动点
	% 都包含数据(第一列)、时间戳(第二列)



%采样频率
fs=10000;
%每周波采样点数
fnumber=fs/50;
%特征矩阵的时间间隔 D_T
D_T=1;
% 特征种类数，这里是1
SP_num=size(SP,2);
%滤波队列初始化
SP_A_filter_queue=zeros(filter_num,SP_num+1);

% 异动点
data_sample_A=[];
% 滤波后
SP_A_filter_list=[];

% 这是原始版本写法。
% 这样做 使得  滤波 后 矩阵 少一个点
% for m= 1: size(SP,1)-1   

    % 滤波
    for m= 1: size(SP,1)
        index=m;
        %新增特征量矩阵：
        SP_A_filter_queue_temp=[SP(m,1),index];

        %滤波队列更新
        SP_A_filter_queue(1,:)=[];
        SP_A_filter_queue=[SP_A_filter_queue;SP_A_filter_queue_temp];

        SP_A_filter_list_temp=zeros(1,length(SP_A_filter_queue_temp));

        %滤波
        % 每一列(属性/特征)取中位数
        for i=1:size(SP_A_filter_queue,2)-1
            SP_A_filter_list_temp(i)= median(SP_A_filter_queue(:,i));
        end

        SP_A_filter_list_temp(end)= index-(filter_num+1)/2 ;

        % 某一列特征 滤波后
        SP_A_filter_list=[SP_A_filter_list;SP_A_filter_list_temp];
        
    % 滤波
    
    % 10-07
        %起始点初始化
        if m==25
            data_1_A=SP_A_filter_list_temp(1);
            index_St=m;
            continue
        end

        %起始点初始化
        if m==26
            data_2_A=SP_A_filter_list_temp(1);
            index_End=m;
            data_sample_A_before=data_2_A;
            data_sample_A_before_temp=SP_A_filter_list_temp;
            %dtime为采样间隔时间
            dtime=D_T*(index_End-index_St);
            data_slope_A=(data_2_A-data_1_A)/dtime;
            s_limit_up1_A=data_2_A+data_slope_A*dtime+K1*data_2_A;
            s_limit_down1_A=data_2_A+data_slope_A*dtime-K1*data_2_A;
            s_limit_up2_A=data_2_A+data_slope_A*dtime+K2*data_2_A;
            s_limit_down2_A=data_2_A+data_slope_A*dtime-K2*data_2_A;

            continue
        end

        %异动点检测
        if m>=27
            I_A_filter=SP_A_filter_list_temp(1);

            if (I_A_filter<=s_limit_down2_A*(1-K2))||(I_A_filter>=s_limit_up2_A*(1+K1))

                %输出(输出当前异动点及异动点之前的点)
                data_sample_A=[data_sample_A;data_sample_A_before_temp;SP_A_filter_list_temp];       
                data_1_A=data_sample_A_before;

                index_St=index-1;
                index_End=index;
                dtime=D_T*(index_End-index_St);
                data_slope_A=(I_A_filter-data_1_A)/dtime;
                s_limit_up1_A=I_A_filter+data_slope_A*dtime+K1*I_A_filter;
                s_limit_down1_A=I_A_filter+data_slope_A*dtime-K1*I_A_filter;
                s_limit_up2_A=I_A_filter+data_slope_A*dtime+K2*I_A_filter;
                s_limit_down2_A=I_A_filter+data_slope_A*dtime-K2*I_A_filter;

            elseif ((I_A_filter<s_limit_down1_A*(1-K1))&&(I_A_filter>s_limit_down2_A*(1-K2)))||((I_A_filter>s_limit_up1_A*(1+K1))&&(I_A_filter<s_limit_up2_A*(1+K2)))
                index_End=index;
                dtime=D_T*(index_End-index_St);
                data_slope_A=(I_A_filter-data_1_A)/dtime;
                k_dtime=D_T;
                s_limit_up1_A=I_A_filter+data_slope_A*k_dtime+K1*I_A_filter;
                s_limit_down1_A=I_A_filter+data_slope_A*k_dtime-K1*I_A_filter;
                s_limit_up2_A=I_A_filter+data_slope_A*k_dtime+K2*I_A_filter;
                s_limit_down2_A=I_A_filter+data_slope_A*k_dtime-K2*I_A_filter;

            else
                index_End=index_St+1;
                dtime=D_T*(index_End-index_St);

                s_limit_up1_A=I_A_filter+data_slope_A*dtime+K1*I_A_filter;
                s_limit_down1_A=I_A_filter+data_slope_A*dtime-K1*I_A_filter;
                s_limit_up2_A=I_A_filter+data_slope_A*dtime+K2*I_A_filter;
                s_limit_down2_A=I_A_filter+data_slope_A*dtime-K2*I_A_filter;                                  
            end        
            %将当前特征信息存入历史信息中
            data_sample_A_before=I_A_filter;        
            data_sample_A_before_temp=SP_A_filter_list_temp;        
        end     
    end
    % 10-07

%删除上报的重复点
data_sample_A=unique(data_sample_A,'rows');
if size(data_sample_A,2) > 0
    data_sample_A=sortrows(data_sample_A,size(data_sample_A,2));
end

% 计算采样率
% srate=sum(sum(data_sample_A(1,:)~=0))/size(SP_filter,1);

% SP_filter, data_sample 多一列index
% 第一列数据，第二列是时间戳
% wave_filtered =  SP_filter(:,1) ;
% abnorms =  data_sample(:,1) ;

wave_filtered =  SP_A_filter_list ;
abnorms =  data_sample_A ;

end
