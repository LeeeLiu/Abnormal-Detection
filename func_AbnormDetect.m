

% func_AbnormDetect���춯����� baseline
% ԭ���ߣ�YSH
% �ڶ����ߣ�LT
	% ��YSH�汾���� ��װ+ɾ����������һ�� ����
	
function [wave_filtered, abnorms] = func_AbnormDetect(SP,K1,K2, filter_num)

% ���룺
    % SP����������(ԭʼ����)����ȡ SP �����������μ� TXT
	% K1,K2, �춯����Ĺؼ���ֵ����Χ�� 0~1
	% filter_num���˲�����ȡ 49
	
% �����
    % wave_filtered���˲�����(ǰ24����������Ϊ0)
    % abnorms���춯��
	% ����������(��һ��)��ʱ���(�ڶ���)



%����Ƶ��
fs=10000;
%ÿ�ܲ���������
fnumber=fs/50;
%���������ʱ���� D_T
D_T=1;
% ������������������1
SP_num=size(SP,2);
%�˲����г�ʼ��
SP_A_filter_queue=zeros(filter_num,SP_num+1);

% �춯��
data_sample_A=[];
% �˲���
SP_A_filter_list=[];

% ����ԭʼ�汾д����
% ������ ʹ��  �˲� �� ���� ��һ����
% for m= 1: size(SP,1)-1   

    % �˲�
    for m= 1: size(SP,1)
        index=m;
        %��������������
        SP_A_filter_queue_temp=[SP(m,1),index];

        %�˲����и���
        SP_A_filter_queue(1,:)=[];
        SP_A_filter_queue=[SP_A_filter_queue;SP_A_filter_queue_temp];

        SP_A_filter_list_temp=zeros(1,length(SP_A_filter_queue_temp));

        %�˲�
        % ÿһ��(����/����)ȡ��λ��
        for i=1:size(SP_A_filter_queue,2)-1
            SP_A_filter_list_temp(i)= median(SP_A_filter_queue(:,i));
        end

        SP_A_filter_list_temp(end)= index-(filter_num+1)/2 ;

        % ĳһ������ �˲���
        SP_A_filter_list=[SP_A_filter_list;SP_A_filter_list_temp];
        
    % �˲�
    
    % 10-07
        %��ʼ���ʼ��
        if m==25
            data_1_A=SP_A_filter_list_temp(1);
            index_St=m;
            continue
        end

        %��ʼ���ʼ��
        if m==26
            data_2_A=SP_A_filter_list_temp(1);
            index_End=m;
            data_sample_A_before=data_2_A;
            data_sample_A_before_temp=SP_A_filter_list_temp;
            %dtimeΪ�������ʱ��
            dtime=D_T*(index_End-index_St);
            data_slope_A=(data_2_A-data_1_A)/dtime;
            s_limit_up1_A=data_2_A+data_slope_A*dtime+K1*data_2_A;
            s_limit_down1_A=data_2_A+data_slope_A*dtime-K1*data_2_A;
            s_limit_up2_A=data_2_A+data_slope_A*dtime+K2*data_2_A;
            s_limit_down2_A=data_2_A+data_slope_A*dtime-K2*data_2_A;

            continue
        end

        %�춯����
        if m>=27
            I_A_filter=SP_A_filter_list_temp(1);

            if (I_A_filter<=s_limit_down2_A*(1-K2))||(I_A_filter>=s_limit_up2_A*(1+K1))

                %���(�����ǰ�춯�㼰�춯��֮ǰ�ĵ�)
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
            %����ǰ������Ϣ������ʷ��Ϣ��
            data_sample_A_before=I_A_filter;        
            data_sample_A_before_temp=SP_A_filter_list_temp;        
        end     
    end
    % 10-07

%ɾ���ϱ����ظ���
data_sample_A=unique(data_sample_A,'rows');
if size(data_sample_A,2) > 0
    data_sample_A=sortrows(data_sample_A,size(data_sample_A,2));
end

% ���������
% srate=sum(sum(data_sample_A(1,:)~=0))/size(SP_filter,1);

% SP_filter, data_sample ��һ��index
% ��һ�����ݣ��ڶ�����ʱ���
% wave_filtered =  SP_filter(:,1) ;
% abnorms =  data_sample(:,1) ;

wave_filtered =  SP_A_filter_list ;
abnorms =  data_sample_A ;

end
