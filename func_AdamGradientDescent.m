

function [fin, x] = func_AdamGradientDescent(f, x,  grad,  a, eps, beta1, beta2)
% Adam：梯度下降 数值优化
% 输入	
    % f - 一个具有解析形式的函数。比如 f = @(X) 2*X ;
    % x - 待优化的变量
    % a - 步长
    % eps - epsilon，一个很小很小的值，放在分母，避免除以0
    % beta1 - Adam经验值参数
    % beta2 - Adam经验值参数
% 输出
	% x：优化一步之后的变量
    % fin：如果值为1，说明梯度太小（平滑），无法继续优化，则停止。否则继续。
    
    persistent k;		% persistent变量 相当于 C/C++ 静态变量
    persistent m;
    persistent v;
    persistent beta1_;
    persistent beta2_;
    
    if isempty(k)
        k = 1;
        m = grad;
        v = grad.^2;
        beta1_ = beta1;
        beta2_ = beta2;
    else
        k = k + 1;
        m = beta1*m + (1-beta1)*grad;
        v = beta2*v + (1-beta2)*grad.^2;
        beta1_ = beta1_ * beta1;
        beta2_ = beta2_ * beta2;
    end
    
    m_hat = m / (1-beta1_);
    v_hat = v / (1-beta2_);
    
    fin = abs(grad) <= 1e-3;
    x = x - a / (sqrt(v_hat) + eps) * m_hat;
end
