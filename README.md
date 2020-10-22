# Abnormal-Detection

### 一、任务目标
1. 异动点检测意思是，对于 一维 时序 波形数据，希望 用 最少的点，来 描绘 出一个 尽量 接近 原始的 波形轮廓  的 曲线。
2. 为什么用最少的点？由于通信带宽受限，我们需要用最少的点，来表示尽量完整的信息。也可以 理解为 数据压缩。
3. 例如：当曲线为一条直线时，只需两个点就可以把这条线画出来；当曲线为一条折线时，需要每个折线的转折点；当曲线为波动的曲线时，需要每一个极值点。

### 二、现有方法：twitter-AnomalyDetection
1. 工业界做异常检测（网络流量分析、金融交易等），一般是这样做的：
	- 无监督（统计方法居多。twitter的AnomalyDetection就是这样做的 ）→  
	- 得到标签 →  得到有标注训练集 → 
	- 有监督（机器学习/深度学习）

2. 机器学习/深度学习 方法，普遍问题是：
    - 无法得到准确标注的 训练集。除非人工标注（耗时耗力）
    - 非常吃数据。如果数据量不够多，很容易 过拟合。

3. 目前，基于无监督的方法（比如统计），最典型的就是  Grubbs Test 。
twitter方法【 [源码-R版本](https://github.com/twitter/AnomalyDetection)，[论文](https://arxiv.org/pdf/1704.07706.pdf)】，核心思想就是 Grubbs Test（很多地方用这个）：
	- 假设 原始数据服从 正态分布，构造 学生T分布  统计量C。其中，x_bar 和 s 表示时序数据X 的 均值 和 标准差。
	        - ![](https://img-1300025586.cos.ap-shanghai.myqcloud.com/abnormal-detection.png)
	- 如果满足：	
		- ![](https://img-1300025586.cos.ap-shanghai.myqcloud.com/abnormal-detection(1).png)
	- 就认为，当前测试的这个点，是异常的。
		- ![](https://img-1300025586.cos.ap-shanghai.myqcloud.com/abnormal-detection(2).png)


4. 但是，推特方法具有 [局限性](https://anomaly.io/anomaly-detection-twitter-r/) ，对于 数据集 要求较高。

### 三、本方法
#### 1. 思想
- 本方法 非常 手工，即，基于规则 + 数值优化 策略。算不上 无监督/有监督。
#### 2. 动机
- 由于 twitter方法 在 这些数据集上 不 work，只好另寻办法了。（一开始，我根据现有方法找案例，套模型，行不通。后来，调整思维方式：根据 现有数据集，想方法，调整方法，才慢慢有了效果…）
#### 3. 本方法步骤
1. 预处理
	- 清洗。工业界数据不是完美的，比如 掉点啊，NaN啊，之类的。如果连续掉点太多 或者	NaN，直接跳过。掉点不多就取前后均值（比较粗暴）。
	- 滤波。保留较大的波动，滤掉较小的波动（噪声）。对于时延较长的小幅度抖动，就可以平滑掉。
2. 异动检测。通过两个阈值K，得到第一次的异动点序列。第一次的序列一般较多，方便步骤3做过滤。
3. 数值优化。对于 第一次的异动点序列，进一步过滤，减少异动点。
	- 考虑 异动点个数比例N  + 波形曲线-异动点曲线 面积 差距  两个因素，设计损失
	- 运用梯度下降策略，优化N
	- 重复 3，直到一定的迭代次数


### 四、实现
1. 先看这个：		
	  `1. 代码+功能说明\先看这里！！代码文件说明.txt`

### 五、数据集+实验效果
 1. 数据集说明：
  - 用户用电 数据（比如 电流、功率）。
  - 在 /data 目录 下
  
 2. 效果
  - 在 /result 目录 下
  - 展示形式都是fig图： 滤波后曲线 上 描点（异动点）

