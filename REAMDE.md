项目简介
FPGAReader 是一个高性能的 FPGA 布局优化器，核心采用改进的 Fiduccia-Mattheyses (FM) 算法。
它读取布局和网表信息，通过迭代移动节点来最小化跨 Die 的互连线（Cut Size），同时满足 Die 的资源
容量约束（RAM/DSP）。

输入输出
输入
    布局文件：支持原始 .pl 或粗化后的 .pl 文件（格式同 Coarsing 工具）。
    网表文件：Verilog 风格的网表文件。
输出
    结果文件 (FPGA_FM_output.place)：
        格式：实例名 y坐标 Die索引 [FIXED]
    记录了 FM 优化后的最终布局位置。
统计 (.csv)：
    记录初始 Cut、最终 Cut、优化提升比例、运行时间以及各 Die 的资源利用率。

关键配置 (代码内常量)
    is_coarsing = true	判断传入的布局文件是否为已经粗化过的文件。
    UPDATE_FREQUENCY = 5000：每移动 5000 次重建一次增益桶（Gain Bucket），防止累积误差。
    NET_FANOUT_THRESHOLD = 100：忽略大于 100 个引脚的网络对邻居增益的影响。
    MAX_SLL_BETWEEN_DIES：Die 间最大允许的 SLL 连接数（软约束）。