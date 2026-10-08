# 文献映射：保留具体区别，不把经典工具计作创新

2026-09-27，已读取以下作者稿HTML正文指定部分。

1. Haimovich, Seeber, Aldana-López, Gómez-Gutiérrez，
   [Differentiator for Noisy Sampled Signals with Best Worst-Case Accuracy](https://arxiv.org/html/2106.05320)。
   阅读§II、§III-A式(7)–(13)、Remark 1、§III-B Theorem 1及后续解释。
   其LP给条件外界，定理证明跨数据最坏误差最优，不声称每条历史都精确。
   我们的一例条件收紧不推翻其定理；同数据对照已实现，并非只与最近样本差分比较。
   混合p/v/a、条件未来位置和方向量不同于该特定估计目标，但这种范围不同不自动
   构成TAC创新；还须证明结构性新结论或算法收益。

2. Haddad与Halder，
   [The Curious Case of Integrator Reach Sets, Part I: Basic Theory](https://arxiv.org/html/2102.11423)。
   阅读§III Theorem 1支持函数公式、相关定义与Appendix D多项式根/绝对积分证明。
   有界输入经积分核的支持函数和多项式分段积分已存在。本包在混合观测条件下
   引入有限乘子并输出原始/对偶证书；不能把无观测支持公式或每段有限根数叫新定理。

3. Nowzari与Cortés，
   [Team-triggered coordination for real-time control of networked cyber-physical systems](https://arxiv.org/html/1410.2298)。
   阅读§IV-A状态/控制承诺、§IV-C/D请求与警告、§VI延迟/丢包处理。
   发送控制集合与按需请求本身已存在。本包使用物理变化率约束下的实际历史，
   没有承诺未来保持控制，也没有重命名该协调协议。

另外打开Haddad等2023 LTI可达集作者稿并核对支持函数推导段；它仍是相关工具，
不计新查新通过证据。近期Motion-Aware CBF编队工作只看了大学仓储摘要，未将
其完整控制定理作为已逐式审核内容，也未声称复现它的实验。

本轮审计结论：**获得可信的连续信息计算模块与正面条件改善证据，未通过主创新
终审。** 不能因“都是集合估计”就逻辑否定全部后续算法，也不能因未找到逐字
相同的混合数据代码就宣布新颖。尚缺严格的最近集合估计/最优恢复算法对照、
请求前决策理论与真正的多机闭环结果。纯对偶化+标准插值安全界不独立立为主线。
