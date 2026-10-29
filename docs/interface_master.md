# Interface Builder 的 interface-master 引擎

在 Interface Builder 下方选择 **interface-master · CSL / DSCL**。
所有新实现都位于本仓库的 `core/interface/`；旁侧 `interface_master`
仓库只作为依赖读取。没有修改其源码，也不调用会在当前目录写日志、CIF 或
POSCAR 的原生高层搜索/生成方法。

## 文献与实现范围

阅读了用户提供的三篇正文：

- Zur 与 McGill (1984)，*Lattice match: an application to heteroepitaxy*：
  二维平移晶格、整数超胞、面积与失配之间的取舍；原有 pymatgen 引擎继续使用 ZSL。
- Stradi 等 (2017)，*Method for determining optimal supercell representation
  of interfaces*：区分面内转动、均匀变形与超胞大小，本引擎显式显示转角和主应变。
- Inoue 等 (2026)，*Systematic modeling of atomic matching in incoherent and
  heterointerfaces*：采用其式 (5)–(6)、连分数/Farey 周期解释及式 (27)–(28)
  的二维共格指数。该指数是几何分类，不能代替结构弛豫或能量计算。

这是一套基于这些思想的有限搜索与建模流程，未宣称复现论文全部图表或补充表。
支持一般有序晶体的指定平面，未自动遍历所有 hkl、搜索重构、生成位错网络、
计算化学键能或预测热力学稳定性。部分占位结构须先确定具体原子占位。

## 平面、转角与整数超胞

用户输入 hkl 或标准六方 hkil，均相对于各自输入晶胞。程序先获取体相原胞，
由 InterfaceMaster 的 `get_primitive_hkl` 转换晶面，再用
`get_pri_vec_inplane` 得到真正的面内原始平移基。附带位点属性参与原胞识别。
SlabGenerator 的定向晶胞可能含多个面内原胞，构建时显式换基并去除完全重复的
周期像，不把偶然选择的较大表面晶胞当作最小 CSL。

两侧法向各自转到 +z，面内第一基矢各自对齐 +x，然后绕 +z 旋转 film。
因此 0° 的参考是这两套所选原始表面基的 a 方向，不是任意约定的体相 [100]。
负 Miller 指数保留法向意义。最终模型由 substrate 的上表面接 film 的下表面；
film 原子不做镜像操作，保留其手性。厚度的“单位晶面”以转换后原胞的
`d_hkl` 为准，是最小厚度要求，不是精确的原子层数。

用列基矢表示未应变的二维基底 S、F，对于采样角 θ：

```text
U = S^-1 R(θ) F
U' = round(q U) / q = P/q       q = 1, …, q_max
F' = S U'
D = F' [R(θ) F]^-1
```

这沿用 InterfaceMaster 近似 CSL 搜索中的有理矩阵思想，但完整遍历用户指定的
角度和分母网格、去除相同 P/q，并保留排序最好的有限个候选，而非找到首个即停止。
每个矩阵采用统一分母 q；不是任意独立分母组合的穷举，也不保证连续角度空间的
全局最优。旋转区间端点都采样；最后一步不足步长时直接取终点。

D 的奇异值减 1 是面内主应变；面积变化为 `|det(D)-1|`。
D 的极分解可能附带小的转动修正。列表同时记录采样角、修正和实际极分解转角；
修正须在用户上限内，实际角也须位于指定区间内。一般晶体不自动假设 60° 或
90° 周期，用户应按其晶体对称性选择搜索区间。

## CSL、DSCL 与 σ

`commensurability.py` 用整数最大公约数、二维 Hermite 标准形和对偶关系
精确计算 `Z² ∩ (P/q)Z²` 与 `Z² + (P/q)Z²`。
交晶格用 `(L0 ∩ L1)* = L0* + L1*` 构造。随后调用 InterfaceMaster 原生
LLL 约化；两侧超胞矩阵由同一个 CSL 反求并核实为整数。

```text
S Hs = F' Hf = C
Is = |det Hs|, If = |det Hf|
σ = sqrt(Is If) = sqrt(A_CSL / A_DSCL)
```

这些等式描述的是施加有限 film 变形后的可公度近似体系；不表示原始两个实数
晶格具有精确 CSL。绘图中点代表 Bravais 平移格点，原子基元在主视图建模后显示。
双晶格图每侧最多画 250 个格点；σ 和面积始终用完整晶格计算。
主图保持 CSL 和 DSCL 的真实相对尺度，附图单独缩放 DSCL 并标出两个基矢长度，
避免面积相差很大时只能看到一小块紫色线条。附图明确标记为独立比例。

选定模型构建时还调用 InterfaceMaster `DSCcalc`，为二维基底附加相同辅助法向，
单独计算二维晶格和，并检查其包含两套平移基且面积与精确结果一致。
若原生浮点运算失败，保留精确二维 HNF 结果，并在模型报告的 `dscl_method`
及界面摘要中说明原因。没有把任意三维 DSC 的投影误当成二维 DSCL。

## 匹配周期与 Farey 解释

对于两套面内原始 a 或 b 向量的长度 `dO`、`dS`，枚举互素正整数 M、N：

```text
Δ1 = (M dO − N dS) / (N dS)
Δ2 = |M dO − N dS| / (M dO + N dS)
```

M 是 film 周期数、N 是 substrate 周期数，`N/M` 逼近 `dO/dS`。
显示低失配候选和主连分数收敛分数，并用 Farey 父分数解释较大周期的分解，
如 `11/13 = 5/6 ⊕ 6/7`。这些是 1D 周期指南；对斜晶格或非零旋转，
两个独立的长度比不足以求解二维匹配，实际超胞仍来自完整二维矩阵。
它们不被称为已弛豫的原子结构单元。

## 面板参数

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `im_angle_start/stop/step` | 0° / 60° / 1° | 有限旋转搜索区间及步长，未自动进行晶体对称约化 |
| `im_rotation_tol` | 0.5° | 有理化允许的附加极分解转角 |
| `im_max_denominator` | 30 | 矩阵统一分母上界；同时用作 1D 周期指南的 M、N 上界 |
| `im_max_strain` | 0.03 | 两个主应变绝对值都不超过 3% |
| `im_max_area_change` | 0.05 | film 面内面积变化不超过 5% |
| `im_max_sigma` | 100 | 广义二维共格指数上界 |
| `im_max_index` | 400 | 任一侧 CSL 包含的面内原胞数上界 |
| `max_area` | 400 Å² | 限制 substrate CSL 及对应未应变 film 超胞面积 |
| `max_matches` | 100 | 完整有限搜索后保留的最好候选数 |
| `im_rank_by` | `sigma` | 按 σ、最大主应变或面积排序；不是按能量排序 |
| `im_shift_basis` | `dscl` | 网格采样使用 DSCL 或匹配超胞，与自定义平移的基矢选择独立 |
| `custom_shift_basis` | `supercell` | 自定义单结构平移使用 `dscl` 或匹配超胞 `supercell`；旧配置缺省时维持 a/b 分数的含义 |
| `system` | 空字符串 | XYZ 的体系标签；为空时使用 substrate-film 的约化化学式 |
| `im_contact_depth` | 1 Å | 选取接触面的近表面原子用于投影配准统计 |
| `im_registry_tolerance` | 0.3 Å | film 接触层原子投影靠近 substrate 投影的判据 |
| `im_min_distance` | 1 Å | 过近接触判据，报告满足条件的 film 原子数，不自动删除 |
| `im_max_atoms` | 50000 | 单个模型原子数上界 |

终止面聚类、厚度、gap、vacuum 和 shift 模式控件与另一个引擎共用；
ZSL 专用的长度/夹角容差、双向比较及 `lowest` 在本引擎下隐藏且不参与计算。
为限制资源，最多 250000 个角度×分母试探，定向 slab 中间结构最多 200000
原子，全部采样帧合计最多 2000000 原子；终止面枚举所用定向晶胞最多
5000 原子、两侧终止面组合最多 10000 对。超限时给出调整参数的错误。

## 原子模型与平移

Substrate 固定，film 通过已接受的变形闭合周期边界，法向不额外施加弹性应变。
自定义 shift 生成一个模型；网格采样生成 Na×Nb 个独立结构。
自定义模式中，`offset_a` 和 `offset_b` 分别是界面显示的 d1、d2 系数：
`translation = d1 * basis_vector_1 + d2 * basis_vector_2`。
`custom_shift_basis=dscl` 使用当前匹配候选的两根 DSCL 基矢，`supercell`
使用匹配超胞的 a/b。允许负系数；输入系数在内部元数据和 JSON 中保留原值。
DSCL 自定义平移和采样都先将系数转换成匹配超胞的 a/b 分数，只有面内坐标回绕。
DSCL 平移等价性首先描述理想平移晶格；实际元素基元和终止面须结合原子结构判断。

## 扩展 XYZ 的构型标签

导出单帧或全部帧时选择 **XYZ · extended format (*.xyz)**。每帧第二行包含
晶格、原子列定义、周期性，以及唯一的自定义属性 `config_type`：

```text
config_type=hZrO-Zr_sigma_6.00_theta_30.0_zspacing_2.0_d1_0.5_d2_0.3
```

格式为 `{system}_sigma_{sigma}_theta_{theta}_zspacing_{zspacing}_d1_{d1}_d2_{d2}`。
sigma 与 theta 来自实际选中的候选，theta 是最终旋转角，不是扫描网格角；
zspacing 是 `gap` 的界面间距（Å），不是外部真空长度或 slab 厚度。
网格采样的 d1/d2 为当前帧的 `(i/Na, j/Nb)`；自定义模式为用户输入系数。
sigma 保留两位小数，theta 保留一位小数。
zspacing、d1、d2 在构型标签、数值控件和当前帧摘要中统一保留一位小数，
包括整数值（如 `2.0`、`0.0`）。这仅影响文本显示；例如 d2 实际为
`0.333333333333` 时显示 `0.3`，建模和 JSON 报告仍使用原始值。
未修改输入内容时，失焦或读取控件不会将显示值写回计算参数。
XYZ 不额外写入 engine、frame、shift_basis、基矢矩阵、终止面或接触诊断等字段。
这些信息保留在程序的内存帧元数据及 JSON 报告中，供可视化和分析使用。
第二行仅保留标准 `Lattice`、`Properties`、`pbc` 及 `config_type`。
重新读取精简 XYZ 时不恢复未导出的 DSCL 基矢叠加。
界面摘要显示当前帧的 config_type，支持选中复制；标签随模型生成，
导出时使用对应帧的元数据，不从可能已改变的控件临时拼接。

每帧计算 film 接触层原子到 substrate 接触层投影的最近距离、投影重合比例，
以及整个 film 到 substrate 的最小三维距离和过近接触的 film 原子数。
周期像仅沿面内方向展开，并用最小奇异值下界判断未遍历的像是否还能改善结果，
不把斜晶格的分数坐标逐分量舍入当作精确最小像。
投影统计是 film→substrate 的定向最近邻，可有多个 film 原子对应同一个
substrate 原子，不要求一对一，也不区分元素种类。最高投影重合比例不等于
最低界面能；不据此修改原子、生成键或自动选择物理基态。

模型沿用主窗口 3D 渲染、相机、结构导出和帧滑块；输出保持 a、b 为面内基矢，
c 为法向，不再重排晶格行和分数坐标。帧元数据设置 `interface_normal_axis=2`，
公共平面绘制函数据此在固定的第三分数坐标上绘制界面，标记为 `fractional c`。
旧帧缺少该字段时继续按第一基矢解释，Bicrystal Builder 的约定保持不变。
平移、σ、转角、两侧指数和接触统计保留在内存及 JSON 报告中；扩展 XYZ
仍只写标准结构字段和 `config_type`。
左侧原子模型在界面平面上以紫色显示一个 DSCL 单元轮廓和 d1、d2 平移向量。
单元居中放置，使用真实长度；其基矢按当前结构的分数坐标保存，随帧切换正确
定位，不跟随 film 平移改变自身的晶格基。该叠加层随界面平面显示选项一起开关。
右侧独立比例附图继续保留，便于看清较小的 DSCL。旧 XYZ 若未保存基矢元数据，
不会猜测或绘制 DSCL；重新构建模型后即可显示。
**Save matching / model report…** 保存搜索候选及当前配置对应的建模报告；
未建模或配置已更改时 `model` 为 null，不会导出过期模型的诊断作为当前结果。

## 运行与发布

源码启动脚本将旁侧 `interface_master` 加入 PYTHONPATH；冻结包沿用已有的
InterfaceMaster 收集规则，并增加本仓库的新模块和 SciPy 空间搜索入口。
SciPy 已存在于本机运行环境，未安装新包；依赖清单记录其直接使用。
此前已对算法进行已知算例与实际结构核查；运行流程保持构建时核对，不在搜索
阶段增加强制检查。本次仅调整显示与恢复原流程，未自动运行测试或重新打包。
