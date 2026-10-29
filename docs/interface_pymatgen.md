# Interface Builder 的 pymatgen 引擎

本模块依据本机 `interface` Python 环境中的 pymatgen 2026.9.24 / pymatgen-core
2026.9.23 源码设计。核心是 `pymatgen.analysis.interfaces.zsl.ZSLGenerator`，
通过 `CoherentInterfaceBuilder` 选择终止面并构建原子界面。

文献：Zur 与 McGill，*Lattice match: an application to heteroepitaxy*，
J. Appl. Phys. 55, 378–386 (1984)，[DOI](https://doi.org/10.1063/1.333084)。
本地参考文献位于用户提供的 `D:\zotero\download\model\Zur-1984\`。
Stradi (2017) 和 Inoue (2026) 用于对照方法范围；本引擎不声称实现这两篇论文的全部算法。
接口说明见 [pymatgen 官方文档](https://pymatgen.org/pymatgen.analysis.interfaces.html)。

## 原生匹配参数

| 原生参数 | 默认值 | 界面含义 |
| --- | --- | --- |
| `max_area` | 400 Å² | 搜索两侧未应变超晶格的面积上界，原生实现使用严格小于上界的倍数 |
| `max_area_ratio_tol` | 0.09 | 对原胞面积比与整数面积倍数比进行预筛选；不是最终面积应变或长度应变阈值 |
| `max_length_tol` | 0.03 | 两个约化基矢的相对长度差容差，0.03 表示 3% |
| `max_angle_tol` | 0.01 | 两套基矢夹角的相对差容差，0.01 表示 1%，不是 0.01° |
| `bidirectional` | false | film→substrate 或反向比较任一满足即可（逻辑 OR） |
| 调用参数 `lowest` | false | 仅保留原生枚举顺序中的第一个匹配；不是能量最优，也不承诺全局最小应变 |

单向长度比较为 `|length_substrate / length_film - 1|`，夹角比较为
`|angle_substrate / angle_film - 1|`。开启双向比较时，也接受反向比值满足容差。
两侧晶面由用户直接指定；软件不额外扫描一个独立的 twist 角度网格。
旋转关系隐含在匹配的二维整数超胞基矢中。

`SubstrateAnalyzer` 还提供自动晶面枚举的 `film_max_miller`、
`substrate_max_miller`，以及可选弹性张量与参考能量输入。这些属于批量基底分析
接口，本模块采用指定晶面加 `CoherentInterfaceBuilder` 的建模流程，未提供
这些额外分析入口，也不将几何失配或 von Mises 应变称为能量。

## 终止面与原子建模参数

| 配置参数 | 默认值 | 含义 |
| --- | --- | --- |
| `substrate_plane` / `film_plane` | (001) / (001) | 各自输入晶胞基底中的接触晶面，不自动标准化输入晶胞 |
| `substrate_termination_ftol` / `film_termination_ftol` | 0.25 Å | 原生 `termination_ftol=(film, substrate)`：区分终止原子平面的距离容差 |
| `label_index` | true | 界面默认保留带编号的终止面标签，避免不同终止面因同名标签被覆盖；原生默认 false |
| `filter_out_sym_slabs` | true | 原生 slab 对称等价过滤；关闭可以探索更多终止面组合 |
| `termination` | null | 搜索后选择 `[film 标签, substrate 标签]`，null 在后端使用第一个组合 |
| `film_thickness` / `substrate_thickness` | 3 / 3 | 原生最小 slab 厚度，实际厚度由完整的 slab 单元决定 |
| `in_layers` | true | true 使用定向晶胞的单位晶面数，须为正整数；false 使用 Å |
| `gap` | 2 Å | 基底与薄膜的原子表面间距 |
| `vacuum_over_film` | 20 Å | 界面结构外侧的真空厚度 |
| `shift_mode` | `custom` | `custom` 自定义平移生成一个结构；`sampling` 从零起点进行网格采样 |
| `offset_a` / `offset_b` | 0 / 0 | 仅自定义模式使用，沿匹配后的 a/b 基矢移动薄膜，使用分数坐标，按周期取模 |
| `system` | 空字符串 | 导出 XYZ 的体系标签；为空则使用 substrate-film 的约化化学式 |

以基底为参考，薄膜面内晶格被强制匹配到基底，遵循 pymatgen 的共格界面构建
策略。这不是双侧弹性分担应变，也不执行原子弛豫或预测最低界面能。

## 应用增加的控制与可视化

- `max_matches=100`：限制保留的原生匹配候选数量；若仍有额外候选，明确提示
  已截断。仅对已保留的候选按薄膜面积、应变排序，不宣称全局最优。
- `match_index=0`：选择零起始候选。界面显示编号从 1 开始。
- 候选图用相同长度比例尺显示两套未应变约化超胞，并将各自第一基矢对齐；
  用于比较胞形和失配，不表示两种 bulk 晶体的绝对取向。
- 超胞矩阵包含约化基矢变换，基底是 pymatgen 生成的原始面内 slab 基矢，
  不是可直接作用在 bulk 三维晶格上的矩阵。原始笛卡尔基矢同时列出。
- 显示两侧面积、两个方向有符号长度差、相对角度差、绝对角度差以及薄膜
  Green–Lagrange 应变对应的 von Mises 应变。
- `offset_grid=[1,1]`：仅网格采样模式使用，生成 `n1*n2` 个结构。每方向不重复
  端点，采样 `(i/n1, j/n2)`，起点固定为 `(0,0)`，不叠加自定义偏移。
  这是匹配超胞采样，不是 DSC/CNID。自定义模式始终只生成一个结构。
- 两种模式仅显示各自的控件，切换时保留另一模式的输入。保存 JSON 时同时
  保留两组值，由 `shift_mode` 决定生效的参数；导出帧数遵循同一规则。
- 早期未保存 `shift_mode` 的 pymatgen 配置，若网格不是 `[1,1]` 且自定义平移
  周期等价于零，则识别为采样模式；网格为 `[1,1]` 时识别为自定义模式。
  若同时设置了非零自定义平移和多点网格，加载时要求显式选择 `shift_mode`，
  避免悄悄改变旧组合平移的含义。
- 搜索和建模均在后台运行。搜索参数变化会失效已有候选；改变几何参数或所选
  候选时复用小型搜索缓存。文件大小或修改时间变化也会失效缓存。
- 配置可以保存为 `interface.json`。原有 InterfaceMaster 专用配置不进行隐式
  转换。interface-master 是独立的 CSL/DSCL 引擎，详见
  [InterfaceMaster/Inoue 说明](interface_master.md)，不回退到 pymatgen 或旧实现。

主窗口原子预览、相机保留、POSCAR/CIF/XYZ 导出及多帧浏览沿用公共功能。
输出保持面内两根基矢为 a、b，堆叠/法向方向为 c，不再循环置换晶格行和分数坐标。
帧元数据通过 `interface_normal_axis=2` 告诉公共平面绘制函数使用第三分数坐标，
界面标记显示为 `fractional c`。旧帧缺少该字段时仍采用第一基矢的约定；
原子笛卡尔位置、匹配结果、平移和应变不因此改变。

## 自定义系数与 XYZ config_type

自定义单结构使用 `translation = d1 * a + d2 * b`，d1/d2 对应保存的
`offset_a`/`offset_b`。本引擎不计算 DSCL，因此始终使用匹配超胞的 a/b；
另一引擎的 DSCL 选择不会改变此处的系数含义。输入系数在元数据中保留原值，
只有构建所用的面内分数坐标按周期回绕。

单帧和多帧均可选择扩展 XYZ 导出。每帧第二行包含：

```text
config_type=Cu-Zr_zsl_match_1_strain_0.012_zspacing_2.0_d1_0.5_d2_0.3
```

格式为 `{system}_zsl_match_{编号}_strain_{von_Mises应变}_zspacing_{界面间距}_d1_{系数}_d2_{系数}`。
zspacing 取 `gap`（Å），不使用外部真空长度。
zspacing、d1、d2 在构型标签、数值控件和当前帧摘要中统一保留一位小数；
例如原始系数 `0.333333333333` 显示为 `0.3`。计算和 JSON 报告保留真实值，
未编辑输入内容时，读取或失焦不会把舍入后的显示值写回参数。
候选编号从 1 开始，应变为无量纲值（不是百分数）。采用 ZSL 的候选及应变标记，
不为其虚构一个未计算的 Sigma 或绝对旋转角。采样帧记录各自的 `(i/Na, j/Nb)`。
第二行的自定义属性仅写 `config_type`，另保留标准结构字段 `Lattice`、
`Properties` 和 `pbc`。基矢、帧编号、真空、终止面等内部元数据不额外写入
XYZ，仍保留在程序和 JSON 报告中。
