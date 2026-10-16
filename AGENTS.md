# `visualize/` 可视化程序说明

本目录是 CCkit 的 Windows 桌面可视化程序。它使用 PySide6 构建界面，
使用 PyVista/VTK 渲染结构，并通过 PyInstaller 生成独立的 onedir 程序。
它不属于 `CCkit <command>` 命令分发器，也不应把 GUI 依赖加入命令行脚本的
公共依赖。

## 当前功能

- **OmniGB**：读取晶体结构文件和 `gbin.json`，调用
  `core/omnigb/omnigb_backend.py` 生成晶界预览或 CSL 扫描结果。
- **Interface Builder**：使用 `core/interface/interface_backend.py` 匹配两种
  材料的接触面，并使用同目录的 `interface_panel.py` 管理独立参数卡片。
- **wyckoff model**：使用 `core/wyckoff/wyckoff_model.py` 生成晶体，并使用
  同目录的 `wyckoff_panel.py` 管理空间群、晶格和位点输入。
- 右侧 OmniGB 卡片包含输入文件、模型摘要、`dp1`/`dp2`/translation 步长、
  build/scan 模式选择、旋转轴/晶面/角度/Sigma/扫描范围控件、JSON 编辑器、
  自动刷新，以及把当前编辑器状态另存为独立 JSON 快照的按钮。
- 左侧使用 PyVista/VTK 显示结构。个人偏好包括原子大小、全局成键范围、
  元素对成键规则、对称晶面颜色/透明度。
- 当前显示结构可以导出为 POSCAR 或 CIF；scan 模式可以导出成功结果 CSV。
  修改 JSON 后，重建预览会保留相机位置、方向和缩放。
- 所有可调数字控件和下拉选项控件都必须禁止鼠标滚轮交互：控件应消费
  `wheelEvent`，不能因滚轮改变数值、进入编辑状态或切换选项；需要修改时
  只能通过点击、键盘或步进按钮完成。

## 架构

`gb_viewer.py` 是应用壳层，负责 Qt 生命周期、文件选择、右侧卡片容器、
后台构建线程、相机状态和渲染。它不应继续为每个新功能复制一套主窗口。

`feature_registry.py` 提供与 Qt 无关的功能注册表：

- `ViewerFeature` 描述稳定的功能 key、显示名称、说明、脚本路径和可用状态；
- `FeatureRegistry` 维护功能顺序，并提供 `register()`、`get()`、`all()` 和
  `available()`；
- `create_default_registry()` 注册已连接的 **OmniGB**、**Interface Builder**
  和 **wyckoff model**，并保留 Surface Builder、Structure Analysis 等空卡片。

以后接入新脚本时，应先在 `feature_registry.py` 注册 `ViewerFeature`，再在
`gb_viewer.py` 中为该 key 添加独立的参数卡片和构建适配器。新功能不应把
OmniGB 的字段、文件名或建模逻辑散落到全局 UI 中。需要多个脚本时，优先让
每个功能实现统一的“输入结构 + 配置 → 预览结构”接口，后台线程只负责调度，
渲染层只接收结构对象。

## 核心代码目录与功能分组

以后开发可视化程序的核心代码，必须放在 `visualize/core/` 下，并按功能
创建独立子目录。当前目录组织为：

```text
core/
├── omnigb/
│   └── omnigb_backend.py
├── interface/
│   ├── interface_backend.py
│   └── interface_panel.py
└── wyckoff/
    ├── wyckoff_model.py
    └── wyckoff_panel.py
```

1. 新功能的建模算法、参数处理和构建适配器统一放在
   `core/<功能名>/`，独立的功能参数卡片模块也放在同一目录。
   全局应用壳层、功能注册表、启动和打包入口放在 `visualize/` 顶层。
2. 从 `Scripts/` 接入现有脚本时，将需要的实现复制到相应功能子目录，
   保留原命令行脚本；可视化程序使用 `core/` 下的实现。
   复制后的模块应采用独立文件名，避免与公开命令脚本重名，干扰
   `CCkit_PATHS` 的递归文件查找。
3. 可视化程序专用的核心模块从旧位置移入相应功能子目录，统一更新调用方，
   旧位置不保留重复实现或兼容转发文件。
4. 功能子目录包含 `__init__.py`，支持明确的模块导入。新增或迁移模块时，
   同步更新 `feature_registry.py`、主窗口导入、功能间依赖和
   `gb_viewer.spec` 的模块及数据收集配置。
5. 源码与独立程序采用一致的核心目录层级；通过文件路径加载的后端需要
   随包保留其源文件和相对目录，使功能间调用在打包后仍能正确定位。

## 新功能的独立界面与布局设计

以后新增可视化功能时，必须为该功能独立设计界面，并遵守以下要求：

1. 每个功能使用独立的参数卡片或面板，管理自己的输入、配置和交互状态；
   切换功能时保留各自填写的内容。共用主窗口、渲染和导出能力。
2. 按该功能的操作流程安排控件，明确区分输入、参数、重复条目、执行按钮
   和结果摘要；标签、单位及必要的输入提示应清楚易读。
3. 布局应紧凑、美观，统一对齐方式、控件间距和内边距。卡片容器应按当前
   显示功能的内容计算尺寸，避免隐藏页面撑高当前卡片或控件被纵向拉开，
   造成大块无意义留白。
4. 晶格参数等独立参数应逐项换行；可添加、删除的原子或位点等重复条目，
   每条使用独立的小边框或卡片，配有清晰的编号和删除操作。
5. 新功能应合理适应窗口缩放和条目增减；内容超出可用高度时通过滚动访问，
   避免控件遮挡、文本截断或主要操作难以找到。
6. 接入功能时同步完成布局设计和静态复查；用户要求界面验证时，检查功能
   切换、窄窗口及添加和删除条目后的显示效果。

## 入口和打包

- `run_gb_viewer.bat`：使用本地 Python 启动源码程序，并自动加入旁侧的
  InterfaceMaster checkout。
- `run_gb_viewer_exe.bat`：启动已打包的独立程序。
- `build_gb_viewer.ps1`：使用指定 Python 环境运行 PyInstaller。
- `gb_viewer.spec`：收集 Qt、VTK、PyVista、pymatgen、InterfaceMaster、
  功能注册表和 `core/` 下的功能模块及其运行时文件。
- 输出程序位于 `dist/CCkitGBViewer/`；复制程序时必须保留整个目录。

打包程序启动后不依赖目标机器上的 Python 或 InterfaceMaster。源码开发时
可使用 `D:\miniconda3\envs\interface\python.exe`，但不要把本机路径写入
程序逻辑或配置文件。

## Python 依赖自动安装与随程序打包

用户已持续授权：以后开发或修改本可视化程序时，遇到功能所需但尚未安装的
Python 软件包，应自动安装，并将其纳入可视化程序的独立发布包，无需再次
询问是否安装。

1. 先确认实际用于运行和打包可视化程序的 Python 环境；使用该解释器的
   `-m pip` 安装缺失的软件包及其必要依赖，避免安装到无关环境。
2. 同步更新 `requirements-gb-viewer.txt`，按兼容性需要记录版本约束。
3. 同步更新 `gb_viewer.spec` 的模块、数据文件和动态库收集配置；有必要时
   更新 `build_gb_viewer.ps1` 的依赖检查，确保软件包进入独立程序的运行时。
4. 新增或补齐依赖后，自动执行将这些依赖纳入 `dist/CCkitGBViewer/` 所需的
   打包。此类依赖打包属于上述持续授权，无需另行确认；交付时说明安装的
   软件包、目标 Python 环境和独立发布包是否已更新。
5. 如联网、环境写入或打包受权限限制，按工具要求申请必要权限；如仍受阻，
   明确报告阻塞原因及尚未完成的安装或打包事项。

## 人工打包步骤

在 Windows PowerShell 中，从仓库根目录执行：

```powershell
cd D:\Github_repositories\CCkit

powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\visualize\build_gb_viewer.ps1 `
  -PythonExecutable "D:\miniconda3\envs\interface\python.exe"
```

构建脚本会使用 `visualize\gb_viewer.spec`，并把结果写入：

```text
D:\Github_repositories\CCkit\visualize\dist\CCkitGBViewer\
```

独立程序是其中的 `CCkitGBViewer.exe`。也可以从仓库根目录用样例文件启动：

```powershell
.\visualize\run_gb_viewer_exe.bat `
  Scripts\model\ZrO2.cif `
  Scripts\model\gbin.json
```

发布或复制程序时必须保留整个 `dist\CCkitGBViewer\` 目录，不能只复制
`.exe`。只有依赖发生变化时才需要使用 `build_gb_viewer.ps1 -Clean`；普通源码
修改直接重新运行上面的命令即可。

## 修改和验证

**除非用户明确要求，开发过程中不得自动运行测试或 `--verify` 检查。
打包仅在用户明确要求，或为落实上述 Python 依赖自动安装与随程序打包规则时
自动执行。**
用户明确要求后，才按照下面的命令执行相应验证或构建；仅修改文档时不需要
运行 Python 检查。

修改前先读取仓库根目录的 `CONVENTIONS.md` 和 `PROJECT_GUIDE.md`，并检查
`git status`。不要提交 `build/`、`dist/`、`__pycache__/`、截图或临时验证
输出。

用户明确要求验证时，至少执行：

```powershell
D:\miniconda3\envs\interface\python.exe -m py_compile visualize\gb_viewer.py visualize\feature_registry.py
D:\miniconda3\envs\interface\python.exe visualize\gb_viewer.py -h
```

有可用样例时，还应运行源码程序的 `--verify` 检查真实窗口和预览；修改打包
相关内容后，应重新构建并对 `dist/CCkitGBViewer/CCkitGBViewer.exe` 做同样的
验证。验证应覆盖新功能卡片、结构渲染、导出和已有 OmniGB 建模路径。
