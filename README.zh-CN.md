# bosanquetFoam（OpenFOAM v2306）

本 solver 从本机原版 `applications/solvers/multiphase/interFoam` 复制后修改。
原版的 VOF、动量方程和压力方程结构均保留，只在 `pEqn.H` 中替换毛细力项。

## 连通性方案

每次 alpha 方程更新后，程序执行以下操作：

1. 用 `alpha.water >= connectivityAlpha` 将单元划分为水/非水；
2. 用 OpenFOAM `regionSplit` 计算全局连通分量；
3. 从 `seedPatches` 邻接的水单元选出主供液水体；
4. 跨 processor 面同步连通关系，串行和并行使用同一全局分量编号；
5. 只在“主水体 + 配置的 cellZone + alpha 界面”上使用 Bosanquet 力；
6. 与 seed 断开的液滴或液团继续使用原版 interFoam CSF 力。

混合公式为：

```
Fcap = (1 - M) Fcsf + M Fbos
Fbos = -interpolate(pCapBos) snGrad(alpha.water)
```

其中 `M` 是 `mainMeniscusMask`。

## 编译和运行

```bash
source /usr/lib/openfoam/openfoam2306/etc/bashrc
wmake
bosanquetFoam -case /path/to/case
```

并行运行：

```bash
decomposePar -case /path/to/case -force
mpirun -np 4 bosanquetFoam -case /path/to/case -parallel
```

## constant/bosanquetProperties

将 `bosanquetProperties.template` 复制为 case 的
`constant/bosanquetProperties`，然后填写：

- `connectivityAlpha`：判定水单元用于连通性计算的 alpha 阈值；
- `seedPatches`：与主储液区接触的一个或多个边界 patch；
- `name`：区域的说明性名称；
- `wallPatch`：提供接触角边界条件的壁面 patch；
- `cellZone`：允许应用 Bosanquet 力的单元区域；
- `radius`：该区域的等效圆管半径，带长度量纲。

`wallPatch` 在 `0/alpha.water` 中必须使用继承自
`alphaContactAngleTwoPhaseFvPatchScalarField` 的接触角边界条件，例如
`constantAlphaContactAngle`。

程序会写出三个诊断场：`mainWaterMask`、`bosanquetRegionMask` 和
`mainMeniscusMask`，可在 ParaView 中直接检查选择结果。

## 已验证

- OpenFOAM.com v2306：编译成功；
- 串行一时间步：成功；
- 4 核并行一时间步：成功，连通性统计与串行相同；
- 测试 case 包含一个接触 seed 的主水体和一个断开的水岛；断开水岛未进入
  `mainWaterMask`。

注意：v2306 的 `-dry-run` simplified mesh 在含 cellZone 的测试中会错误报告
zone 不存在；正常运行不受影响。
