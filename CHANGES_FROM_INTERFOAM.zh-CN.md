# 相比原版 interFoam 的修改

对照基准：OpenFOAM.com v2306 的
`applications/solvers/multiphase/interFoam`。

## 保持原样的部分

- 原版 VOF alpha 输运、MULES、速度方程和 PIMPLE 循环；
- 原版 `p_rgh` 压力方程及非正交修正流程；
- 原版重力项、湍流模型、动态网格和 `fvOptions`；
- 与主供液水体断开的液滴仍使用原版 CSF 表面张力。

## bosanquetFoam.C

- 程序名由 `interFoam` 改为 `bosanquetFoam`；
- 增加接触角、表面张力、`regionSplit` 和并行面同步相关头文件；
- 在原版 `createFields.H` 后读取 Bosanquet 配置并建立诊断场；
- 初始化时以及每次 alpha 方程更新后重新检测水体连通性。

## createBosanquetFields.H（新增）

- 读取 `constant/bosanquetProperties`；
- 用 `alpha.water >= connectivityAlpha` 判定参与连通性检测的水单元；
- 用 OpenFOAM `regionSplit` 建立全局连通分量；
- 用 `seedPatches` 确定主储液区/供液柱；
- 同步 processor 边界面，在并行计算中合并跨分区水体；
- 按每个 region 的 `wallPatch`、`cellZone` 和 `radius` 计算
  `pCapBos = 2 sigma cos(theta) / R`；
- 写出 `pCapBos`、`mainWaterMask`、`bosanquetRegionMask` 和
  `mainMeniscusMask` 供 ParaView 检查。

## readBosanquetContactAngle.H（新增）

- 从对应 `alpha.water` 壁面接触角边界条件读取接触角；
- 对该 wall patch 按面面积加权平均；
- 如果 wall patch 不是接触角边界条件则立即报错，避免静默使用错误参数。

## pEqn.H

原版毛细力：

```text
Fcsf = mixture.surfaceTensionForce()
```

新增 Bosanquet 力：

```text
Fbos = -interpolate(pCapBos) * snGrad(alpha.water)
```

实际送入原版 `phig` 的毛细力：

```text
Fcap = (1 - mainMeniscusMask) * Fcsf
     + mainMeniscusMask * Fbos
```

因此，只有位于配置 cellZone 内、并且属于 seed 连通主水体的界面使用
Bosanquet；其余界面继续使用原版 CSF。`phig` 之后的压力求解代码未改动。

## Make

- 增加 `meshTools` 头文件和链接库，以使用 `regionSplit` 和 `syncTools`；
- 输出程序改为 `$FOAM_USER_APPBIN/bosanquetFoam`。

## Case 新增要求

- `constant/bosanquetProperties`；
- 至少一个有效 `seedPatches`；
- 每个 Bosanquet region 对应一个非空 `cellZone`；
- 每个 `wallPatch` 在 `0/alpha.water` 中使用接触角边界条件。
