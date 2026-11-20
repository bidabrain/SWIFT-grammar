# LRD_base_boost 崩溃：BH subgrid_mass 变成 NaN（"zero subgrid mass"）

记录日期：2026-11-20
源码：SWIFT-LRD，branch `LRD`，revision `747c2025`（运行用的 swift_mpi 就是这个版本编出来的）
配置：`--enable-mpi --enable-fof --with-subgrid=SPIN_JET_EAGLE-XL --with-hydro=sphenix --with-kernel=wendland-C2`
（**没有** `--enable-debugging-checks`，所以所有 `SWIFT_DEBUG_CHECKS` 检查都没启用）

## 1. 现象

运行目录：`grammar:/gpfs/dbi224/LRD/LRD_base_boost`，参数文件 `eagle_spin_jet.yml`，单 rank，32 线程。

三次作业（1318937、1323607、1325598）都在同一个地方崩溃：

```
black_holes/SPIN_JET/black_holes.h:black_holes_prepare_feedback():1013:
Black hole with id 31264365 tried to evolve spin with zero (or less) subgrid mass.
```

- 每次都跑了约 8.8 小时（墙钟约 31600–31740 s），崩在 a ≈ 0.07187（z ≈ 12.915），step ≈ 43848。
- 都是同一个 BH（id 31264365），说明问题可以确定性地复现。
- sbatch 没加 `--restart`，所以每次都从 IC 重新开始。
- 和 MPI、内存、节点都无关（MaxRSS 约 23 GB）。

## 2. 结论：subgrid_mass 是 NaN，SPIN_JET 的检查拦不住 NaN

- `black_holes.h:1009`：`if (bp->subgrid_mass > 0.) {...} else error(...)`，崩溃就是在这里报出来的。
- `black_holes.h:1094`：质量更新之后的保护写的是 `if (bp->subgrid_mass < 0.)`，触发时打 warning（输出到 stderr），并把质量重置成种子质量。
  **但 stderr 里一条 "negative mass" warning 都没有。** 所以进入第 1009 行时，质量不是负数，而是 **NaN**（极小概率是正好 0）。NaN 在 `<0` 和 `>0` 两处比较里都是 false。
- 同样拦不住 NaN 的地方还有：
  - spin 截断（`black_holes.h:1022` 附近的 `spin_final > 0.998 / < -0.998 / fabsf<0.01`）：NaN 直接写进 `bp->spin`
  - `black_hole_spinup_rate()` 开头的范围检查（`black_holes_spin.h:969`）
  - `src/minmax.h` 里的 `min`/`max` 宏：遇到 NaN 返回另一边的值。例如 `accr_rate = min(accr_rate, f_Edd*Edd)` 时，如果 accr_rate 是 NaN，会直接变成 Eddington 上限；`black_holes_compute_timestep()` 里 `max(bh_timestep, time_step_min)` 遇到 NaN 会返回 time_step_min
- 还没确定的是：**哪个量最先变成 NaN 或 inf**。需要复现才能定位（见第 5 节）。

### 已排除的路径

- **LRD 改动（提交 `69c7d370e`）**：只是在 `in_lrd_phase` 时跳过 jet 和 energy reservoir 的累积，再加几个诊断量（`eddington_fraction_unsuppressed`、`in_lrd_phase`、`lrd_confined_mass`），没有改 subgrid_mass、spin 或吸积率。
- **"jet 效率之和 >1，把质量吃成负数"**：按源码公式（float32）算 slim disc：
  - spin 的平衡值约 0.1–0.3（spinup_rate 在 a≈0.3 以上为负）
  - 在这个范围内，`1 - ε_r - ε_j - ε_w` 大于 0.7（f_Edd 取 4.6、20、100 都是如此）
  - 只有 a 大于约 0.9 时才会小于 0，而那条路径有 1094 行的保护和 warning
- **rho_gas == 0**：`prepare_feedback` 开头会直接 return。
- **合并**：snapshot 里 NumberOfMergers=0，而且合并的 spin 函数（`black_holes_spin.h:1293`）对 ≤0 的质量有单独的 error。

## 3. 崩溃前这个 BH 的状态（snap_0005.hdf5，a=0.06737，z=13.84）

全盒只有 2 个 BH。id 31264365 的数据如下（内部单位：1e10 Msun、Mpc、km/s，时间单位约 978 Gyr）：

| 量 | 值 |
|---|---|
| SubgridMasses | 1.666e-5（约 1.7e5 Msun）；DynamicalMasses 1.76e-3 |
| AccretionModes | 2（slim disc） |
| EddingtonFractions / Unsuppressed | 4.63 / 4.63 |
| AccretionBoostFactors | 30 |
| AccretionRates | 2.30 |
| Spins | 0.355 |
| Rad / Jet / Wind eff | 0.0045 / 0.075 / 0.123 |
| LRDConfinedFlags | 1（处于 LRD 阶段） |
| Energy / JetReservoir | 1010 / 856（AGN 事件数为 0） |
| TimeBins | 40 |
| FormationScaleFactors | 0.0651（种子 1e5 Msun） |
| GasDensities / SoundSpeeds / T | 39099 / 0.362 / 25903 K |

## 4. 关键线索：BH 时间步塌缩到下限

`timesteps.txt`：从 a≈0.068 一直到崩溃，**每一步只更新 1 个 gpart 和 1 个 bpart**（gas=0），时间步固定在 bin 34（dt≈6.5e-10，约 640 yr），连续约 43,800 步。8.8 小时几乎全耗在这里。

- `SPINJETAGN:minimum_timestep_yr: 1000` 约为 1.02e-9 内部单位，向下取整到 bin 34，也就是说 **BH 时间步被截在最小值**（`black_holes.h:49` 的 `black_holes_compute_timestep`，`max(bh_timestep, props->time_step_min)`）。
- 在 snap_0005 时这个 BH 还在 bin 40（约 4e-8）。按 dt_heat、dt_jet、dt_ang_mom 估算也是约 3e-8，一致。
  之后时间步缩小了 64 倍以上，说明 mdot 在短时间内暴涨，很可能顶到了 `max_eddington_fraction: 100`（Bondi ∝ M²、boost_alpha=30，在 LRD 阶段又没有反馈自我调节）。
- 时间步被截断后，每步的 spin 和质量变化不再受 1% 容差约束（`black_holes.h:1423` 起的 dt_ang_mom 逻辑）。NaN 最可能在这段时间产生。
- 另外：LRD 阶段不往 reservoir 里存能量，但 `compute_timestep` 仍然按 Energy_rate 和 Jet_rate 来限制 dt_heat 和 dt_jet。这些反馈实际上不会发生，却在限制时间步，这是一个可以改的不一致（主要影响性能）。

相关参数（`eagle_spin_jet.yml` 的 SPINJETAGN 部分）：
`with_boost_factor: 1, boost_alpha: 30, boost_beta: 0, max_eddington_fraction: 100, use_krumholz: 1, with_krumholz_vorticity: 0, use_multi_phase_bondi: 0, use_nibbling: 1, include_slim_disk: 1, use_jets_in_thin_disc: 1, accretion_efficiency_mode: Variable, use_lrd_confinement: 1, lrd_confinement_eddington_fraction: 1.0, minimum_timestep_yr: 1000`；TimeIntegration `dt_min: 1e-10`。

## 5. 下一步

1. **复现**：`restart/swift_000000.rst`（22:47 写入）离崩溃（23:03）只有约 16 分钟墙钟。
   用 `./swift_mpi --restart ...` 接着跑。注意：如果改了 bpart 结构，或者开了 debugging-checks 导致结构布局变化，旧 restart 就读不了，只能从 IC 重跑（约 8.8 小时）或从更早的 snapshot 起跑。
2. **加诊断**（不改 bpart 结构，保证 restart 兼容）：在 `black_holes_prepare_feedback()` 的以下位置用 `isfinite()` 检查，一旦不是有限值，就打印这个 BH 的全部输入（M_sub、rho_gas、cs、v_gas、σ、mdot（cap 前后）、Edd、f_Edd、accretion_efficiency、模式、spin、各 eff、spinup_rate、delta_m_0、delta_m_real、dt），然后 `error()`：
   - Bondi/Krumholz 的 accr_rate 算出来后（`black_holes.h` 约 780 行）
   - boost 之后、cap 之前（约 845 行）
   - 重新计算 accretion_efficiency 和 eddington_fraction 之后（约 890 行）
   - spin_final 算出来后（约 1012 行）
   - `subgrid_mass += delta_m_real` 之后（约 1075 行）
   - `black_holes_compute_timestep()` 的返回值
3. **可能的修复**（定位之后再选）：
   - 1094 行改成 `if (!(bp->subgrid_mass > 0.f))`，并对 NaN 单独 error 或修复
   - spin 截断前先检查 `isfinite(spin_final)`
   - 吸积率和 Eddington 比做 `isfinite` 保护（避免 min 宏把 NaN 吞成上限）
   - 时间步：LRD 阶段不用 dt_heat/dt_jet 限制；或者被截断时按比例限制本步的 Δm/M 和 Δa
   - 物理参数：boost_alpha=30 加 f_Edd 上限 100，会让种子 BH 在 LRD 阶段失控增长，需要评估是否合理

## 6. 运行环境（grammar）

- 环境脚本 `~/gcc_ompi5`：gnu14/14.2.0、openmpi5/5.0.7、hwloc/2.12.0、libfabric/1.18.0、gsl/2.8、hdf5/1.14.6、ucx/1.18.0-mt，加上 `~/local/omp5/{fftw3,GKlib,METIS,ParMETIS}`
- 作业脚本 `/gpfs/dbi224/LRD/LRD_base_boost/mpirun.sbatch`（1 节点、1 rank、32 线程、`--exclusive`）

## 7. 诊断结果（作业 1328342，从 restart 起跑，2026-11-20）

带诊断的 `swift_mpi`（工作区里 `black_holes.h` 的 `BH_CHECK_FINITE`）从 `restart/` 续跑，墙钟 697 s 后在 step 28994、a=0.07053 被拦下：

```
BH 31264365: non-finite bp->angular_momentum_direction[0] = -nan (black_holes.h:1093)
  M_sub=1.846301e+19 M_dyn=5.017666e-02 mdot=6.639668e+25 f_Edd=100 f_Edd_raw=1.138214e+20
  acc_eff=1 mode=2 lrd=1 spin=0.0924 eps_r=6.2e-06 eps_j=6.5e-03 eps_w=0.144 dt=6.364e-10
```

**根因不是数值 bug，而是 subgrid_mass 的指数失控增长，最后 float32 溢出：**

1. M_sub = 1.85e19 内部单位 ≈ **1.8e29 Msun**（M_dyn 只有 5e8 Msun）。f_Edd 一直顶在 `max_eddington_fraction: 100`。
2. 以 100 倍 Eddington 吸积时质量 e 折叠时间 ≈ t_Sal/(100·(1-ε)) ≈ (0.061 × 450 Myr)/85 ≈ 0.32 Myr。
   从 1.67e-5 到 1.85e19 是 ln(1.1e24) ≈ 55 个 e 折叠 ≈ 18 Myr；
   而 ~29000 步 × 626 yr ≈ 18 Myr。**完全吻合**：就是纯粹的 Eddington 上限指数增长。
3. NaN 的产生点：`black_hole_angular_momentum_magnitude()`（`black_holes_spin.h:127`）用 float 算
   `M_sub² · a · G/c`，M_sub² ≈ 3.4e38 ≥ FLT_MAX → inf；随后 `ang_mom_total / modulus` = inf/inf = NaN，
   NaN 再扩散到 spin、subgrid_mass，最终在原版里表现为 "zero subgrid mass"。
4. subgrid_mass 的增长不受可用气体约束（nibbling 只限制 M_dyn），所以没有任何东西阻止它。
   LRD 阶段又没有反馈，Bondi ∝ M² 的失控无法自我调节（f_Edd_raw 达到 1e20）。

结论：第 5 节里的 NaN 保护只能把崩溃变成更早/更清楚的报错，**真正需要改的是物理设定**：
- 在 LRD 阶段（或全局）限制 f_Edd 上限（例如 max_eddington_fraction 设为 ~10 或更低），或降低 boost_alpha；
- 和/或限制 M_sub 不超过 M_dyn（或某个与可用气体挂钩的上限）；
- 和/或让 LRD 阶段有一个退出条件（例如持续时间或累计质量上限），而不是只看 f_Edd_raw。
