# Walkthrough: THz-Unified-Optimizer Fixes

All four tasks described in the implementation plan have been successfully executed and verified.

## Changes Implemented

1. **Excluded Series 4 & 5**: Added a check in [run_overnight_pipeline.py](file:///c:/THz-Unified-Optimizer/scripts/run_overnight_pipeline.py) for the number of angular points. Datasets with fewer than 5 points are skipped.
2. **Fixed D Bound**: Changed the upper bound for the `D` parameter in [optimizer_2d.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/optimizer_2d.py) to `p_fixed_um - 0.5` (which evaluates to $15.0$ µm since $P=15.5$ µm), avoiding the degenerate $D=P$ solution.
3. **Convergence Logging**:
   - Integrated callback logging into the 2D optimizer [optimizer_2d.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/optimizer_2d.py) (logs loss and parameters every 100 iterations) and added warning logs when parameters end up on the boundaries.
   - Added logging for the 1D stages in [optimizer_1d.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/optimizer_1d.py) along with boundary warnings.
4. **Split Amplitude & Phase Loss**: Implemented a weighted combination of amplitude MSE and phase MSE in the 2D objective function in [optimizer_2d.py](file:///c:/THz-Unified-Optimizer/unified_optimizer/optimizer_2d.py) with weights $W_{\text{amp}}=1.0$ and $W_{\text{phase}}=0.1$, using proper unwrapping via `arctan2`.

---

## Verification Results

### Excluded Series 4 and 5
The output log shows that `series4` and `series5` were skipped correctly:
```text
>>> ОБРАБОТКА СЕРИИ: series4 <<<
Загружено точек по углам: 1
[WARNING] [series4] ПРОПУСК: только 1 угловых точек (требуется >= 5). Серия исключена.

>>> ОБРАБОТКА СЕРИИ: series5 <<<
Загружено точек по углам: 1
[WARNING] [series5] ПРОПУСК: только 1 угловых точек (требуется >= 5). Серия исключена.
```

### D Bounds
The results in [overnight_results.json](file:///c:/THz-Unified-Optimizer/results/overnight_results.json) confirm that the parameter `D_eff_um` for 2D optimization is no longer degenerate ($D=15.5$ µm):
* **356att**: $D = 4.40$ µm
* **series1**: $D = 4.32$ µm
* **series2**: $D = 4.93$ µm
* **series3**: $D = 4.34$ µm
* **Global_Average**: $D = 4.93$ µm

### Convergence Logger Output
In `results/overnight_execution.log`, we can see detailed optimization logs:
```text
  [1D stage1] fun=0.0020, nit=23, x=[16.5, 9.675, 0.0, 1.014]
  [1D stage2] fun=128.8421, nit=10, x=[17.0, 9.175, 2.066, 3.0]
  [1D stage3] fun=113.9418, nit=26, x=[17.2, 8.975, 1.907, 3.0, 0.955, 0.0]
  [WARNING]   [1D] ВНИМАНИЕ: параметр 'P_um'=17.2000 на границе (np.float64(16.8), np.float64(17.2))!
...
  [2D iter= 100] loss=0.002635 | D_um=4.7634, loss_factor=0.4004, gamma=0.1014, angle_offset=0.0079, tau_ps=0.0285
  [2D iter= 200] loss=0.002241 | D_um=4.9326, loss_factor=0.3788, gamma=0.2456, angle_offset=0.0079, tau_ps=0.0337
  [2D iter= 300] loss=0.002157 | D_um=4.9331, loss_factor=0.4025, gamma=0.1001, angle_offset=0.0076, tau_ps=0.0346
  [2D] Завершено: success=True, nit=370, nfev=641
```
