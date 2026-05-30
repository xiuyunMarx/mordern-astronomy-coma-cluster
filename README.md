# DESI Coma 数据处理与可视化

本目录包含 `desi-coma.csv` 及若干用于筛选、聚类和可视化 DESI Coma 数据的脚本。当前 CSV 使用的主要列为：

- `mean_fiber_ra`: RA，单位为 degree
- `mean_fiber_dec`: Dec，单位为 degree
- `z`: redshift

当前 `desi-coma.csv` 已经被清洗为只保留 `z <= 0.08` 的数据点。

## 环境依赖

脚本需要 Python 3，并依赖：

```bash
pip install matplotlib numpy scikit-learn
```

## 脚本说明

### `plot_mean_fiber_ra_dec.py`

功能：绘制 `mean_fiber_ra` 与 `mean_fiber_dec` 的二维散点图。可通过 `--z-low` 和 `--z-high` 限制 redshift 范围，只显示指定 `z` 区间内的数据点。

默认输入：

```text
desi-coma.csv
```

默认输出：

```text
desi-coma-mean-fiber-ra-dec-scatter.png
```

基本用法：

```bash
python3 plot_mean_fiber_ra_dec.py
```

限制 `z` 范围：

```bash
python3 plot_mean_fiber_ra_dec.py --z-low 0.013 --z-high 0.033
```

指定输入输出：

```bash
python3 plot_mean_fiber_ra_dec.py \
  --input desi-coma.csv \
  --output desi-coma-mean-fiber-ra-dec-scatter.png \
  --z-low 0.013 \
  --z-high 0.033
```

### `plot_z_histogram.py`

功能：绘制 `z` 的 histogram。默认只保留 `z <= 0.08` 的数据点。

默认输入：

```text
desi-coma.csv
```

默认输出：

```text
desi-coma-z-histogram.png
```

基本用法：

```bash
python3 plot_z_histogram.py
```

修改 bins 数量：

```bash
python3 plot_z_histogram.py --bins 120
```

修改 `z` 上限：

```bash
python3 plot_z_histogram.py --max-z 0.05
```

指定输入输出：

```bash
python3 plot_z_histogram.py \
  --input desi-coma.csv \
  --output desi-coma-z-histogram.png \
  --bins 80 \
  --max-z 0.08
```

### `gmm_coma_center.py`

功能：对 `desi-coma.csv` 的 `mean_fiber_ra`、`mean_fiber_dec`、`z` 做三维 GMM/EM 聚类。用户需要提供初始 Coma center，脚本会把这个 center 作为一个 Gaussian component 的初始均值，并在 EM 迭代过程中优化该 center。

默认输入：

```text
desi-coma.csv
```

默认输出：

```text
desi-coma-gmm-center-cluster.csv
desi-coma-gmm-center-trace.csv
```

其中：

- `desi-coma-gmm-center-cluster.csv`: 优化后的 Coma component 对应的数据点
- `desi-coma-gmm-center-trace.csv`: 每次 EM 迭代后的 center、selected count、mean probability 等轨迹信息

基本用法：

```bash
python3 gmm_coma_center.py \
  --center-ra 194.746077 \
  --center-dec 27.914923 \
  --center-z 0.02325
```

常用参数：

- `--components`: GMM component 数量，默认 `4`
- `--max-iter`: EM 最大迭代次数，默认 `200`
- `--center-tol`: center 收敛阈值，单位是标准化后的特征空间，默认 `1e-4`
- `--min-prob`: 写入输出 CSV 的最低 posterior probability，默认 `0.0`
- `--covariance-type`: GMM covariance 类型，默认 `full`
- `--random-state`: 随机种子，默认 `42`

提高成员筛选置信度：

```bash
python3 gmm_coma_center.py \
  --center-ra 194.746077 \
  --center-dec 27.914923 \
  --center-z 0.02325 \
  --min-prob 0.8
```

修改 component 数：

```bash
python3 gmm_coma_center.py \
  --center-ra 194.746077 \
  --center-dec 27.914923 \
  --center-z 0.02325 \
  --components 5
```

指定输出文件：

```bash
python3 gmm_coma_center.py \
  --center-ra 194.746077 \
  --center-dec 27.914923 \
  --center-z 0.02325 \
  --output desi-coma-gmm-center-cluster.csv \
  --trace-output desi-coma-gmm-center-trace.csv
```


### `visualize_coma_cluster.py`

功能：可视化已经筛选出的 Coma cluster CSV。默认读取 `gmm_coma_center.py` 生成的 `desi-coma-gmm-center-cluster.csv`，输出一张包含 RA/Dec 散点图和 `z` histogram 的 PNG。RA/Dec 散点图按 `z` 着色，并用黑色 `x` 标出 cluster 中心。

默认输入：

```text
desi-coma-gmm-center-cluster.csv
```

默认输出：

```text
desi-coma-gmm-center-cluster-visualization.png
```

基本用法：

```bash
python3 visualize_coma_cluster.py
```

指定输入输出：

```bash
python3 visualize_coma_cluster.py \
  --input desi-coma-gmm-center-cluster.csv \
  --output desi-coma-gmm-center-cluster-visualization.png
```

调整 histogram bins 和散点大小：

```bash
python3 visualize_coma_cluster.py --bins 80 --point-size 6
```

## 工作流

1. 绘制 redshift 分布：

```bash
python3 plot_z_histogram.py
```

2. 查看 RA/Dec 空间分布：

```bash
python3 plot_mean_fiber_ra_dec.py --z-low 0.013 --z-high 0.033
```

3. 用初始 Coma center 运行 GMM/EM 聚类：
这一步建议初始的 gmm components 数量为5, 请根据实际情况调整. coma cluster 理论上应该是一个椭圆形.
```bash
python3 gmm_coma_center.py \
  --center-ra 194.746077 \
  --center-dec 27.914923 \
  --center-z 0.02325 \
  --components 5 
```

4. 检查输出：

```text
desi-coma-gmm-center-cluster.csv
desi-coma-gmm-center-trace.csv
```


5. 可视化筛选出的 Coma cluster：

```bash
python3 visualize_coma_cluster.py
```
