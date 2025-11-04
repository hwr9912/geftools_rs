import tifffile as tiff
import numpy as np
import os
import h5py
import gc
from stereo.tools.cell_cut import CellCut

cc = CellCut(cgef_out_dir="out/")
bgef_path = "out/Y00855N1.bgef"
mask_path = "out/Y00855N1_ssDNA_regist_mask.tif"
tmp_mask = "out/Y00855N1_ssDNA_regist_mask_roi_nocomp.tif"

# 1) 从 bGEF 读取 bin1 的坐标范围（表达的最小/最大 x,y）
with h5py.File(bgef_path, "r") as f:
    ds = f["/geneExp/bin1/expression"]
    minx = int(ds.attrs["minX"])
    miny = int(ds.attrs["minY"])
    maxx = int(ds.attrs["maxX"])
    maxy = int(ds.attrs["maxY"])

print("GEF ROI (x,y,w,h):", (minx, miny, maxx, maxy))

# 2) 读整幅 mask
mask = tiff.imread(mask_path)
if mask.ndim == 3:
    mask = mask[..., 0]
H_mask, W_mask = mask.shape
print("Mask shape:", mask.shape)

# —— 计算按 GEF ROI 的期望尺寸（含端）
exp_h = (maxy - miny) + 1
exp_w = (maxx - minx) + 1

# —— 边界检查
if not (0 <= minx <= maxx < W_mask and 0 <= miny <= maxy < H_mask):
    raise ValueError(
        f"GEF ROI 越界: minx={minx}, miny={miny}, maxx={maxx}, maxy={maxy}, "
        f"而 mask 尺寸是 ({H_mask},{W_mask})"
    )

# —— 严格按 [min, max]（含端）裁剪；切片右边界 +1 变成右开
mask_roi = mask[miny:maxy+1, minx:maxx+1]
print("ROI shape (should match exp_h, exp_w):", mask_roi.shape, "expected:", (exp_h, exp_w))

# —— 安全写 TIFF（无压缩、释放句柄）
mask_roi_u8 = np.ascontiguousarray(mask_roi.astype(np.uint8))
tiff.imwrite(tmp_mask, mask_roi_u8, compression=None)

# 5) 调用 CellCut
out_path = cc.cell_cut(bgef_path=bgef_path, mask_path=tmp_mask)
print("cell_cut 输出：", out_path)
