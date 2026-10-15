"""波形光学示例库生成器 —— 公共工具（仅标准库）。

本模块提供：
  * 物理常量与解析公式（用于把「预期现象」里的数字真正算出来，而不是猜）；
  * Config / 场景 (scene) 构造助手；
  * 采样-分辨率自检（元件特征尺寸必须 >= 若干像素，否则内核里根本表示不出来）；
  * Registry：收集示例、按 (类别, 文件名) 去重、按配置内容去重。

物理约定全部来自 docs/PHYSICS.md：
  * 时间因子 exp(-i w t)，沿 +z 传播的平面波 e^{+ikz}；
  * 薄透镜 t = exp(-i k r^2/(2f))；镜面曲率相位 exp(-i k curv r^2)；
  * propagate 的距离是沿光束方向的**路程**（折返后仍取正值，相位继续累加）。
"""

import hashlib
import json
import math
import os

# --------------------------------------------------------------------------
# 物理常量 / 常用波长
# --------------------------------------------------------------------------
PI = math.pi

WL = {
    "405": 405e-9,     # 紫（半导体激光）
    "450": 450e-9,     # 蓝
    "473": 473e-9,     # 蓝
    "500": 500e-9,     # 青绿
    "532": 532e-9,     # 绿（倍频 Nd:YAG）
    "550": 550e-9,     # 绿（视见函数峰值）
    "589": 589e-9,     # 钠黄双线平均
    "600": 600e-9,     # 橙
    "633": 632.8e-9,   # HeNe 红
    "650": 650e-9,     # 红
    "700": 700e-9,     # 深红
    "780": 780e-9,     # 近红外（Ti:S）
    "850": 850e-9,     # 近红外（VCSEL）
    "1064": 1064e-9,   # Nd:YAG
    "1550": 1550e-9,   # 光纤通信
}

# --------------------------------------------------------------------------
# 解析公式（SI）
# --------------------------------------------------------------------------


def k_of(wl):
    return 2.0 * PI / wl


def fringe_spacing(wl, z, d):
    """双缝/双光束干涉条纹间距 Δx = λ z / d（小角近似）。"""
    return wl * z / d


def sinc_first_zero(wl, z, a):
    """单缝（宽 a）夫琅禾费第一暗点位置 x = λ z / a。"""
    return wl * z / a


def airy_radius(wl, f, D):
    """圆孔艾里斑第一暗环半径 r = 1.22 λ f / D。"""
    return 1.22 * wl * f / D


def airy_peak(wl, f, P, D):
    """无像差圆孔聚焦的轴上峰值强度 I = |P| (π R²) / (λ² f²) / π R² ... 见 docs/PHYSICS.md §9：
    Strehl 参考值 I_ideal = P·πR²/(λ²z²)，R=D/2。"""
    R = D / 2.0
    return P * PI * R * R / (wl * wl * f * f)


def talbot_distance(wl, period):
    """泰伯自成像距离 z_T = 2 Λ² / λ。"""
    return 2.0 * period * period / wl


def zone_radius(m, wl, f):
    """波带片第 m 环外半径 r_m² = m λ f + (m λ / 2)²。"""
    return math.sqrt(m * wl * f + (m * wl / 2.0) ** 2)


def zone_count(radius, wl, f):
    """外半径 radius 内包含的波带数（由 r_m² = mλf + (mλ/2)² 反解）。"""
    a = (wl / 2.0) ** 2
    b = wl * f
    m = (-b + math.sqrt(b * b + 4 * a * radius * radius)) / (2 * a)
    return int(math.floor(m))


def bessel_core_radius(wl, beta):
    """截断贝塞尔光束中心亮斑半径 ≈ 2.405 / (k sin β)。"""
    return 2.405 / (k_of(wl) * math.sin(beta))


def gauss_divergence(wl, w0):
    """高斯光束远场半角 θ = λ /(π w0)。"""
    return wl / (PI * w0)


def rayleigh_range(wl, w0):
    """瑞利长度 z_R = π w0² / λ。"""
    return PI * w0 * w0 / wl


def waist_at(wl, w0, z):
    """w(z) = w0 sqrt(1 + (z/z_R)²)。"""
    zr = rayleigh_range(wl, w0)
    return w0 * math.sqrt(1.0 + (z / zr) ** 2)


def gouy(wl, w0, z):
    """古伊相位 ζ(z) = atan(z / z_R)。"""
    return math.atan(z / rayleigh_range(wl, w0))


def lens_focus(a, R1, R2, n):
    """薄透镜造镜者公式（供描述用）。"""
    return 1.0 / ((n - 1.0) * (1.0 / R1 - 1.0 / R2))


def grating_order_pos(q, f, wl, period):
    """光栅 q 级在傅里叶透镜焦面的位置 x_q = q f λ / Λ。"""
    return q * f * wl / period


def grating_angle(q, wl, period):
    """光栅方程 sinθ_q = q λ / Λ（法向入射）。"""
    return math.asin(min(1.0, abs(q) * wl / period)) * (1 if q >= 0 else -1)


def malus(power, theta):
    """马吕斯定律 I = I0 cos²θ。"""
    return power * math.cos(theta) ** 2


def fp_fsr(wl, L, n=1.0):
    """法布里-珀罗自由光谱范围 Δλ = λ²/(2 n L)。"""
    return wl * wl / (2.0 * n * L)


def fp_finesse(R):
    """法布里-珀罗精细度 F = π√R/(1-R)。"""
    return PI * math.sqrt(R) / (1.0 - R)


def fp_transmission_peak(R):
    """法布里-珀罗峰值透射率（无吸收）→ 1。"""
    return 1.0


def wedge_fringe_spacing(wl, alpha, n=1.0):
    """楔形（夹角 alpha）等厚条纹间距 Δx = λ / (2 n alpha)。"""
    return wl / (2.0 * n * alpha)


def fresnel_number(a, wl, z):
    """菲涅耳数 N_F = a²/(λ z)。"""
    return a * a / (wl * z)


def edge_fringe_scale(wl, z):
    """直边衍射条纹尺度 ~ sqrt(λ z / 2)。"""
    return math.sqrt(wl * z / 2.0)


def tilt_shift(z, theta):
    """倾斜平面波的横移 z·tanθ。"""
    return z * math.tan(theta)


def speckle_size(wl, z, D):
    """散斑平均尺寸 ~ λ z / D。"""
    return wl * z / D


# --------------------------------------------------------------------------
# 元件 / 光源 构造助手
# --------------------------------------------------------------------------

def pos(x=0.0, y=0.0, z=0.0):
    return {"x": float(x), "y": float(y), "z": float(z)}


def vec(x=0.0, y=0.0, z=1.0):
    return {"x": float(x), "y": float(y), "z": float(z)}


def comp(cid, ctype, p, yaw=0.0, pitch=0.0, roll=0.0, shape=None, params=None, label=None):
    """定位元件。yaw/pitch/roll 为 0 时省略，保持 JSON 精简。"""
    d = {"id": cid, "type": ctype, "pos": pos(*p)}
    if label:
        d["label"] = label
    if yaw:
        d["yaw"] = float(yaw)
    if pitch:
        d["pitch"] = float(pitch)
    if roll:
        d["roll"] = float(roll)
    if shape is not None:
        d["shape"] = shape
    if params:
        d["params"] = params
    return d


def source(sid, stype, params, p, direction=(0.0, 0.0, 1.0), wavelength=None,
           group=None, label=None):
    d = {"id": sid, "type": stype, "params": params, "pos": pos(*p), "dir": vec(*direction)}
    if label:
        d["label"] = label
    if wavelength is not None:
        d["wavelength"] = float(wavelength)
    if group is not None:
        d["group"] = group
    return d


# —— 形状（shape） ——

def circle(r):
    return {"kind": "circle", "params": {"radius": float(r)}}


def square(w):
    return {"kind": "square", "params": {"width": float(w)}}


def rectangle(w, h):
    return {"kind": "rectangle", "params": {"width": float(w), "height": float(h)}}


def slit(w, h):
    return {"kind": "slit", "params": {"width": float(w), "height": float(h)}}


def double_slit(separation, width, height):
    return {"kind": "double_slit",
            "params": {"separation": float(separation), "width": float(width),
                       "height": float(height)}}


def ellipse(a, b):
    return {"kind": "ellipse", "params": {"a": float(a), "b": float(b)}}


def ring(rin, rout):
    return {"kind": "ring", "params": {"rin": float(rin), "rout": float(rout)}}


def triangle(r):
    return {"kind": "triangle", "params": {"radius": float(r)}}


def polygon(r, sides):
    return {"kind": "polygon", "params": {"radius": float(r), "sides": int(sides)}}


def star(r, inner, points):
    return {"kind": "star", "params": {"radius": float(r), "inner": float(inner),
                                       "points": int(points)}}


def cross(length, width):
    return {"kind": "cross", "params": {"length": float(length), "width": float(width)}}


def superellipse(a, b, order):
    return {"kind": "superellipse", "params": {"a": float(a), "b": float(b),
                                               "order": float(order)}}


def custom(vertices):
    return {"kind": "custom", "params": {"vertices": vertices}}


# --------------------------------------------------------------------------
# Config 构造
# --------------------------------------------------------------------------

DIFFICULTY = ("入门", "进阶", "高级")


class Registry(object):
    """示例收集器：文件名与配置内容双重去重。"""

    def __init__(self):
        self.items = []          # [(filename, meta_dict, config_dict)]
        self._names = set()
        self._bodies = set()

    def add(self, filename, category, difficulty, description, expected,
            config, tags=None, theory=None, sample=False):
        assert filename.endswith(".json"), filename
        assert difficulty in DIFFICULTY, difficulty
        base = filename.split("/")[-1]
        if base in self._names:
            raise ValueError("duplicate filename: " + filename)
        body = json.dumps({k: v for k, v in config.items() if k not in ("meta", "name")},
                          sort_keys=True)
        digest = hashlib.sha1(body.encode()).hexdigest()
        if digest in self._bodies:
            raise ValueError("duplicate config body: " + filename)
        self._names.add(base)
        self._bodies.add(digest)
        meta = {
            "name": config["name"],
            "file": filename,
            "zod": {
                "category": category,
                "difficulty": difficulty,
                "description": description,
                "expected": expected,
                "tags": list(tags or []),
            },
        }
        if theory:
            meta["zod"]["theory"] = theory
        if sample:
            meta["sample"] = True
        self.items.append((filename, meta, config))

    def counts(self):
        c = {}
        for _f, m, _c in self.items:
            cat = m["zod"]["category"]
            c[cat] = c.get(cat, 0) + 1
        return c


def build_config(name, size, width, wl, components, sources, *, polarized=False,
                 method="asm", evanescent="decay", meta=None, sample=False,
                 category=None, difficulty=None, description=None, expected=None,
                 tags=None, theory=None, bandlimit=(0.9, 0.05)):
    """组装一个完整、可直接 POST /api/simulate 的 Config JSON。"""
    cfg = {
        "name": name,
        "grid": {"size": int(size), "width": float(width)},
        "wavelength": float(wl),
        "polarized": bool(polarized),
        "method": method,
        "evanescent": evanescent,
        "evanescent_limit": 0,
        "backward_regularize": False,
        "tikhonov_alpha": 0,
        "bandlimit": {"fraction": float(bandlimit[0]), "sigma": float(bandlimit[1])},
        "sources": sources,
        "scene": {"components": components},
    }
    if meta is None:
        meta = {
            "zod": {
                "category": category,
                "difficulty": difficulty,
                "description": description,
                "expected": expected,
                "tags": list(tags or []),
            }
        }
        if theory:
            meta["zod"]["theory"] = theory
    if sample:
        meta["sample"] = True
    cfg["meta"] = meta
    return cfg


# --------------------------------------------------------------------------
# 常用布局
# --------------------------------------------------------------------------

SRC_Z = -0.05          # 所有示例统一把光源放在 z=-50 mm
DET_YAW = PI           # 面向 -z 入射来的光（探测器法线朝向来光）


def plane_source(wl=None, power=1e-3, tilt_x=0.0, tilt_y=0.0, z=SRC_Z,
                 group=None, label=None, sid="src"):
    params = {"power": power}
    if tilt_x:
        params["tilt_x"] = float(tilt_x)
    if tilt_y:
        params["tilt_y"] = float(tilt_y)
    return source(sid, "plane", params, (0.0, 0.0, z), wavelength=wl,
                  group=group, label=label)


def gauss_source(waist, wl=None, power=1e-3, x=0.0, y=0.0, z=SRC_Z,
                 tilt_x=0.0, tilt_y=0.0, polarization=None, group=None,
                 label=None, sid="src"):
    params = {"power": power, "waist": float(waist)}
    if x or y:
        params["x"] = float(x)
        params["y"] = float(y)
    if tilt_x:
        params["tilt_x"] = float(tilt_x)
    if tilt_y:
        params["tilt_y"] = float(tilt_y)
    if polarization:
        params["polarization"] = polarization
    return source(sid, "gaussian", params, (0.0, 0.0, z), wavelength=wl,
                  group=group, label=label)


def lg_source(waist, p, l, wl=None, power=1e-3, z=SRC_Z, label=None, sid="src"):
    return source(sid, "laguerre_gaussian",
                  {"power": power, "waist": float(waist), "p": int(p), "l": int(l)},
                  (0.0, 0.0, z), wavelength=wl, label=label)


def hg_source(waist, m, n, wl=None, power=1e-3, z=SRC_Z, label=None, sid="src"):
    return source(sid, "hermite_gaussian",
                  {"power": power, "waist": float(waist), "m": int(m), "n": int(n)},
                  (0.0, 0.0, z), wavelength=wl, label=label)


def bessel_source(beta, radius, wl=None, power=1e-3, z=SRC_Z, label=None, sid="src"):
    return source(sid, "bessel",
                  {"power": power, "beta": float(beta), "radius": float(radius)},
                  (0.0, 0.0, z), wavelength=wl, label=label)


def spherical_source(radius, converging=False, wl=None, power=1e-3, z=SRC_Z,
                     label=None, sid="src"):
    return source(sid, "spherical",
                  {"power": power, "radius": float(radius), "converging": bool(converging)},
                  (0.0, 0.0, z), wavelength=wl, label=label)


def detector(cid, z, shape=None, label=None, params=None, x=0.0, y=0.0, yaw=DET_YAW):
    """探测器：z 可以是标量（放在 (x,y,z)）或三元组 (x,y,z)。"""
    if isinstance(z, (tuple, list)):
        x, y, z = z
    return comp(cid, "sensor", (x, y, z), yaw=yaw, shape=shape, params=params, label=label)


def det_rect(cid, z, w, h, label=None, **kw):
    return detector(cid, z, shape=rectangle(w, h), label=label, **kw)


def det_circle(cid, z, r, label=None, **kw):
    return detector(cid, z, shape=circle(r), label=label, **kw)


def mirror(cid, p, reflectivity=1.0, curvature=0.0, tilt_x=0.0, tilt_y=0.0,
           yaw=0.0, shape=None, label=None):
    params = {"reflectivity": float(reflectivity)}
    if curvature:
        params["curvature"] = float(curvature)
    if tilt_x:
        params["tilt_x"] = float(tilt_x)
    if tilt_y:
        params["tilt_y"] = float(tilt_y)
    return comp(cid, "mirror", p, yaw=yaw, shape=shape, params=params, label=label)


def lens(cid, p, f, aperture=None, yaw=0.0, label=None, x=None, y=None):
    params = {"f": float(f)}
    if x is not None:
        params["x"] = float(x)
    if y is not None:
        params["y"] = float(y)
    sh = None
    if aperture:
        sh = circle(aperture)
    return comp(cid, "lens", p, yaw=yaw, shape=sh, params=params, label=label)


def cyl_lens(cid, p, f, rotation=0.0, aperture=None, yaw=0.0, label=None):
    params = {"f": float(f), "rotation": float(rotation)}
    sh = circle(aperture) if aperture else None
    return comp(cid, "cylindrical_lens", p, yaw=yaw, shape=sh, params=params, label=label)


def aperture(cid, p, shape, edge_sigma=None, rotation=None, yaw=0.0, label=None,
             parallel=None, x=None, y=None):
    params = {}
    if edge_sigma:
        params["edge_sigma"] = float(edge_sigma)
    if rotation:
        params["rotation"] = float(rotation)
    if parallel is not None:
        params["parallel"] = bool(parallel)
    if x is not None:
        params["x"] = float(x)
    if y is not None:
        params["y"] = float(y)
    return comp(cid, "aperture", p, yaw=yaw, shape=shape, params=params or None, label=label)


def beam_splitter(cid, p, reflectivity=0.5, yaw=-PI / 4, aperture=None,
                  phase=None, label=None):
    params = {"reflectivity": float(reflectivity)}
    if phase:
        params["phase"] = float(phase)
    sh = circle(aperture) if aperture else None
    return comp(cid, "beamsplitter", p, yaw=yaw, shape=sh, params=params, label=label)


# --------------------------------------------------------------------------
# 采样自检
# --------------------------------------------------------------------------

def check_sampling(dx, features, minimum=4.0):
    """断言：每个特征尺寸至少占 `minimum` 个像素。

    内核网格 dx = W/N；小于一个像素的结构根本不存在，2 像素只能勉强表示，
    4–5 像素才有可用的对比度（硬边无软化时）。生成器在产出前强制自检。
    """
    bad = [(nm, v, v / dx) for nm, v in features.items() if v / dx < minimum]
    if bad:
        raise AssertionError("under-sampled feature(s) dx=%.3g m: %s" % (dx, bad))


def window_ok(W, needed):
    """断言：全部需要的半宽都在窗口内（留 5% 余量）。"""
    half = W / 2.0 * 0.95
    bad = [(nm, v) for nm, v in needed.items() if v > half]
    if bad:
        raise AssertionError("outside window W=%.3g m: %s" % (W, bad))


WINDOWS = (0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.016, 0.02, 0.03, 0.04)


def pick_window(needed_half, candidates=None):
    """选一个最小的标准窗口，使 needed_half 落在其中（留 5% 余量）。"""
    for W in (candidates or WINDOWS):
        if needed_half <= W / 2.0 * 0.95:
            return W
    raise AssertionError("no standard window fits half-width %.4g m" % needed_half)


def pick_grid(feature, W, minimum=3.0, sizes=(512, 1024, 2048)):
    """选一个最小的标准网格 N，使 feature 至少占 minimum 个像素。"""
    for N in sizes:
        if feature / (W / float(N)) >= minimum:
            return N
    raise AssertionError("feature %.4g m too small for W=%.3g m at N<=%d" % (feature, W, sizes[-1]))


def grid_of(size, width):
    return size, width, width / size


# --------------------------------------------------------------------------
# 写盘
# --------------------------------------------------------------------------

def write_all(registry, outdir, index_name="index.json"):
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    index = {}
    for filename, meta, cfg in registry.items:
        path = os.path.join(outdir, filename)
        d = os.path.dirname(path)
        if d and not os.path.isdir(d):
            os.makedirs(d)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(cfg, fh, ensure_ascii=False, indent=2, sort_keys=False)
            fh.write("\n")
        index[filename] = meta
    with open(os.path.join(outdir, index_name), "w", encoding="utf-8") as fh:
        json.dump(index, fh, ensure_ascii=False, indent=2, sort_keys=True)
        fh.write("\n")
    return index


def wave_configs_only(registry):
    """library/wave/ 顶层（非子目录）的示例。"""
    return [x for x in registry.items if "/" not in x[0]]


def quantum_configs(registry):
    return [x for x in registry.items if x[0].startswith("quantum/")]
