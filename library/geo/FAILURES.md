# library/geo 生成失败清单

由 `go run ./tools/genlib` 生成（确定性输出，无时间戳）；校验与索引由 `go run ./tools/libcheck` 完成。

## 计数口径

- 登记候选（进入库的候选条目）：480 = 写入 397 + 候选跳过 83
- 构造失败（求解或布局失败，未产生候选，未登记）：9
- 合并（参数指纹与已登记条目相同，不计入候选）：1
- 同名文件改名（参数不同，不计入失败）：0
- 失败记录合计：92 = 构造失败 9 + 候选跳过 83（逐条不重复）

## 候选跳过（已登记，测量 / 断言 / 描述复核未通过，未写入 JSON）

1. singlet-pcm-f25-nbk7.json: assertion: virtual-focus RMS 1.35046mm exceeds half the aperture 2.5mm
2. singlet-pcm-f50-nbk7.json: assertion: virtual-focus RMS 2.70092mm exceeds half the aperture 5mm
3. singlet-pcm-f100-nbk7.json: assertion: virtual-focus RMS 5.40184mm exceeds half the aperture 10mm
4. singlet-pcm-f150-nbk7.json: assertion: virtual-focus RMS 8.10277mm exceeds half the aperture 15mm
5. singlet-pcm-f200-nbk7.json: assertion: virtual-focus RMS 10.8037mm exceeds half the aperture 20mm
6. singlet-pcm-f300-nbk7.json: assertion: virtual-focus RMS 16.2055mm exceeds half the aperture 30mm
7. singlet-pcm-f500-nbk7.json: assertion: virtual-focus RMS 27.0092mm exceeds half the aperture 50mm
8. singlet-pcm-f1000-nbk7.json: assertion: virtual-focus RMS 34.0066mm exceeds half the aperture 60.024mm
9. singlet-pcm-f25-nsf5.json: assertion: virtual-focus RMS 1.37996mm exceeds half the aperture 2.5mm
10. singlet-pcm-f50-nsf5.json: assertion: virtual-focus RMS 2.75993mm exceeds half the aperture 5mm
11. singlet-pcm-f100-nsf5.json: assertion: virtual-focus RMS 5.51986mm exceeds half the aperture 10mm
12. singlet-pcm-f150-nsf5.json: assertion: virtual-focus RMS 8.27979mm exceeds half the aperture 15mm
13. singlet-pcm-f200-nsf5.json: assertion: virtual-focus RMS 11.0397mm exceeds half the aperture 20mm
14. singlet-pcm-f300-nsf5.json: assertion: virtual-focus RMS 16.5596mm exceeds half the aperture 30mm
15. singlet-pcm-f500-nsf5.json: assertion: virtual-focus RMS 27.5993mm exceeds half the aperture 50mm
16. singlet-pcm-f1000-nsf5.json: assertion: virtual-focus RMS 34.26mm exceeds half the aperture 60.024mm
17. singlet-pcmconc-f50-nbk7.json: assertion: virtual-focus RMS 2.96419mm exceeds half the aperture 5mm
18. singlet-pcmconc-f100-nbk7.json: assertion: virtual-focus RMS 5.92839mm exceeds half the aperture 10mm
19. singlet-pcmconc-f200-nbk7.json: assertion: virtual-focus RMS 11.8568mm exceeds half the aperture 20mm
20. singlet-pcmconc-f400-nbk7.json: assertion: virtual-focus RMS 23.7135mm exceeds half the aperture 40mm
21. singlet-menn-f50-nbk7.json: assertion: virtual-focus RMS 2.61095mm exceeds half the aperture 5mm
22. singlet-menn-f100-nbk7.json: assertion: virtual-focus RMS 5.22192mm exceeds half the aperture 10mm
23. singlet-menn-f200-nbk7.json: assertion: virtual-focus RMS 10.4438mm exceeds half the aperture 20mm
24. singlet-menn-f400-nbk7.json: assertion: virtual-focus RMS 20.8876mm exceeds half the aperture 40mm
25. singlet-menn-f50-f2.json: assertion: virtual-focus RMS 2.78809mm exceeds half the aperture 5mm
26. singlet-menn-f100-f2.json: assertion: virtual-focus RMS 5.57616mm exceeds half the aperture 10mm
27. singlet-menn-f200-f2.json: assertion: virtual-focus RMS 11.1523mm exceeds half the aperture 20mm
28. singlet-menn-f400-f2.json: assertion: virtual-focus RMS 22.3046mm exceeds half the aperture 40mm
29. achromat-dls-f50-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
30. achromat-dls-f60-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
31. achromat-dls-f150-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
32. achromat-solve-f100-f2p8-nbk7-f2.json: assertion: only 7 of 37 sampled rays reach the image plane; only 7 of 37 sampled off-axis rays reach the image plane
33. achromat-solve-f100-f4-nbk7-f2.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
34. petzval-f100-f5.json: assertion: only 19 of 37 sampled off-axis rays reach the image plane
35. petzval-f150-f5.json: assertion: only 19 of 37 sampled off-axis rays reach the image plane
36. petzval-f200-f5.json: assertion: only 19 of 37 sampled off-axis rays reach the image plane
37. petzval-f300-f5.json: assertion: only 15 of 37 sampled off-axis rays reach the image plane
38. petzval-f400-f5.json: assertion: only 15 of 37 sampled off-axis rays reach the image plane
39. petzval-f500-f5.json: assertion: only 15 of 37 sampled off-axis rays reach the image plane
40. symmetric-f50-f4.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
41. symmetric-f100-f4.json: measure: lens: no ray reached the image plane
42. symmetric-f150-f4.json: measure: lens: no ray reached the image plane
43. symmetric-f200-f4.json: measure: lens: no ray reached the image plane
44. symmetric-f300-f4.json: measure: lens: no ray reached the image plane
45. macro-f50-f4.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
46. macro-f100-f4.json: measure: lens: no ray reached the image plane
47. eyepiece-huygens-fe20.json: measure: lens: no ray reached the image plane
48. eyepiece-huygens-fe25.json: measure: lens: no ray reached the image plane
49. eyepiece-huygens-fe32.json: measure: lens: no ray reached the image plane
50. eyepiece-huygens-fe40.json: measure: lens: no ray reached the image plane
51. eyepiece-huygens-fe50.json: measure: lens: no ray reached the image plane
52. eyepiece-ramsden-fe15.json: measure: lens: no ray reached the image plane
53. eyepiece-ramsden-fe20.json: measure: lens: no ray reached the image plane
54. eyepiece-ramsden-fe25.json: measure: lens: no ray reached the image plane
55. eyepiece-ramsden-fe32.json: measure: lens: no ray reached the image plane
56. eyepiece-ramsden-fe40.json: measure: lens: no ray reached the image plane
57. eyepiece-abbe-fe6.json: measure: lens: no ray reached the image plane
58. eyepiece-abbe-fe10.json: measure: lens: no ray reached the image plane
59. eyepiece-abbe-fe12.5.json: measure: lens: no ray reached the image plane
60. eyepiece-abbe-fe15.json: measure: lens: no ray reached the image plane
61. eyepiece-abbe-fe20.json: measure: lens: no ray reached the image plane
62. eyepiece-abbe-fe25.json: measure: lens: no ray reached the image plane
63. eyepiece-huygens-pup5-fe20.json: measure: lens: no ray reached the image plane
64. eyepiece-huygens-pup5-fe32.json: measure: lens: no ray reached the image plane
65. eyepiece-abbe-short-fe9.json: measure: lens: no ray reached the image plane
66. eyepiece-abbe-short-fe16.json: measure: lens: no ray reached the image plane
67. micro-achr-40x.json: assertion: on-axis spot RMS 5.27607mm exceeds half the focal length; field spot RMS 5.30986mm exceeds the focal length
68. micro-achr-60x.json: assertion: on-axis spot RMS 31.8611mm exceeds half the focal length; field spot RMS 32.004mm exceeds the focal length
69. micro-lister-10x.json: assertion: only 19 of 37 sampled rays reach the image plane; only 15 of 37 sampled off-axis rays reach the image plane
70. micro-lister-20x.json: assertion: only 7 of 37 sampled rays reach the image plane; only 7 of 37 sampled off-axis rays reach the image plane
71. micro-compound-20x.json: assertion: virtual-focus RMS 2.15196mm exceeds half the aperture 4mm
72. micro-compound-40x.json: measure: only 2 of 21 rays reach the last surface
73. mirror-cassegrain-corrector-f800.json: assertion: only 4 of 21 traced rays survive
74. mirror-cassegrain-corrector-f1200.json: assertion: only 4 of 21 traced rays survive
75. mirror-fieldlens-f400.json: measure: only 0 of 21 rays reach the last surface
76. mirror-fieldlens-f800.json: measure: only 0 of 21 rays reach the last surface
77. relay-1to1-f50.json: measure: lens: no ray reached the image plane
78. relay-1to1-f100.json: measure: lens: no ray reached the image plane
79. relay-1to1-f200.json: measure: lens: no ray reached the image plane
80. relay-4f-f100-f100.json: assertion: virtual-focus RMS 7.68603mm exceeds half the aperture 12.5mm
81. telephoto-f150-ratio0p85.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
82. telephoto-f300-ratio0p75.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
83. telephoto-f500-ratio0p7.json: measure: lens: no ray reached the image plane

## 构造失败（求解或布局失败，未产生候选）

1. bend k=1.3: solveEFL: no sign change of the power residual for s in [1e-5.11, 1e0.89] (target -100 mm)
2. bend k=2: solveEFL: no sign change of the power residual for s in [1e-5.36, 1e0.64] (target -100 mm)
3. micro lister 40x: EFL residual 1.37e-06 after freezing (target 4, got 4.000005460802599)
4. reducer galilean 0.5x: solveAfocal: no sign change on [-30,210]
5. reducer galilean 0.25x: solveAfocal: no sign change on [-45,175]
6. reducer galilean 0.125x: solveAfocal: no sign change on [-52.5,157.5]
7. relay1to1 f=50: solveEFL: no sign change of the power residual for s in [1e-1.52, 1e4.48] (target 50 mm)
8. relay1to1 f=100: solveEFL: no sign change of the power residual for s in [1e-1.52, 1e4.48] (target 100 mm)
9. relay1to1 f=200: solveEFL: no sign change of the power residual for s in [1e-1.52, 1e4.48] (target 200 mm)

## 合并（参数完全相同，保留先出现的条目）

1. mirror-newtonian-f1500-f5.json 与 mirror-parabola-f1500-f4.json 参数指纹相同（同一光学设计），合并

## 同名改名（参数不同，文件名加 -altN 后缀）

（无）


## 本批未覆盖（留给后续批次）

- 非球面（偶次/奇次多项式、Forbes Q 多项式）：非球面单透镜、施密特校正板、施密特-卡塞格林 A4 校正板已从本批移除（geo 有偶次非球面，但本批只做球面与解析二次曲面）
- 里奇-克雷蒂安（双曲面，两镜圆锥常数联立求解）与格利高里式：本批移除。RC 需要圆锥求解；格利高里式的副镜共轭符号未核实，实测焦点与设计位置不符（焦点守卫可检出），待修正后再加入
- 自由曲面（XY 多项式）、泽尼克面、DOE / 二元光学、热光学（dn/dT、热膨胀）、偏振（琼斯矩阵）、杂散光（二次反射）：本批未涉及
- 倾斜 / 偏心（坐标断点）：牛顿式以光轴展开的主镜表示（未建 45° 折轴镜）；卡塞格林的中心遮挡未建模
- 公差进阶（灵敏度、逆灵敏度、补偿器）、膜系联合评价：本批未涉及
