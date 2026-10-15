# library/geo 生成失败清单

由 `go run ./tools/genlib` 生成（确定性输出，无时间戳）；校验与索引由 `go run ./tools/libcheck` 完成。

## 计数口径

- 登记候选（进入库的候选条目）：676 = 写入 624 + 候选跳过 52
- 构造失败（求解、布局或质量门槛失败，未产生候选，未登记）：13
- 合并（参数指纹与已登记条目相同，不计入候选）：1
- 同名文件改名（参数不同，不计入失败）：0
- 失败记录合计：65 = 构造失败 13 + 候选跳过 52（逐条不重复）

## 候选跳过（已登记，测量 / 断言 / 描述复核未通过，未写入 JSON）

1. achromat-dls-f50-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
2. achromat-dls-f60-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
3. achromat-dls-f150-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
4. achromat-solve-f100-f2p8-nbk7-f2.json: assertion: only 7 of 37 sampled rays reach the image plane; only 7 of 37 sampled off-axis rays reach the image plane
5. achromat-solve-f100-f4-nbk7-f2.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane
6. petzval-f100-f5.json: assertion: only 19 of 37 sampled off-axis rays reach the image plane
7. petzval-f150-f5.json: assertion: only 19 of 37 sampled off-axis rays reach the image plane
8. petzval-f200-f5.json: assertion: only 19 of 37 sampled off-axis rays reach the image plane
9. petzval-f300-f5.json: assertion: only 15 of 37 sampled off-axis rays reach the image plane
10. petzval-f400-f5.json: assertion: only 15 of 37 sampled off-axis rays reach the image plane
11. petzval-f500-f5.json: assertion: only 15 of 37 sampled off-axis rays reach the image plane
12. symmetric-f50-f4.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
13. symmetric-f100-f4.json: measure: lens: no ray reached the image plane
14. symmetric-f150-f4.json: measure: lens: no ray reached the image plane
15. symmetric-f200-f4.json: measure: lens: no ray reached the image plane
16. symmetric-f300-f4.json: measure: lens: no ray reached the image plane
17. macro-f50-f4.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
18. macro-f100-f4.json: measure: lens: no ray reached the image plane
19. eyepiece-huygens-fe20.json: measure: lens: no ray reached the image plane
20. eyepiece-huygens-fe25.json: measure: lens: no ray reached the image plane
21. eyepiece-huygens-fe32.json: measure: lens: no ray reached the image plane
22. eyepiece-huygens-fe40.json: measure: lens: no ray reached the image plane
23. eyepiece-huygens-fe50.json: measure: lens: no ray reached the image plane
24. eyepiece-ramsden-fe15.json: measure: lens: no ray reached the image plane
25. eyepiece-ramsden-fe20.json: measure: lens: no ray reached the image plane
26. eyepiece-ramsden-fe25.json: measure: lens: no ray reached the image plane
27. eyepiece-ramsden-fe32.json: measure: lens: no ray reached the image plane
28. eyepiece-ramsden-fe40.json: measure: lens: no ray reached the image plane
29. eyepiece-abbe-fe6.json: measure: lens: no ray reached the image plane
30. eyepiece-abbe-fe10.json: measure: lens: no ray reached the image plane
31. eyepiece-abbe-fe12.5.json: measure: lens: no ray reached the image plane
32. eyepiece-abbe-fe15.json: measure: lens: no ray reached the image plane
33. eyepiece-abbe-fe20.json: measure: lens: no ray reached the image plane
34. eyepiece-abbe-fe25.json: measure: lens: no ray reached the image plane
35. eyepiece-huygens-pup5-fe20.json: measure: lens: no ray reached the image plane
36. eyepiece-huygens-pup5-fe32.json: measure: lens: no ray reached the image plane
37. eyepiece-abbe-short-fe9.json: measure: lens: no ray reached the image plane
38. eyepiece-abbe-short-fe16.json: measure: lens: no ray reached the image plane
39. micro-achr-40x.json: assertion: on-axis spot RMS 5.27607mm exceeds half the focal length; field spot RMS 5.30986mm exceeds the focal length
40. micro-achr-60x.json: assertion: on-axis spot RMS 31.8611mm exceeds half the focal length; field spot RMS 32.004mm exceeds the focal length
41. micro-lister-10x.json: assertion: only 19 of 37 sampled rays reach the image plane; only 15 of 37 sampled off-axis rays reach the image plane
42. micro-lister-20x.json: assertion: only 7 of 37 sampled rays reach the image plane; only 7 of 37 sampled off-axis rays reach the image plane
43. micro-compound-40x.json: measure: only 2 of 21 rays reach the last surface
44. mirror-fieldlens-f400.json: measure: only 0 of 21 rays reach the last surface
45. mirror-fieldlens-f800.json: measure: only 0 of 21 rays reach the last surface
46. relay-1to1-f50.json: measure: lens: no ray reached the image plane
47. relay-1to1-f100.json: measure: lens: no ray reached the image plane
48. relay-1to1-f200.json: measure: lens: no ray reached the image plane
49. telephoto-f150-ratio0p85.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
50. telephoto-f300-ratio0p75.json: assertion: on-axis spot RMS = 0 is not positive and finite; only 1 of 37 sampled rays reach the image plane; field spot RMS = 0 is not positive and finite; only 1 of 37 sampled off-axis rays reach the image plane
51. telephoto-f500-ratio0p7.json: measure: lens: no ray reached the image plane
52. achromat-asphere-s1-a4-f100-nsk16-nsf5.json: assertion: only 19 of 37 sampled rays reach the image plane; only 19 of 37 sampled off-axis rays reach the image plane

## 构造失败（求解/布局/质量门槛失败，未产生候选）

1. bend k=1.3: solveEFL: no sign change of the power residual for s in [1e-5.11, 1e0.89] (target -100 mm)
2. bend k=2: solveEFL: no sign change of the power residual for s in [1e-5.36, 1e0.64] (target -100 mm)
3. micro lister 40x: EFL residual 1.37e-06 after freezing (target 4, got 4.000005460802599)
4. reducer galilean 0.5x: solveAfocal: no sign change on [-30,210]
5. reducer galilean 0.25x: solveAfocal: no sign change on [-45,175]
6. reducer galilean 0.125x: solveAfocal: no sign change on [-52.5,157.5]
7. relay1to1 f=50: solveEFL: no sign change of the power residual for s in [1e-1.52, 1e4.48] (target 50 mm)
8. relay1to1 f=100: solveEFL: no sign change of the power residual for s in [1e-1.52, 1e4.48] (target 100 mm)
9. relay1to1 f=200: solveEFL: no sign change of the power residual for s in [1e-1.52, 1e4.48] (target 200 mm)
10. achromat asphere N-BK7/F2 f=100 F/2.8: marginal height: geo: ray escape at surface 1
11. achromat asphere N-SK16/N-SF5 f=50 F/5: marginal height: geo: ray escape at surface 1
12. achromat asphere N-SK16/N-SF5 f=150 F/5: marginal height: geo: ray escape at surface 1
13. achromat asphere N-SK16/N-SF5 f=100 F/2.8: marginal height: geo: ray escape at surface 1

## 合并（参数完全相同，保留先出现的条目）

1. mirror-newtonian-f1500-f5.json 与 mirror-parabola-f1500-f4.json 参数指纹相同（同一光学设计），合并

## 同名改名（参数不同，文件名加 -altN 后缀）

（无）


## 本批未覆盖（留给后续批次）

- 自由曲面（XY 多项式离轴面，缺口 G3）与泽尼克面（G4）：依赖尚未合并的面型能力，本批未生成（若某例需要该能力则跳过并记录于此）
- 施密特校正板与施密特-卡塞格林的解析 A4 校正板：需要校正板设计解，本次未做
- 里奇-克雷蒂安（双曲面，两镜圆锥常数联立求解）与格利高里式：RC 需要圆锥求解；格利高里式的副镜共轭符号未核实
- DOE / 二元光学、热光学（dn/dT、热膨胀）、偏振（琼斯矩阵）、杂散光（二次反射）：本批未涉及
- 卡塞格林的中心遮挡、牛顿式的完整 45° 折轴镜筒：本批只做单折轴镜折转与折叠臂透镜的折转光路
- 公差进阶（灵敏度、逆灵敏度、补偿器）、膜系联合评价：本批未涉及
