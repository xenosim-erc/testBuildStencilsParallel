

# Results for testBuildStencilsParallel_v2

- Looping is performed 300 times (6 integration points per face)
  - In each loop field D is exchanged between processors
- Cubic interpolation (p=3) is assumed:
  - 20 cells in stencil in 2D
  - 100 cells in stencil in 3D
- Stencils are build using layer technique
- Speed-up is calculated by dividing time with time in serial
- Total time - looping time = stencil construction
- Looping time is max looping time over all processors

## 2D case

<img src="run_Apple_M4_Pro_2026-01-14_23-06-54/strongScaling_hex_2D.pdf" width="600">

<div style="display:flex; flex-wrap:wrap; gap:20px; justify-content:center;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_400.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_1600.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_6400.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_25600.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_102400.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_409600.pdf" style="width:45%;">
</div>

# 3D case

<img src="run_Apple_M4_Pro_2026-01-14_23-06-54/strongScaling_hex_3D.pdf" width="600">

Same diagram with limited y axis:

<img src="run_Apple_M4_Pro_2026-01-14_23-06-54/strongScaling_hex_3D_limitedY.pdf" width="600">

<div style="display:flex; flex-wrap:wrap; gap:20px; justify-content:center;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_1000.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_4096.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_15625.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_64000.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_262144.pdf" style="width:45%;">
  <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_1000000.pdf" style="width:45%;">
</div>



### **Some conclusions:**

- Building the stencil using indexedOctree is slower. In 2D indexedOctree is extremely slow. For diagrams in paper i was using idexedOctree - according to these results we will improve timings a lot with layer technique.

  Timings for 2D case - (**left octree, right layer technique**)

  <div style="display:flex; flex-wrap:wrap; gap:20px; justify-content:center;">
    <img src="run_Apple_M4_Pro_2026-01-15_10-47-06/timings_hex_2D_25600.pdf" style="width:45%;">
    <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_2D_25600.pdf" style="width:45%;">
  </div>

  Timings for 3D case - (**left octree, right layer technique**)

  <div style="display:flex; flex-wrap:wrap; gap:20px; justify-content:center;">
    <img src="run_Apple_M4_Pro_2026-01-15_10-54-00/timings_hex_3D_15625.pdf" style="width:45%;">
    <img src="run_Apple_M4_Pro_2026-01-14_23-06-54/timings_hex_3D_15625.pdf" style="width:45%;">
  </div>

- I am observing some unexplained anomalies in 3D, with extreme speed-ups or slow-downs.

