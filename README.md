# LBM MATLAB
## Latviski
Šeit ir atrodami 8 MATLAB kodi režģa Bolcmaņa metodes lietojumiem plūsmas un konvekcijas-difūzijas-reakcijas problēmu risināšanai.
Visi kodi ir pārbaudīti uz MATLAB versijas R2019a.

- [Teslas ventiļa simulācija](BGK_tesla_final.m)
- [Siltumvadīšanas vienādojuma risināšana](Heat_transfer.m)
- [Degšanas procesa simulācija](BGK_combustion.m)
- [Siltas masas celšanās simulācija](BGK_buoyancy.m)
- Cilindra aptecēšanas simulācija
    - [Naivā implementācija](benchmark_naive.m)
    - [Optimizētā implementācija nr. 1](benchmark_cpu.m)
    - [Optimizētā implementācija nr. 2](cpu_optimized_v2.m)
    - [`gpuArray` implementācija](benchmark_gpu.m)

## English
Here You can find 8 MATLAB codes for lattice Boltzmann method applications in solving flow and convection-diffusion-reaction problems.
All codes are tested on MATLAB version R2019a.

- [Tesla valve simulation](BGK_tesla_final.m)
- [Heat transfer equation solving](Heat_transfer.m)
- [Simulation of a combustion process](BGK_combustion.m)
- [Simulation of buoyancy caused by warm mass](BGK_buoyancy.m)
- Flow around a cylinder
    - [Naive implementation](benchmark_naive.m)
    - [Optimized implementation no. 1](benchmark_cpu.m)
    - [Optimized implementation no. 2](cpu_optimized_v2.m)
    - [`gpuArray` implementation](benchmark_gpu.m)