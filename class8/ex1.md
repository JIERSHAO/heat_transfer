对于微元，仅考虑扩散项控制方程： ${\rho c}\frac{\partial T}{\partial t} = \frac{\partial }{\partial x}\left( {\lambda \frac{\partial T}{\partial x}}\right)$

---
对于第$i$个微元，依据惯例标记，时空积分，有：
$${\int }_{x }{\int }_{t}\rho c\frac{\partial T}{\partial x}{dt}{dx} = {\int }_{t}{\int }_{x}\frac{\partial }{\partial x}\left( {\lambda \frac{\partial T}{\partial x}}\right) {dxdt }$$ $${\int }_{x}\left( {{\left( \rho cT\right) }_{i}^{n + 1} - {\left( \rho cT\right) }_{i}^{n}}\right) {dx} = {\int }_{t}\left( {\lambda \frac{\partial T}{\partial x}}\right) {|}_{w}^{e}{dt}$$

阶梯化：
$$\left( {{\left( \rho cT\right) }_{i}^{n + 1} - {\left( \rho cT\right) }_{i}^{n}}\right) {\Delta x} = {\Delta t} \cdot  \left( {{\lambda }_{i, e}\frac{\partial {T}_{e}}{\partial {x}_{e}} - {\lambda }_{i, w}\frac{\partial {T}_{w}}{\partial {x}_{w}}}\right)$$

中心差分：
$$\left( {{\left( \rho cT\right) }_{i}^{n + 1} - {\left( \rho cT\right) }_{i}^{n}}\right) {\Delta x} = {\Delta t} \cdot  \left( {{\lambda }_{i, x}\frac{{T}_{i + 1} - {T}_{i}}{\Delta x} - {\lambda }_{i, w}\frac{{T}_{i} - {T}_{i - 1}}{\Delta x}}\right)$$

---
对任意空间上，取惯常标记，求和：
$$\mathop{\sum }\limits_{{i = 1}}^{N}\left( {{\left( \rho cT\right) }_{i}^{n + 1} - {\left( \rho cT\right) }_{i}^{n}}\right) {\Delta x} = \mathop{\sum }\limits_{{i = 1}}^{N}{\Delta t}\left( {{\lambda }_{i,e}\frac{T_{i + 1} - T_{i}}{\Delta x} - {\lambda }_{i,w}\frac{T_{i} - T_{i-1}}{\Delta x}}\right)$$
因为 ${\lambda }_{i,e}\frac{{T}_{i + 1} - {T}_{i}}{\Delta x} = {\lambda }_{i+1,w}\frac{{T}_{i + 1} - {T}_{i}}{\Delta x}$，有：
$${\left( \rho cT\right) }_{L}^{n + 1} - {\left( \rho cT\right) }_{L}^{n} = {\Delta t}\left( {{\lambda }_{N,e}\frac{T_{i + 1} - T_{i}}{\Delta x} - {\lambda }_{1,w}\frac{T_{i} - T_{i-1}}{\Delta x}}\right)$$
$${\left( \rho cT\right) }_{L}^{n + 1} - {\left( \rho cT\right) }_{L}^{n} = {\Delta t}\left( {{\lambda }_{out}\frac{\partial {T}_{out}}{\partial x} - {\lambda }_{in}\frac{\partial {T}_{in}}{\partial x}}\right)
$$

