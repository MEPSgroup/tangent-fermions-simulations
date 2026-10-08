# Tangent Fermions — Spin-1 Simulations

Mathematica notebooks for simulating pseudospin-1 (Lieb lattice / dice lattice) fermions discretized using the **tangent fermion** (TF) scheme, which eliminates lattice fermion doubling.

## Repository structure

```
tangent-fermions-simulations/
├── super-Klein-tunneling/
│   └── tangent_fermions_spin1_wavepacket-1.nb   # barrier tunneling simulation
├── Zitterbewegung/
│   └── tangent_fermions_spin1_wavepacket.nb      # ZB trembling-motion simulation
└── notes/
    └── ZB_simulation_status.md                   # physics notes, results, bibliography
```

## Physics background

### Tangent fermion discretization

The tangent fermion scheme replaces the continuum momentum $k_\alpha$ with the lattice approximant

$$k_\alpha \;\to\; \frac{2}{a}\tan\!\left(\frac{a k_\alpha}{2}\right)$$

leading to a **generalized eigenproblem** $H_\mathrm{TF}\,\Psi = E\,P_\mathrm{TF}\,\Psi$ with no fermion doubling.
See [Pacholski et al. (2021)](https://doi.org/10.21468/SciPostPhys.11.5.105) and [Beenakker et al. (2023)](https://doi.org/10.1002/andp.202200506).

### Pseudospin-1 Hamiltonian

The low-energy Hamiltonian is $H = v_F(S_x k_x + S_y k_y)$, where $S_\alpha$ are the $3\times 3$ spin-1 matrices.
The spectrum has three bands: $E_+$, $E_\mathrm{flat}=0$, $E_-$.

### Super-Klein tunneling

At $E = V_0/2$, transmission through a rectangular barrier is $T=1$ for **all** angles and **all** barrier widths — a unique consequence of the flat band.
Optimal simulation parameters: `Nx=240, Ny=160, σ=16, k0x = 2 ArcTan[0.25]`, yielding $T \approx 0.994$.

### Zitterbewegung (ZB)

For an initial state $\frac{1}{\sqrt{2}}(|E_+\rangle + |E_\mathrm{flat}\rangle)$ the centroid oscillates at $\omega_\mathrm{ZB} = E_+ = 2\tan(k_{0x}/2)$, **half** the graphene value — a direct signature of the flat band.

## Safety note

`ParallelDo` is **never** used in these notebooks — it can crash Mathematica on some systems. All loops use `Do`.

## References

- Stacey (1982), PRD — original TF construction  
- Pacholski, Groth, Beenakker (2021) — generalized eigenproblem formulation  
- Beenakker et al. (2023), Annalen der Physik — TF review  
- Urban et al. (2011), PRB — super-Klein tunneling prediction  
- Biswas & Ghosh (2017), J. Phys.: CM — ZB in α-T₃ / pseudospin-1  
- Goldman, Urban, Bercioux (2011), PRA — Lieb lattice topology  

## Authors

D. Bercioux (DIPC) — dario.bercioux@dipc.org
