# Spin-1 TF Zitterbewegung — Simulation & Literature Notes

---

## 1. Numerical verification

### Bugs fixed in original notebook (`tangent_fermions_spin1_wavepacket.nb`)

**Bug 1 — `applyKernels` `ArrayReshape` (wrong Transpose)**
```mathematica
(* WRONG: *)  Transpose[ArrayReshape[psiK3, {3, nx*ny}]]
(* FIXED: *)  ArrayReshape[psiK3, {nx*ny, 3}]
```
Evidence: drift halved (20.5 vs 41 expected sites) with the wrong version.

**Bug 2 — Fourier convention**
```mathematica
(* WRONG: Fourier for real→k, InverseFourier for k→real *)
(* FIXED: InverseFourier for real→k, Fourier for k→real *)
```
Reason: Mathematica's `Fourier` uses $e^{+2\pi i k n/N}$ (positive sign = physics inverse DFT). Evidence: centroid drifted LEFT with correct $k_0$ before fix.

### Verified results (Wolfram Language Evaluator)

| Quantity | Theory | Measured | Ratio |
|---|---|---|---|
| Drift (E₊ eigenstate, 700 steps, dt=0.05) | 41.01 sites | 40.64 sites | 0.991 ✓ |
| ZB frequency — damped sine fit to ⟨Sx⟩ | 0.1318 c/t | 0.1300 c/t | 0.986 ✓ |
| ZB frequency — 8× zero-pad FFT | 0.1318 c/t | 0.1284 c/t | 0.974 |
| Final norm | 1.0 | 1.0000000 | — ✓ |

### Physical parameters at $k_0 = (\pi/4,\, 0)$

- TF dispersion: $E_{\rm TF}(k_0) = 2\tan(\pi/8) = 0.8284$ rad/t
- Group velocity: $v_g = \sec^2(\pi/8) = 1.1716$ sites/t
- ZB period: $T_{\rm ZB} = 2\pi/E_{\rm TF} = 7.584$ time-units
- ZB amplitude in $\langle x\rangle$: ~0.065 sites (small — see §2)
- ZB amplitude in $\langle S_x\rangle$: ~0.006
- ZB decay time: $\tau_{\rm fit} \approx 22$ time-units

### Initial state for ZB measurement

Equal superposition of E₊ and E_flat eigenstates of $S_x$:
```mathematica
psiPlus = N[(1/2){1., Sqrt[2]+0.I, 1.}];   (* Sx eigenvalue +1 *)
psiFlat = N[(1/Sqrt[2]){1., 0.+0.I, -1.}]; (* Sx eigenvalue  0 *)
spinor0ZB = Developer`ToPackedArray[(1/Sqrt[2])(psiPlus + psiFlat)];
```

### Code file

`ZB_verified_code.m` — complete Mathematica code (Sections 1 and 2), ready to paste into the notebook.

### SAFETY constraint

`ParallelDo` is **NEVER** used — it crashed the computer. All loops use `Do`.

---

## 2. Physics of Zitterbewegung — discussion notes

### Why two bands and not three?

The two-band picture is a **choice of initial state**, not a mathematical constraint. Three eigenstates exist — $E_+$, $E_{\rm flat} = 0$, $E_-$ — and a general initial state excites all three. The interference pairs and their ZB frequencies are:

| Pair | $\omega_{\rm ZB}$ |
|---|---|
| $E_+ \leftrightarrow E_{\rm flat}$ | $E_+ = 2\tan(k_{0x}/2)$ |
| $E_{\rm flat} \leftrightarrow E_-$ | $E_+$ (by symmetry) |
| $E_+ \leftrightarrow E_-$ | $2E_+$ |

With all three bands populated, the centroid oscillates at two incommensurable frequencies simultaneously. We chose the E₊ + E_flat superposition for three reasons:

1. **Single clean frequency** — one ZB tone, sharp FFT peak, clean damped-sine fit.
2. **Flat-band signature** — mixing E₊ with E_flat gives $\omega_{\rm ZB} = E_+$, exactly half the graphene value (E₊ + E₋ gives $2E_+$). This factor of ½ is a direct, measurable consequence of pseudospin-1 topology and has no analogue in graphene.
3. **Physical interpretability** — E_flat has zero group velocity everywhere (perfectly flat band), so the picture is clean: one component sprints, one is frozen, the centroid trembles between them.

### E₊ + E₋ ZB is graphene-like

In graphene (spin-½), ZB involves E₊ and E₋ — both components rush away from each other at $\pm v_F$, the centroid stays near $x_0$, oscillating at $\omega_{\rm ZB} = 2E$. The same physics reappears in pseudospin-1 if one prepares the E₊ + E₋ superposition: same qualitative picture, $\omega_{\rm ZB} = 2E_+$, nothing topologically new.

### The flat band changes the picture

The flat band has zero group velocity at **every** $k$-point. A wave packet built from flat-band states does not propagate. This has two consequences when mixing E₊ + E_flat:

**ZB frequency is halved.** $\omega_{\rm ZB} = E_+ - 0 = E_+$ instead of $2E_+$. The literature confirms this: in the $\alpha$-T₃ interpolation (graphene at $\alpha=0$, dice/pseudospin-1 at $\alpha=1$), the ZB frequency of a rim-site wave packet at $\alpha=1$ is exactly half the $\alpha=0$ (graphene) value [Bis17].

**ZB migrates from position to spin.** Because the flat band carries no current — the off-diagonal velocity matrix element $\langle E_{\rm flat}\lvert v_x \rvert E_+\rangle$ vanishes at $k_y = 0$ — the ZB oscillation in $\langle x\rangle$ is nearly zero (0.065 sites in our simulation, arising only from off-axis $k_y$ contributions of the Gaussian wave packet). The trembling survives instead in the internal spin degree of freedom: $\langle S_x\rangle(t)$ oscillates cleanly at $\omega_{\rm ZB}$, while the position channel is almost frozen. This is absent in the graphene case, where both components propagate and position and spin tremble together.

### Summary: what is new vs graphene

| Feature | Graphene (spin-½) | Pseudospin-1 (E₊+E₋) | Pseudospin-1 (E₊+E_flat) |
|---|---|---|---|
| $\omega_{\rm ZB}$ | $2E$ | $2E_+$ | $E_+$ — **halved** |
| Both components move? | Yes | Yes | No — E_flat frozen |
| ZB in $\langle x\rangle$ | Large | Large | **Nearly zero** |
| ZB in spin | Yes | Yes | **Dominant channel** |
| Flat-band role | None | None | **Central** |

---

## 3. Bibliographic research

*(Four deep searches completed in Undermind workspace "Tangent Fermions Spin-1 — Literature".)*

### 3a. Foundational references — pseudospin-1 / dice / Lieb lattice

**[Gol11]** N. Goldman, D. F. Urban, D. Bercioux — *Topological Phases for Fermionic Cold Atoms on the Lieb Lattice*, PRA 83, 063601 (2011).  
The key reference establishing the low-energy spin-1 Dirac Hamiltonian on the Lieb lattice, the single Dirac cone at the M point, and the Kane–Mele topological phases. Direct ancestor of this project.

**[Ber09]** D. Bercioux et al. — *Massless Dirac-Weyl Fermions in a T3 Optical Lattice*, PRA (2009).  
Original T₃ optical-lattice proposal establishing the Dirac-cone-plus-flat-band spectrum for pseudospin-1.

**[Dor11]** B. D'ora, J. Kailasvuori, R. Moessner — *Lattice generalization of the Dirac equation to general spin and the role of the flat band*, PRB 195422 (2011).  
Systematic spin-S generalization; clarifies how the flat band modifies transport and response in integer-pseudospin systems.

**[Urb11b]** D. F. Urban et al. — *Barrier transmission of Dirac-like pseudospin-one particles*, PRB 115136 (2011).  
Canonical prediction of all-angle super-Klein tunneling in the dice lattice.

### 3b. Tangent fermions

**[Sta82]** R. Stacey — *Eliminating lattice fermion doubling*, PRD (1982). Original construction.

**[Pac21]** M. J. Pacholski, C. Groth, C. Beenakker — *Generalized eigenproblem without fermion doubling for Dirac fermions on a lattice* (2021). Practical modern formulation (generalized eigenproblem $H\psi = E P\psi$) — the framework we use directly.

**[Bee23]** C. Beenakker et al. — *Tangent Fermions: Dirac or Majorana Fermions on a Lattice Without Fermion Doubling*, Annalen der Physik (2023). Most complete recent reference.

**[Vel22]** A. Vela et al. — *Massless Dirac Fermions on a Space-Time Lattice with a Topologically Protected Dirac Cone*, Annalen der Physik (2022). Topological protection of the Dirac cone in the TF discretization.

### 3c. ZB in pseudospin-1 / α-T₃ — key references

**[Bis17]** T. Biswas, T. K. Ghosh — *Dynamics of a quasiparticle in the α-T₃ model: role of pseudospin polarization and transverse magnetic field on Zitterbewegung*, J. Phys.: Condens. Matter (2017). **The central reference for ZB in pseudospin-1.**
Key findings directly relevant to our work:
- Rim-site (off-hub) initial state at $\alpha=1$ (dice/pseudospin-1): single ZB frequency, equal to **half** the $\alpha=0$ (graphene) frequency. Confirms our $\omega_{\rm ZB} = E_+$ instead of $2E_+$.
- Intermediate $0 < \alpha < 1$: two distinct ZB frequencies (beating). Three-band multi-frequency ZB appears generically.
- Hub-site initial state: single ZB frequency throughout the $\alpha$ interpolation.
- With transverse magnetic field: persistent (non-decaying) oscillations, more complex multi-frequency structure.
- Without field (our case): ZB is transient, decaying on a timescale $\tau \sim \sigma/\delta v$ where $\sigma$ is the packet width and $\delta v$ the velocity difference between bands.

**[Dav09]** G. Dávid, J. Cserti — *General theory of Zitterbewegung*, PRB 121417 (2009).  
Establishes that ZB is generically multifrequency for momentum-dependent multiband Hamiltonians. Framework for interpreting the three-band case.

**[Ori22]** D. O. Oriekhov, V. Gusynin — *Optical conductivity of semi-Dirac and pseudospin-1 models: Zitterbewegung approach*, PRB (2022).  
Computes optical conductivity via velocity correlators (ZB approach) for gapped dice and Lieb models. Discusses flat-band contribution thresholds. Closest paper to ZB in the Lieb lattice specifically.

**[Cse10]** J. Cserti, G. Dávid — *Relation between Zitterbewegung and the charge conductivity, Berry curvature, and the Chern number of multiband systems*, PRB 201405 (2010).  
Connects ZB velocity correlators to topological invariants. Relevant for understanding why ZB changes at a topological phase transition.

**[She22]** X. Shen, Y.-Q. Zhu, Z. Li — *Link between Zitterbewegung and topological phase transitions*, PRB (2022).  
ZB as a probe of topological phase transitions in spin-J systems.

**[Zaw11]** W. Zawadzki, T. Rusin — *Zitterbewegung (trembling motion) of electrons in semiconductors: a review*, J. Phys.: Condens. Matter (2011).  
Comprehensive review; good background reference.

### 3d. Lieb lattice wave-packet dynamics / conical diffraction

**[Die16]** F. Diebel et al. — *Conical Diffraction and Composite Lieb Bosons in Photonic Lattices*, PRL (2016). Photonic experiment observing pseudospin-dependent conical diffraction. Closest experimental analogue to ZB dynamics in the Lieb lattice.

**[Ley12]** D. Leykam, O. Bahat-Treidel, A. Desyatnikov — *Pseudospin and nonlinear conical diffraction in Lieb lattices*, PRA 031805 (2012). Theory of pseudospin-1 conical diffraction; distinct wave-packet patterns from different pseudospin initial states.

**[Tai15]** S. Taie et al. — *Coherent driving and freezing of bosonic matter wave in an optical Lieb lattice*, Science Advances (2015). Cold-atom realization; flat-band localization observed.

**[Muk14]** S. Mukherjee et al. — *Observation of a Localized Flat-Band State in a Photonic Lieb Lattice*, PRL (2014). First photonic observation of the flat-band localized state.

**[Vic14]** R. A. Vicencio et al. — *Observation of Localized States in Lieb Photonic Lattices*, PRL (2014). Companion photonic flat-band observation.

### 3e. Open question — gap in the literature

None of the Undermind searches found a paper that directly simulates the **centroid ZB oscillation** $\langle x\rangle(t)$ for a wave packet in the Lieb lattice (or its continuum pseudospin-1 limit) on a real-space lattice. The photonic experiments observe conical diffraction (the spatial intensity pattern), not centroid dynamics. [Ori22] computes the optical conductivity via ZB velocity correlators but not wave-packet trajectories. **This is the calculation performed in this project and represents a gap in the existing literature.**

---

## 4. Super-Klein tunneling simulation notes

### Optimal parameters (verified T ≈ 0.994)

```mathematica
Nx = 240; Ny = 160;
sigma = 16.0;
x0 = Round[0.20 Nx];   (* = 48, leaves 3σ clearance to barrier *)
y0 = Round[0.50 Ny];
xBarL = Round[0.45 Nx];  (* = 108 *)
xBarR = Round[0.55 Nx];  (* = 132 *)
k0x = 2 ArcTan[0.25];   (* E0 = 0.5 = V0/2 exactly → super-Klein regime *)
k0y = 0.0;
V0val = 1.0;
dt = 0.2; Nsteps = 550;
```

**Result:** T = 0.99421 ✓

### Why T < 1 (T ≈ 0.984 for σ=8)

The 2D Gaussian wave packet has a spread in $k_y$. Components with $k_y \neq 0$ have
$E = \sqrt{k_{0x}^2 + k_y^2} > V_0/2$ and transmit with T < 1.
Analytical integral over the Gaussian $k_y$ distribution gives $\langle T\rangle \approx 0.986$, matching the simulation.

**Fix:** use larger σ (smaller $k_y$ spread). Scaling rule: $N_x \geq 10\sigma$, $x_0 \geq 3\sigma$ before the barrier.

### Grid scaling rule

| σ | Nx (min) | Nsteps |
|---|---|---|
| 8 | 120 | 350 |
| 16 | 240 | 550 |
| 32 | 500 | 1100 |

For σ=32, wave packet center must be at least 96 sites before the barrier (3σ clearance).

---

*Last updated: October 2026*
