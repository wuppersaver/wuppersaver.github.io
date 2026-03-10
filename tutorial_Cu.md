---
layout: post
title: Photoemission Tutorial
---
# OptaDOS Photoemission Tutorial (Cu(100), CASTEP 23

---

## 0. Overview

This tutorial walks through the photoemission module of OptaDOS using a Cu(100) slab as the main example, with optional extensions. It is structured in layers:

1. Main Cu(100) tutorial  
   - Build and converge a Cu(100) slab with pymatgen.  
   - Run a CASTEP 23 spectral task that writes all files needed for OptaDOS photoemission.  
   - Run a single‑energy 3‑step photoemission calculation in OptaDOS (`task : photoemission`).  
   - Generalise to a photon energy sweep (`task : photo_energy_sweep`).  
   - Run a 1‑step model for the same system.  
   - Produce a binding‑energy curve (EDC) for comparison to ARPES‑style data.

2. Optional / advanced tutorials  
   - Layer definition with `photo_layers_tops`.  
   - Automatic slab z‑limit helper (computing `photo_slab_min` / `photo_slab_max`).  
   - MgO(6 layers) / Ag(16 layers) interface and layer‑dependent IMFPs.  
   - ARPES‑style maps (`ekin_ptrans_map`, `const_bindenergy_p_map`).  
   - Field emission and Schottky barrier lowering.  
   - DS‑like emission model.

Throughout, we will:

- Highlight minimal user‑set parameters vs defaults / optional tweaks.  
- Include questions a new user might ask and answer them inline.  
- Provide Python/Matplotlib snippets for plotting; more complete scripts are planned for GitHub (TODO links).

---

## 1. Building a Cu(100) slab with pymatgen

### 1.1. Why a slab and how thick?

Photoemission is a surface phenomenon, but the electronic structure deep inside the slab should be bulk‑like. A good slab must:

- Have enough layers that the central layers reproduce bulk PDOS.  
- Have sufficient vacuum to avoid spurious interaction between periodic images and to define a clean work function.

For Cu(100), this tutorial assumes:

- 16 layers of Cu (symmetric slab: 8 layers per side).  
- Vacuum thickness ≈ 30 Å along the surface normal (here taken as z).

Question: Can I use fewer layers?  
Answer: We can, but must check that the innermost layer PDOS matches bulk Cu reasonably. If not, QE and MTE will not be representative of a semi‑infinite surface.

---

### 1.2. Generating Cu(100) slab with pymatgen

Python context: we will use `pymatgen` to start from bulk fcc Cu and build a (100) slab.

- We assume to start from **bulk fcc Cu** with lattice parameter $a \approx 3.615\ \mathrm{\AA}$ (adjust to our preferred value or relaxed DFT value).
- The **(100) surface normal** will be taken along **Cartesian z** so that `photo_slab_min` and `photo_slab_max` refer to the $z$‑direction in OptaDOS.
- We want roughly **16 atomic layers** of Cu and about **30 Å of vacuum**.

```python
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator, generate_all_slabs

# 1. Define bulk Cu (fcc, 1 atom / cell)
a = 3.615  # approximate Cu lattice parameter in Angstrom
latt = Lattice.cubic(a)
bulk_cu = Structure(latt, ["Cu"], [[0, 0, 0]])

# 2. Build a (100) slab with 16 atomic layers and ~30 A vacuum
miller_index = (1, 0, 0)
min_slab_thickness = 16 * (a / 2.0)   # for fcc (100), planes spaced by a/2
min_vacuum_thickness = 30.0

slabgen = SlabGenerator(
    initial_structure=bulk_cu,
    miller_index=miller_index,
    min_slab_size=min_slab_thickness,
    min_vacuum_size=min_vacuum_thickness,
    center_slab=True,        # symmetric slab
    in_unit_planes=True,
)

slabs = slabgen.get_slabs()
slab = slabs[0]

# 3. Export to a CIF or directly to CASTEP.cell format
slab.to(fmt="cif", filename="Cu_100_slab.cif")
```

### 1.3 Choosing `photo_slab_min` and `photo_slab_max`

**Practical notes:**

OptaDOS needs to know the **vertical extent of the actual material slab** (excluding vacuum) to:

- Divide the slab into **layers (boxes)** along $z$.
- Correctly compute layer volumes and electron escape distances.

The key steps are:

1. Extract the **atomic $z$‑coordinates** from our final structure (`Cu-out.cell` or the CIF).
2. Determine the **minimum** and **maximum** atomic $z$:  
   $z_{\min}^{\text{atoms}},\ z_{\max}^{\text{atoms}}.$
3. Choose an **atomic radius** $r_{\mathrm{vdW}}$ for Cu (a simple working value is $1.4\ \mathrm{\AA}$).
4. Define the slab limits as:
   - `photo_slab_min`:
     $z_{\min}^{\text{slab}} = z_{\min}^{\text{atoms}} - r_{\mathrm{vdW}}$
   - `photo_slab_max`:
     $z_{\max}^{\text{slab}} = z_{\max}^{\text{atoms}} + r_{\mathrm{vdW}}$

This ensures that:

- All Cu atoms lie strictly **between** `photo_slab_min` and `photo_slab_max`.
- A small buffer beyond the outermost atomic positions is included, representing the “outer edge” of the material.

A minimal helper script to compute these from the slab CIF is:

```python
from pymatgen.io.cif import CifParser
import numpy as np

parser = CifParser("Cu_100_slab.cif")
structure = parser.get_structures()[0]

# Extract all z-coordinates (Cartesian)
z_coords = [site.coords[2] for site in structure]
z_min_atoms = min(z_coords)
z_max_atoms = max(z_coords)

# Choose an approximate Cu radius
r_vdw = 1.4  # Angstrom

photo_slab_min = z_min_atoms - r_vdw
photo_slab_max = z_max_atoms + r_vdw

print("z_min_atoms:", z_min_atoms)
print("z_max_atoms:", z_max_atoms)
print("Recommended photo_slab_min:", photo_slab_min)
print("Recommended photo_slab_max:", photo_slab_max)
```

We can then copy these numeric values into our `Cu.odi`:

```text
photo_slab_min : <value_from_script>
photo_slab_max : <value_from_script>
```

---

### 1.4. Convergence with layer number

Before trusting photoemission results, check:

- The **PDOS** of the **innermost Cu layer** (from CASTEP) and compare to bulk Cu PDOS (i.e. from a bulk calculation).

If the central layer PDOS is significantly different from bulk, increase the number of layers (e.g. 20, 24, …) and re‑check.

This is done during the **CASTEP stage** (see below), not in OptaDOS.

---

## 2. CASTEP 23: spectral run for photoemission

### 2.1. Required files for 3‑step model (Cu primitive slab)

For the 3‑step model (without explicit surface transmission probability), OptaDOS expects these files from CASTEP:

- `Cu.bands` – eigenvalues and k‑points.
- `Cu.dome_bin` – band gradients (needed for adaptive broadening, etc.).
- `Cu.ome_bin` – optical matrix elements.
- `Cu.pdos_bin` – atomic / orbital projection weights.
- `Cu-out.cell` – final cell and atomic positions.

**Not needed for this main tutorial:**

- `*.fem_bin` – only needed for 1‑step model.
- `*.tmprob_bin` – not produced in CASTEP 23.
- `*.gkgrid_bin` – not produced in CASTEP 23; used only in supercell / ARPES‑supercell tutorial.

---

### 2.2. CASTEP `.param` highlights

A sketch of the relevant parameters (simplified):

```text
task                 : spectral
cut_off_energy       : 500 eV              ! example
xc_functional        : PBE
spin_polarised       : false               ! Cu here assumed non-magnetic

# k-point grid: 15x15x1 or higher
kpoint_mp_grid       : 15 15 1
kpoint_mp_offset     : 0.0 0.0 0.0

# Spectral outputs
spectral_task        : optics              ! triggers spectral infrastructure
pdos_calculate_weights : true              ! needed for Cu.pdos_bin
```

> **Question:** *Is 15×15×1 enough?*  
> **Answer:** It’s a **minimum**. Photoemission quantities converge slowly with $k$‑point sampling, so use **as dense as our computing resources allow**: 15×15×1, 21×21×1, 27×27×1, etc.

---

### 2.3. Band count

CASTEP’s spectral task will typically **add conduction bands**. We usually do not need extra settings, but if we want to be safe, we can request ~50% extra bands:

```text
perc_extra_bands          : 50    ! 50% extra bands
```

> **Question:** *What happens if I have too few bands?*  
> **Answer:** The photoemission matrix elements required for transitions above the Fermi level may be missing, truncating QE at higher photon energies.

---

### 2.4. Run CASTEP & check outputs

Run:

```bash
castep Cu
```

Check:

- That the spectral run completed successfully.
- That the files `Cu.bands`, `Cu.dome_bin`, `Cu.ome_bin`, `Cu.pdos_bin`, `Cu-out.cell` exist.

If something failed, inspect CASTEP’s `.err` and `.castep` output, fix, and rerun before moving to OptaDOS.

---

## 3. Minimal OptaDOS `.odi` for single‑energy 3‑step model (Cu)

Create `Cu.odi` in the same directory as the CASTEP outputs.

### 3.1. Minimal parameters we must set

These are **not optional**; we must choose them explicitly.

```text
# --- Task and model ---
task                 : photoemission
photo_model          : 3step            # 3-step model, Bloch final states

# --- Photon energy ---
photo_photon_energy  : 5.0             # example: 5.0 eV

# --- Slab geometry (in Angstrom) ---
photo_slab_min       : <value>         # e.g. from z_min_atoms - r_vdw
photo_slab_max       : <value>         # e.g. from z_max_atoms + r_vdw

# --- Work function (in eV) ---
photo_work_function  : <W>             # obtained from potential, literature
                                       # or Fall-Binggeli-Baldereschi method

# --- IMFP (inelastic mean free path) ---
photo_imfp_choice    : const
photo_imfp_value     : <lambda_Cu>     # e.g. from IMFP literature, in Angstrom
```

**Typical values:**

- `photo_slab_min`, `photo_slab_max`: from Section 1.3 (use helper script).
- `photo_work_function`:
  - Must approximate the true work function $W$ of Cu(100).
  - Could be:
    - Extracted from the **electrostatic potential** in the slab (DFT).
    - Taken from **experiment** (literature).
    - Determined via the method of **Fall, Binggeli, Baldereschi**.
  - note that literature values are spread out depending on source and method, so it is difficult to define a "correct" value
- `photo_imfp_value`: choose a **Cu IMFP** at the relevant electron energy from literature (e.g. NIST, literature tables).

> **Question:** *What if I don’t know the exact IMFP?*  
> **Answer:** Use a reasonable literature value for electrons with kinetic energy similar to your experiment (a few eV above threshold). Later you can test sensitivity by varying `photo_imfp_value` by ±5–10%.

---

### 3.2. Recommended defaults and useful options

These can be left at default, but it’s better to make them explicit.

```text
# --- General OptaDOS options ---
broadening           : adaptive
efermi               : optados        # let OptaDOS recompute Fermi energy
energy_unit          : eV

# --- Use / ignore tmprob (surface transmission) ---
photo_use_tmprob     : false
# CASTEP 23 does not provide *.tmprob_bin.
# Setting this to false is a reasonable approximation:
# Expect QE values to differ by ~10% compared to a full transmission model.

# --- Temperature (for Fermi-Dirac occupations) ---
photo_temperature    : 298.0          # K, room temperature

# --- Binding energy Gaussian broadening ---
# Baseline suggestion: 1/2 k_B T
# At T = 298 K, 1/2 k_B T ≈ 0.0129 eV.
photo_bindenergy_broadening : 0.0129

# --- Angular windows (degrees) ---
photo_theta_min      : 0.0
photo_theta_max      : 90.0
photo_phi_min        : 0.0
photo_phi_max        : 90.0
# Currently only one quadrant [0,90°] is considered for phi

# --- Momentum binning for maps (used later) ---
photo_pmat_bin_width : 0.01          # 1/Å
# For quick runs: 0.01 1/Å is acceptable; for higher resolution 0.005;
# 0.001 is much more expensive.
```
> **Question:** *How does `photo_temperature` relate to `photo_bindenergy_broadening`?*  
> **Answer:**  
> - `photo_temperature` enters **Fermi–Dirac occupations** $n(i, \mathbf{k}, T, s)$ and emission probability broadening in
$$\Theta'(\mathbf{k}_\perp) =
\begin{cases}
1, &
\mathbf{k}_\perp(i,j,s) \ge 0 \\[4pt]
\delta\bigl(E_v - E(j,\mathbf{k},s)\bigr), &
\mathbf{k}_\perp(i,j,s) < 0
\end{cases}$$
> - `photo_bindenergy_broadening` is the **Gaussian energy broadening width** used when building energy histograms (e.g. EDCs). A good starting point is $\frac{1}{2}k_B T$, but it can also be increased to simulate additional broadening (instrumental resolution, phonons, disorder).



---

### 3.3. Running the 3‑step calculation

Run serial OptaDOS (adjust executable name to our build):

```bash
optados.x86_64 Cu
```

**Parallel run:**

- Always run parallel if possible, but **do not exceed** the number of k‑points with the number of MPI processes:

  ```bash
  mpirun -np 8 optados.mpi.x86_64 Cu
  ```

  if there are 105 k‑points this is fine; but we are running with 256 processes on 105 k‑points, it will break parallelisation with possibly unexpected behaviour.

> **Question:** *How do I debug crashes?*  
> **Answer:**  
> - Check **`Cu.opt_err`** first: missing files, format mismatches, etc.  
> - Then check the end of **`Cu.odo`** for more context. Often you’ll see messages about missing `*.bands`, `*.pdos_bin`, etc.

---

### 3.4. Interpreting the main output

Key places in `Cu.odo`:

1. **Fermi energy analysis**

   Look for:

   ```text
   +----------------------------- Fermi Energy Analysis --------------------+
   | From Adaptive broadening                                             |
   |   Spin Component : 1 occupation between  ...                         |
   |         Fermi energy (Adaptive broadening) :   ... eV                |
   +-----------------------------------------------------------------------+
   ```

   Check that:

   - The number of electrons is close to the expected integer.
   - The Fermi level is reasonable.

2. **Photoemission parameters block**

   Example (from manual, adapted):

   ```text
   +----------------------- PHOTOEMISSION PARAMETERS -----------------------+
   |  Photoemission Model                       :  3-Step Model             |
   |  Photoemission Final State                 :  Bloch State              |
   |  Photon Energy              (eV)           :  5.0000                   |
   |  Work Function              (eV)           :  4.50                     |
   |  Slab Max Z-Coord.          (Ang)          :  <value>                  |
   |  Slab Min Z-Coord.          (Ang)          :  <value>                  |
   |  IMFP Constant              (Ang)          :  <lambda_Cu>              |
   |  Bulk cutoff dist. (int. multiple of IMFP) :  10.0                     |
   |  Electric Field Strength    (V/Ang)        :  0.0000                   |
   |  Smearing Temperature       (K)            :  298.0                    |
   |  Transverse Momentum Scheme                :  crystal                  |
   |  Theta    - min -           (deg)          :  0.00                     |
   |  Theta    - max -           (deg)          :  90.00                    |
   |  Phi      - min -           (deg)          :  0.00                     |
   |  Phi      - max -           (deg)          :  90.00                    |
   |  Binding Energy Broad. Width (eV)          :  0.0129                   |
   +-----------------------------------------------------------------------+
   ```

3. **Layer ordering and volumes**

   We'll see a table mapping **atoms → layers → z‑coordinates**, and the number of explicit layers used in the calculation (typically half the slab, i.e. 8 layers for a 16‑layer symmetric slab).

4. **Bulk approximation info**

   The code builds a “bulk extension” of the slab by repeating the innermost layer up to a depth of `photo_bulk_cutoff * lambda`. Default is **10×IMFP**, which is usually more than plenty.

5. **IMFP values**

   - Within this tutorial the IMFP is assumed to be a constant, i.e. `photo_imfp_choice : const` with a single `photo_imfp_value`, $$\lambda_\mu^{\mathrm{eff}} = \lambda$$ for all layers. 
   - A second option is the `photo_imfp_choice : layers`. This option requires the user to supply a list of IMFPs for each explicitly calculated layer (i.e. normally # layers / 2). The effective IMFP for each layer is then calculated as:

   $$
   \lambda_\mu^{\mathrm{eff}} =
   \frac{\sum_k t_k\,\lambda_k}{\sum_k t_k}
   $$

   where $t_k$ are the thicknesses of segments of different materials the electron must traverse and $\lambda_k$ their respective IMFPs.  

6. **Final QE and MTE**

   At the end:

   ```text
   | Total Quantum Efficiency (electrons/photon):  ...                     |
   | Weighted Mean Transverse Energy (eV):          ...                     |
   +-----------------------------------------------------------------------+
   ```

---

## 4. Photon energy sweep: `task : photo_energy_sweep`

Once the single‑energy calculation works, we can switch to a sweep to obtain QE(ω) and MTE(ω) curves.

### 4.1. `.odi` changes

Replace:

```text
task                : photoemission
photo_photon_energy : 5.0
```

by:

```text
task                : photo_energy_sweep
photo_photon_min    : 4.0      # eV, example
photo_photon_max    : 6.0      # eV, example
jdos_spacing        : 0.05     # eV; must divide (max - min) exactly
```

**Constraint:**
The difference between the min and the max must be an integer multiple of the `jdos_spacing`, or more formally
$$
\frac{\texttt{photo\_photon\_max} - \texttt{photo\_photon\_min}}{\texttt{jdos\_spacing}}
\in \mathbb{Z}
$$

For the example: $(6.0 - 4.0)/0.05 = 40$.

> **Question:** *Can I keep `photo_photon_energy` set?*  
> **Answer:** No. For `photo_energy_sweep`, **do not** set `photo_photon_energy`. If it is set, OptaDOS will error and stop.

All other parameters (slab geometry, work function, IMFP, etc.) remain as before.

---

### 4.2. Running and reading sweep outputs

Run:

```bash
optados.x86_64 Cu
```

The output will include:

- A set of QE and MTE values per photon energy step.
- Layer‑resolved contributions as before.

---

### 4.3. Plotting QE(ω) and MTE(ω) with Python

**Note:** These code snippets are **illustrative**; full plotting utilities will be uploaded to GitHub (TODO: link like `https://github.com/<user>/optados-photoemission-examples`). The first part exports the relevant values to a file, the python snippet then plots it.

```bash
awk '
  /Photon Energy/ { photon = $(NF-1) }
  /Total Quantum Efficiency \(electrons\/photon\):/ { tqe = $(NF-1) }
  /Weighted Mean Transverse Energy \(eV\):/ { wmte = $(NF-1) }
  END { printf "%s,%s,%s\n", photon, tqe, wmte }
' Cu.odo > Cu_photo_sweep.dat`
```
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Cu_photo_sweep.dat")
photon_energy = data[:, 0]   # eV
QE = data[:, 1]              # electrons / photon
MTE = data[:, 2]             # eV

fig, ax1 = plt.subplots()

ax1.set_xlabel("Photon energy (eV)")
ax1.set_ylabel("QE (electrons / photon)")
ax1.semilogy(photon_energy, QE, 'b-')   # log scale for QE
ax1.grid(True, which='both', axis='y')

ax2 = ax1.twinx()
ax2.set_ylabel("MTE (eV)")
ax2.plot(photon_energy, MTE, 'r-')

plt.title("Cu(100) photoemission: QE and MTE vs photon energy")
plt.tight_layout()
plt.show()
```

**Comparing to experiment:**

- Often, experimental QE is plotted against **excess energy**:
  $$
  E_{excess} = \hbar\omega - W
  $$
- If the DFT work function differs from experiment, using the excess energy as the energy coordinate will hide the work function differences.

---

## 5. 1‑step model for Cu(100)

The 1‑step model treats the final state as a **free electron with an exponential decay** inside the solid before escape.

### 5.1. Additional input from CASTEP: FEM

We must have `Cu.fem_bin` containing free-electron matrix elements and the `energy_info` block:

- `fem_energy_info(1)` – number of photon energy steps (dimensionless).
- `fem_energy_info(2)` – minimum photon energy (eV).
- `fem_energy_info(3)` – photon energy step size (eV).
- `fem_energy_info(4)` – assumed Fermi energy (eV).
- `fem_energy_info(5)` – assumed work function (eV).

All energies are in **eV**.

> **Question:** *Why is this important?*  
> **Answer:** The FEM matrix elements are computed for a fixed set of photon energies with a specific Fermi energy and work function. All these parameters determine the electron's energy after excitation/emission. CASTEP must calculate the matrix elements with the values chosen for OptaDOS to ensure the transition probabilities match.

---

### 5.2. `.odi` changes for 1‑step

Starting from our working 3‑step `.odi`, change:

```text
photo_model       : 1step
photo_use_tmprob  : false           # tmprob is a 3-step concept; not needed here
```

Keep:

- `photo_slab_min/max`
- `photo_work_function`
- `photo_imfp_choice` / `photo_imfp_value`
- `photo_temperature`
- `photo_bindenergy_broadening`

Then run as before:

```bash
optados.x86_64 Cu
```

---

### 5.3. Comparing 3‑step vs 1‑step

We can:

- Run a **single photon energy** calculation with both models.
- Compare:

  - QE at that energy.
  - MTE at that energy.

Differences are expected due to different treatment of final states and electron transport.

---

## 6. Binding‑energy curve (EDC)

A binding‑energy curve (EDC) is often used in ARPES; OptaDOS can generate something analogous.

### 6.1. `.odi` option

Add:

```text
photo_output : bindenergy_curve
```

This will:

- Compute a **QE vs binding energy** histogram.
- Apply **Gaussian broadening in energy** using `photo_bindenergy_broadening`.
- Use bins of width $0.001\ \mathrm{eV}$ internally.

Angular selection is controlled by `photo_theta_min/max`, `photo_phi_min/max`. For example, to simulate a narrow angular acceptance:

```text
photo_theta_min : 0.0
photo_theta_max : 30.0
photo_phi_min   : 0.0
photo_phi_max   : 30.0
```

---

### 6.2. Plotting the EDC

Assuming the file is named `Cu_bindenergy_curve.dat`:

Columns (typical):

1. Binding energy $E_B$ relative to $E_F$ (eV).
2. QE intensity.

Python snippet:

```python
import numpy as np
import matplotlib.pyplot as plt

edc = np.loadtxt("Cu_bindenergy_curve.dat")
E_B = edc[:, 0]  # eV, binding energy (E - E_F)
I   = edc[:, 1]  # QE intensity

plt.figure()
plt.semilogy(E_B, I, 'k-')
plt.xlabel("Binding energy (eV, relative to $E_F$)")
plt.ylabel("QE intensity (arb. units)")
plt.title("Cu(100) EDC (photo_output: bindenergy_curve)")
plt.grid(True, which='both', axis='y')
plt.show()
```

Use:

- **Log scale** for intensity.
- Linear scale for energy.

---

## 7. Common pitfalls and checklists

### 7.1. Geometry and slab

- **Wrong `photo_slab_min/max`**:
  - Too small → clipping atoms.
  - Too large → including vacuum as material.
- **Too few layers**:
  - Central PDOS not bulk‑like; check with CASTEP PDOS.

### 7.2. Missing or inconsistent files

- Forgot `PDOS_CALCULATE_WEIGHTS : TRUE` → no `*.pdos_bin`.
- No `*.ome_bin` or `*.dome_bin` → photoemission tasks cannot run.
- Incorrect or missing `*.fem_bin` if using 1‑step, the OptaDOS will check for file incompatibility.

Always check `Cu.opt_err` first; it often tells us exactly what file is missing or incompatible.

### 7.3. Numerical choices

- `jdos_spacing`:
  - Too small → longer runtime and memory.
  - For energy sweeps and optics, start around $0.01\text{–}0.05\ \mathrm{eV}$ keeping in mind the integer number of steps for `photo_energy_sweep`.
- `photo_pmat_bin_width`:
  - Too fine (e.g. $0.001\ \mathrm{\AA^{-1}}$) greatly increases computation time. Start with $0.01\ \mathrm{\AA^{-1}}$ to get a general picture.
- `efermi`:
  - OptaDOS has good algorithms to calculate $E_F$. Prefer to use `efermi : optados` unless there is have a good reason, like slight misalignment in the 1-step model input.

### 7.4. Parallelism

- **Do not exceed** the number of k‑points with the number of MPI processes.
- Use as many cores as memory allows and as is reasonable for the machine.
