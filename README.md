# OptaDOS Photoemission Tutorials

This repository contains a set of **Markdown tutorials** for the OptaDOS photoemission module using **CASTEP 23**.

The tutorials are organised as separate markdown files and are intended to be served via **GitHub Pages** as a small documentation site. [Main Page](https://wuppersaver.github.io/)

## Usage

- The links associated with each live Tutorial lead to the page version of the markdown file.
- Python plotting snippets are included inline; more complete example scripts will be added within a dedicated repo (planned).

---

## Core Tutorials

1. **Cu(100) Slab Construction and Geometry, 3-Step Photoemission Model, Photon Energy Sweep and EDCs, 1‑Step Photoemission Model**
   
   [Link](https://wuppersaver.github.io/tutorial_Cu.html)

   Status: Complete  

Contents:
   - Building a Cu(100) slab with pymatgen  
   - Choosing slab thickness and vacuum  
   - Determining `photo_slab_min` / `photo_slab_max` from atomic z‑coordinates
   - Required CASTEP settings (`spectral_task`, PDOS weights, optics)  
   - Required files for OptaDOS (bands, dome_bin, ome_bin, pdos_bin, out.cell)
   - Minimal and recommended `Cu.odi` parameters  
   - Running OptaDOS and interpreting QE and MTE  
   - Common pitfalls and checks
   - `task : photo_energy_sweep` and constraints  
   - QE(ω) and MTE(ω) curves  
   - Binding‑energy curves (`photo_output : bindenergy_curve`)
   - FEM inputs (`*.fem_bin` and energy_info)  
   - Switching `photo_model : 1step`  
   - Comparison between emission models

---

## Optional Tutorials

2. **Explicit Layer Definition (`photo_layers_tops`)**  

   Link
   
   Status: WIP
   
   - Using `photo_layers_tops` and `photo_slab_middle` for custom layer boundaries  
   - When to prefer explicit layer definitions

3. **MgO(6) on Ag(16) Interface System**  

   Link
   
   Status: WIP
   
   - Heterostructure geometry and CASTEP setup  
   - Layer‑dependent IMFPs with `photo_imfp_choice : layers`  
   - Distinguishing MgO vs Ag contributions to QE and MTE

4. **ARPES‑Style Maps (Energy–Momentum and Momentum Maps)**  

   Link
   
   Status: WIP
   
   - `photo_output : ekin_ptrans_map` and `const_bindenergy_p_map`  
   - Mapping OptaDOS output to ARPES observables

5. **Field Emission and Schottky Barrier Lowering**  

   Link
   
   Status: WIP
   
   - `photo_elec_field` and effective work function  
   - Impact on QE and threshold behaviour

6. **DS‑Like Photoemission Model**
    
    Link
   
    Status: WIP
   
    - `photo_model : ds_like_pe`  
    - Fast, DOS‑based estimates of QE and MTE

---


## Contributions and Issues

If you:

- Find mistakes,
- Have suggestions for clarifications, or
- Want to contribute additional examples (e.g. different surfaces or materials),

please open an issue or submit a pull request to this repository.
