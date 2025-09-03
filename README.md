# 4th-Year-Project
# Filament Bundle Model of Optic Nerve Damage

This project investigates the mechanical instabilities of the **optic nerve bundle**, focusing on how progressive damage and recovery rules influence deformation and buckling under intraocular pressure (IOP). The aim is to better understand mechanisms relevant to **glaucoma-related vision loss**.

---

## Project Overview
- Built upon a **discrete elastic rod model** implemented in Fortran for energy minimisation.  
- Developed **Python modules** to:  
  - Implement **new lamina phase damage rules** (progressive weakening of bending stiffness).  
  - Add **recovery rules** that restore stiffness under specific curvature conditions.  
  - Create **visualisations** (heatmaps, tables, deformation plots) to illustrate damage and collapse patterns.  

---

## My Contributions
- Designed and implemented **new lamina phase damage rules** in Python.  
- Added **visualisation techniques** (heatmaps, statistical tables, graphs).  
- Analysed simulation results and documented findings in the project report.  

---

## Credits
- **Professor Chris Prior (Durham University)** — provided the base **Fortran code** for the discrete elastic rod model and the **shell script** to run the energy minimisation.  
- **Callum Patel (2025)** — authored the project report and implemented the Python components described above.  

---

## How to Run
1. Clone this repository:  
   ```bash
   git clone https://github.com/callumpatel1/NewGmin.git
2. Compile and run Fortran minimisation code using the provided shell script:
   ./make_runs.sh
3. Apply Python damage/recovery rules:
   python adjust_damage.py
4. Generate visualisations using the provided Python scripts.

## Visualisations
Python scripts generate:
- Heatmaps of bending stiffness distribution.
- Line plots of rod deformation under successive iterations.
- Tables showing stiffness degradation statistics.

##License
This repository is for academic and research purposes. Please contact the author before reuse.

