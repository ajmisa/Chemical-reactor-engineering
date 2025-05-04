# Transesterification in a Rotor–Stator Spinning Disc Reactor (RS-SDR)

This project models the conversion of waste frying oil to fatty acid methyl esters (FAMEs) via transesterification in a rotor–stator spinning disc reactor (RS-SDR). The work was completed as part of the *Advanced Reactor Engineering* course at TU/e. The study investigates the effect of reactor design, mass transfer, and catalyst configuration on productivity, with comparisons to conventional stirred-tank reactors.

## Project Overview

The RS-SDR offers intensified mass transfer and interfacial area, making it suitable for liquid–liquid and three-phase reactions. This work examines transesterification reactions using three different configurations:

### Case 1: Homogeneous Reaction System
- Catalyst: NaOH dissolved in methanol
- Reaction occurs only in the MeOH phase
- TG is the continuous phase; MeOH is dispersed
- Reactor modeled as a CSTR for both phases
- Reaction network includes reversible steps from TG to glycerol

### Case 2: Slurry-Based Heterogeneous System
- Catalyst: Na₂SiO₃ nanoparticles suspended in the TG phase
- Reaction occurs on solid surfaces within a three-phase system
- Agglomeration modeled based on energy dissipation and particle-particle interactions
- Impact of slurry concentration and turbulence analyzed

### Case 3: Coated Catalyst on Reactor Disc
- Catalyst: Na₂SiO₃ nanoparticle layer coated on the rotor
- Reaction assumed to occur on catalyst layer immersed in MeOH–glycerol phase
- Evaluation of coating thickness, porosity, tortuosity, and effective surface area
- Balance between productivity and catalyst utilization considered

## Topics Covered

- Reaction kinetics: 6-step reversible reaction network
- Mass transfer and diffusivity: Wilke–Chang correlation
- Droplet size and energy dissipation: Kolmogorov–Hinze model
- Reactor hydrodynamics: rotation speed, gap size, Reynolds number effects
- Productivity analysis under different design configurations
- Scale-up for national processing demand (based on 44 million kg/year oil input)

## Key Methods and Correlations

- Kinetic constants for TG, DG, MG ↔ FAME + glycerol
- Energy dissipation rate as a function of rotation speed
- Mass transfer modeled via Sherwood correlations
- Agglomeration modeled based on turbulent breakup forces
- Productivity optimization based on design and throughput constraints

## Repository Structure

- `homogeneous/`: Modeling of homogeneous reaction system
- `slurry/`: Slurry-based catalytic transesterification with Na₂SiO₃ particles
- `coated/`: Disc-coated catalyst system design and modeling
- `Advanced_Reactor_Engineering_Final_Assignment.pdf`: Full project report

## Course Information

Advanced Reactor Engineering – ACRE  
Eindhoven University of Technology – Q2 2025

## Author

Adam Jordani Misa  
MSc Chemical Engineering – TU/e  
Email: aj.misa@outlook.com
