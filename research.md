---
layout: page
title: Research
---


### Molecular Dynamics Simulations Guided by X-ray Scattering Signal

<img style="float: right; max-width: 40%; min-width: 400px; width: 500px" src="./assets/images/LAO.png" />

**Keywords: Molecular Dynamics, CUDA C, GPU, X-ray Scattering**

For the Github repo, check <a href="https://github.com/darrenjhsu/XSNAMD" target="_blank">here</a>.

This is my main Ph.D. project where we try to guide the protein motions in molecular dynamics (MD) simulations using the data we collected, which is time-resolved X-ray solution scattering signal (TRXSS) of proteins. For example, we try to drive the 1LST structure to the 2LAO structure using only the TRXSS signal difference between the two structures.

The TRXSS signal reflects changes in the atomic distance distribution, which can be modeled with the Debye formula. However Debye formula is very expensive (proportional to number of scatterers squared) so I developed GPU code to do parallel computing, giving it a 12,000x boost in efficiency.

See our paper here: <a href="https://aip.scitation.org/doi/full/10.1063/5.0007158" target="_blank">Integrating solvation shell structure in experimentally driven molecular dynamics using x-ray solution scattering data</a>.



### Finding the Unfolded State of a Protein

<!--During my Ph.D. I have been working with Prof. Lin Chen at Northwestern University on tracking and characterizing disordered protein structures when they are subject to external excitation, such as light, temperature-jump, pH-jump, or chemical reduction. I did a large part of my research at the Advanced Photon Source at Argonne National Laboratory, where I utilized the BioCARS beamline for time-resolved X-ray solution scattering, 11-ID-D beamline for X-ray transient absorption spectroscopy, and DND-CAT for small-angle X-ray scattering.-->

<img style="float: left; max-width: 40%; min-width: 400px; width: 500px" src="./assets/images/BLA.png" />

**Keywords: Data-driven Simulation, X-ray Scattering, Molecular Dynamics**

With the above method established, I decided to examine an intermediate unfolded state (*molten globule*) of the alpha-lactalbumin.
A temperature jump is applied to the protein solution, causing the protein the unfold. We then use TRXSS to track its unfolding in real time.
The signal is then input to the MD simulation, generating a different set of structures (the cyan ones in the figure; the purple one is the crystal structure) that are compatible with the data, the kinetic model, and the traditional descriptions of the *molten globule* state.
This is then our best guess of what the protein looks like in this state, and the modeling allows us to see the structure and gain insight into the system.

See our paper here: <a href="https://aip.scitation.org/doi/full/10.1063/5.0039194" target="_blank">Unfolding bovine α-lactalbumin with T-jump: Characterizing disordered intermediates via time-resolved x-ray solution scattering and molecular dynamics simulations</a>.


### Side Ph.D. projects 

- Designed an aluminum nitride sample holder that will not be evaporated by the IR-laser we use to generate temperature jumps.
- Probed structural dynamics of platinum and copper complexes using X-ray free electron lasers (XFEL). Multiple papers are in prep.
- Self-assembly of biomaterials as a function of temperature and pH.
- Code contribution to <a href="https://github.com/dleshchev/pytrx" target="_blank">**pytrx**</a>, a package for processing time-resolved X-ray scattering data.

