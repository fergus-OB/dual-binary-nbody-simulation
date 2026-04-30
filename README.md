# Dual Binary N-Body Simulation

A MATLAB numerical simulation of the **HD 74438 2+2 quadruple star system**, modelling two interacting binary pairs using Newtonian gravitational dynamics, ODE integration, barycentric coordinates, orbital visualisation, and energy-conservation diagnostics.

This project was completed as part of a Financial Mathematics final-year dynamical systems project and has been cleaned for public portfolio use.

## Project Overview

HD 74438 is a hierarchical quadruple star system made up of two close binary pairs, AB and CD, which orbit each other. Because all four bodies interact gravitationally, the system provides a useful case study in orbital stability, secular dynamics, and chaotic behaviour.

The project simulates the system as a Newtonian four-body problem and explores whether the observed orbital structure can be reproduced numerically.

## Objectives

- Model HD 74438 as a four-body gravitational system.
- Generate initial conditions for two inner binaries and their wider outer orbit.
- Integrate the equations of motion using MATLAB ODE solvers.
- Validate simulated orbits against known semi-major axes, eccentricities, and orbital periods.
- Visualise inner binary motion, barycentric outer orbits, long-term rosette patterns, and homoclinic-like encounters.
- Monitor numerical reliability using relative energy-error diagnostics.

## Mathematical Model

The simulation treats each star as a point mass and integrates the Newtonian N-body equations:

```math
\ddot{\mathbf r}_i =
-G \sum_{j \ne i} m_j
\frac{\mathbf r_i - \mathbf r_j}{|\mathbf r_i - \mathbf r_j|^3},
\quad i = 1,\dots,4
```

## Repository Structure

```text
dual-binary-nbody-simulation/
├── src/
│   ├── hd74438_setup.m
│   ├── thesis.m
│   ├── thesis2.m
│   ├── Thesis3.m
│   └── rewrite.m
├── outputs/
│   ├── baseline-orbits/
│   ├── long-term-dynamics/
│   ├── chaos-diagnostics/
│   └── exploratory/
├── report/
│   └── dual_binary_systems_report.pdf
├── docs/
├── README.md
├── .gitignore
└── LICENSE
```
