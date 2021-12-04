## Project5 - Double slit experiment
Sara Pernille Jensen, Alec Elías Sigurðarson, Håkon Olav Torvik

This repository contains our code for Project 5 of FYS4150 fall 2021.

# Double slit experiment
We have modelled the double slit experiment by solving the 2D Schrödinger equation for an initial wave packet and constant potential. The PDE was discretized using the Crank Nicolson method, such that the time evolution can be found by solving Au(t+dt) = Bu(t) where A and B are matrices, and u(t) a vector representing the quantum state at time t.

# Code
The code was implemented in Julia. Python and/or Julia was used for plotting. See [julialang.org](julialang.org) or [github.com/JuliaLang/julia](github.com/JuliaLang/julia) for instructions for installing julia.

# TODO
- Skrive funksjonen som lager potensialet
- endre simulate og main.jl slik at vi ikke trenger å lagre hele tilstanden, og heller kun det vi trenger for hver oppgave. Må da kanskje lage en stor, liknende funksjon for hver oppgave 🤷 Er ikke et problem slik det er nå, så dette er ikke veldig prioritert.
- Skrive plottekode
- optimalisere? Fiksa en bug, så nå funker det, men tar ~20min å kjøre for 7a.
