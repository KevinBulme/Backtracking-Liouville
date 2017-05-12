# Backtracking-Liouville
Development of an advanced test-particle simulation tool for high-resolution velocity distribution function in space plasma

In the research group of Prof. Usui and Associate Prof. Miyake in Kobe university, EMSES, electromagnetic spacecraft environment simulator, has been developed and applied to the simulation studies on spacecraft-plasma interactions and other various phenomena occurring in the space environment. In the simulation output, velocity distribution function in each region is very important because it enables us to evaluate the stability of plasma environment. Because of the limited number of plasma particles treated in the simulation, however, it is difficult to obtain the high-resolution function in the EMSES simulations.
In this project, we newly develop an advanced test-particle simulation tool to obtain high-resolution velocity distribution function at arbitrary point in the simulation domain. We plan to apply the back-tracking Liouville method (see the reference) to the test simulation in which particle trajectories are calculated using the fields obtained from in the EMSES simulation. This new simulation tool becomes a very powerful tool to identify the unstable region and source of plasma instabilities in the simulation domain.

Reference
Richard Marchand, Test-Particle Simulation of Space Plasmas, Commun. Comput. Phys., Vol 8, No.3, pp.471-483, 2010.
