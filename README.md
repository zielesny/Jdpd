# Jdpd
Jdpd - An open Java Simulation Kernel for Molecular Fragment Dissipative Particle Dynamics (DPD)

Jdpd is an open Java simulation kernel for Molecular Fragment Dissipative Particle Dynamics (DPD) with parallelizable force calculation, efficient caching options and fast property calculations. It is characterized by an interface and factory-pattern driven design for simple code changes and may help to avoid problems of polyglot programming.
Detailed input/output communication, parallelization and process control as well as internal logging capabilities for debugging purposes are supported. The kernel may be utilized in different simulation environments ranging from flexible scripting solutions up to fully integrated “all-in-one” simulation systems like [MFsim](https://github.com/zielesny/MFsim).

Jdpd uses the [Apache Commons Math](https://commons.apache.org/proper/commons-math/) and [Apache Commons RNG](http://commons.apache.org/proper/commons-rng/) libraries and is published as open source under the GNU General Public License version 3. This repository comprises the Java bytecode libraries (including the Apache Commons Math and RNG libraries), the Javadoc HTML documentation and the Netbeans source code packages including Unit tests.

Jdpd has been described in the scientific literature (the final manuscript *2018 - van den Broek - Jdpd - Final Manucsript.pdf* is added to the repository) and used for a biomolecular DPD study of cyclotide/membrane interactions (see references below).

See text file *JdpdVersionHistory.txt* for a version history with detailed information.

# Important bug fixes
- Jdpd releases **prior to version 1.6.0.0** contained a severe error when restarting a job with continued calculation of the radius of gyration: This led to an illegal file output so that the radius of gyration data of the restarted job result could not be evaluated.
- Jdpd releases **prior to version 1.5.0.0** contained a severe error in the radius of gyration calculation method which could affect the position of particles in the simulation box and lead to wrong and unphysical "particle jumps".
- Jdpd releases **prior to version 1.3.0.0** contained a severe mapping error regarding electrostatics force calculations which could affect uncharged particles.

# References
- [K. van den Broek, H. Kuhn and A. Zielesny, _Jdpd - An open Java Simulation Kernel for Molecular Fragment Dissipative Particle Dynamics_, Journal of Cheminformatics (2018), 10:25](https://doi.org/10.1186/s13321-018-0278-7)
- [K. van den Broek, M. Epple, L. S. Kersten, H. Kuhn and A. Zielesny, _Quantitative Estimation of Cyclotide-Induced Bilayer Membrane Disruption by Lipid Extraction with Mesoscopic Simulation_, Journal of Chemical Information an Modeling (2021), 61, 3027-3040](https://doi.org/10.1021/acs.jcim.1c00332) ([Link to ChemRxiv preprint](https://doi.org/10.26434/chemrxiv.14135783.v1))

# Acknowledgements
The support of [GNWI - Gesellschaft für naturwissenschaftliche Informatik mbH](http://www.gnwi.de) is gratefully acknowledged.
