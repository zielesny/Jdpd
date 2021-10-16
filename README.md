# Jdpd
Jdpd - An open Java Simulation Kernel for Molecular Fragment Dissipative Particle Dynamics (DPD)

Jdpd is an open Java simulation kernel for Molecular Fragment Dissipative Particle Dynamics (DPD) with parallelizable force calculation, efficient caching options and fast property calculations. It is characterized by an interface and factory-pattern driven design for simple code changes and may help to avoid problems of polyglot programming.
Detailed input/output communication, parallelization and process control as well as internal logging capabilities for debugging purposes are supported. The kernel may be utilized in different simulation environments ranging from flexible scripting solutions up to fully integrated “all-in-one” simulation systems like [MFsim](https://github.com/zielesny/MFsim).

Jdpd uses the [Apache Commons RNG](http://commons.apache.org/proper/commons-rng/) libraries and is published as open source under the GNU General Public License version 3. This repository comprises the Java bytecode libraries (including the Apache Commons RNG libraries), the Javadoc HTML documentation and the Netbeans source code packages including Unit tests.

Jdpd is described in the scientific literature (see reference below), the final manuscript (*2018 - van den Broek - Jdpd - Final Manucsript.pdf*) is added to the repository.

See text file *JdpdVersionHistory.txt* for a version history with detailed information.

# Important bug fix
Jdpd releases **prior to version 1.3.0.0** contained a severe mapping error regarding electrostatics force calculations and should be replaced by a later release (i.e. version 1.3.0.0 or higher).

# Reference
[K. van den Broek, H. Kuhn and A. Zielesny, _Jdpd - An open Java Simulation Kernel for Molecular Fragment Dissipative Particle Dynamics_, Journal of Cheminformatics (2018), 10:25](https://doi.org/10.1186/s13321-018-0278-7)

# Acknowledgements
The support of [GNWI - Gesellschaft für naturwissenschaftliche Informatik mbH](http://www.gnwi.de) is gratefully acknowledged.
