# Jdpd
Jdpd - An open Java Simulation Kernel for Molecular Fragment Dissipative Particle Dynamics

Jdpd is an open Java simulation kernel for Molecular Fragment Dissipative Particle Dynamics with parallelizable force calculation, efficient caching options and fast property calculations. It is characterized by an interface and factory-pattern driven design for simple code changes and may help to avoid problems of polyglot programming.
Detailed input/output communication, parallelization and process control as well as internal logging capabilities for debugging purposes are supported. The new kernel may be utilized in different simulation environments ranging from flexible scripting solutions up to fully integrated “all-in-one” simulation systems.

The Jdpd library uses the [Apache Commons RNG](http://commons.apache.org/proper/commons-rng/) libraries and a [PCG pseudorandom generator implementation for Java](https://github.com/alexeyr/pcg-java) and is published as open source under the GNU General Public License version 3. This repository comprises the Java bytecode libraries (including the Apache Commons RNG libraries), the Javadoc HTML documentation and the Netbeans source code packages including Unit tests.

Jdpd is described in the scientific literature (see reference below), the [final manuscript](https://github.com/zielesny/Jdpd/blob/master/2018%20-%20van%20den%20Broek%20-%20Jdpd%20-%20Final%20Manucsript.pdf) is added as a PDF document to the repository.

# Reference
[K. van den Broek, H. Kuhn and A. Zielesny, _Jdpd - An open Java Simulation Kernel for Molecular Fragment Dissipative Particle Dynamics_, Journal of Cheminformatics (2018), 10:25](https://doi.org/10.1186/s13321-018-0278-7)

# Acknowledgements
The support of [GNWI - Gesellschaft für naturwissenschaftliche Informatik mbH](http://www.gnwi.de) is gratefully acknowledged.
