# protein-unfolding
modelling protein unfolding using monte carlo method
Monte Carlo simulation is a technique used to study how a model responds to randomly generated inputs. It typically involves a three-step process:
1. Randomly generate “N” inputs. (sometimes called scenarios).
2. Run the simulation for each of the “N” inputs. Simulations are run on a computerized model of the system being analysed.
3. Aggregate and assess the outputs from the simulations.

## Documentation of Code
Monte Carlo simulation is commonly used to compute several pathways in understanding thermodynamic mechanisms. Denaturation of protein or unfolding of proteins can be viewed analogously as a phase change problem from the thermodynamic point of view. A simulation run is a series of random steps in conformation space, each perturbing some degrees of freedom of the molecule. A step is accepted with a probability that depends on the change in the value of an energy function. Proteins are assumed to be two-dimensional structures in a lattice. The amino acids occupy the lattice points and the covalent amide bonds the lattice edge. Each run at a particular temperature ( here T = 1/Kb, where Kb is Boltzmann constant) consisted of 1000,000 steps. In the first step, a folded 16-mer protein with only those interactions present in the folded structure have favourable interaction energy- the energy of any other non-covalent interaction is considered zero.

In the subsequent steps, the structure was chosen by picking a particular link randomly and giving it a rotation (clockwise or anticlockwise) again randomly, if corner lattice points or else undergo a crankshaft move(2D transformation or rotation in 3D). In Corner move, the corner most amino acid is rotated clockwise or anticlockwise randomly with respect to the adjacent amino acid residue.
Corner move:

(x,y) -> new lattice point, (A(i),B(i)) -> corner most amino acid lattice point, (A(a),B(a)) -> adjacent amino acid lattice point.

x= (cos(theta)*(A(i)-A(a))- sin(theta)*(B(i)-B(a))+A(a))

y= (sin(theta)*(A(i)-A(a))- cos(theta)*(B(i)-B(a))+B(a))

Crankshaft move:

(x,y) -> new lattice point, (A(i),B(i)) -> any one of middle amino acids lattice point

x= A(i-1)+A(i+1)-A(i)

y= B(i-1)+B(i+1)-B(i)

the new structure was accepted after checking it for steric hindrance or overlapping of the
new point with the other lattice points of protein. Assuming the energy of the new and old
structures as E1 and E2 respectively, the probability condition parameter w is given by

W = exp(-(E2-E1)/Kb*T)

where Kb and T are the Boltzmann’s constant and the temperature respectively. The new the step was accepted for folding if W was greater than 1. If W was less than 1, the step was accepted only if a random number, generated from a uniform distribution between the interval 0 and 1, was greater than W. The energy at each step is calculated based on the interaction between only those pairs of interactions present in the native folded protein structure. Every interaction is assigned interaction energy of -1.5 or 1.

## Result and Discussion:
• As the energy of interaction becomes less negative the protein unfolds quickly as the activation energy for the next move decreases.

• Energy vs Monte Carlo steps (iteration) plot when noncovalent native bond interaction energy is -1.5 and -1:
Protein prefers low energy states than higher energy states for more negative energy of interaction that is -1.5 here. The graph is spread throughout the energy states.

While in the case of energy -1, the graph is more populated at higher energy and less at lower energy states

• As the non-covalent interaction energy increases it becomes less feasible to break the interaction hence the graph is less populated at higher energy states for higher energy of interaction i.e. -1.5

• Probability energy function given by W = exp(-(E2-E1)/Kb*T) decreases with increase in interaction energy. Implies the probability of acceptance of move for interaction energy -1.5 is less than for interaction energy -1.

• After achieving zero interaction energy state protein can again start to refold and can repeat the process for multiple times. Here we simulated folding/unfolding of protein for 1 million iterations.

## Conclusion
We are able to visualize a folded protein with given non-covalent bond interaction energy getting unfolded by picking any lattice point randomly and performing the feasible move. Accepting or rejecting the move on the basis of its total energy of interaction. This Monte Carlo simulation method of protein is based on Metropolis algorithm
