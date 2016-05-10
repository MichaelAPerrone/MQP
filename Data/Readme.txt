This folder contains the data collected from simulations of waves propagating in a medium with time-varying properties.
They represent the use of different material parameters
They were done in a space of just 500 points, so there is some granularity and spurious behavior of the numerical method in the calculation of energy, which in part can be taken care of by running a larger simulation. It may also be possible to change the derivative stencil I am using for a more accurate or well-behaved derivative approximation. The overall correct expected behavior of the system is observed though.

1. DefaultPlots are settings used by Billy that are known to work
They have slight reflection/impedance mismatch and wave speeds of 0.6 and 1.1 in the two materials

2. SanityCheck are settings without impedance mismatch or variable wave speed.
We should expect no change in energy besides small numerical deviation, and that is exactly what we observe.

3. NoImpedanceMismatch has no mismatch in impedance. This means waves will not reflect at spatial boundaries.
Surprisingly it has much the same curve as the default case, which just has small impedance mismatch.
This suggests the system is stable against reflection.

4. LargerImpedanceMismatch was especially surprising: There is still a large amount of energy added to the system, suggesting the system can be quite stable against reflection.

5. In ImpedanceAndWaveSpeedDifferrence I used a slower wave speed for the slow material and a faster wave speed for the fast material, in hopes that it would amplify the energy further and stabilize it against reflection, effectively by trapping the waves in the slow material until the switch in material properties occurs. It displayed some strange behavior in the energy measurement: The wave initially propagates as expected and gets stuck in the slow regions and then amplified, but once the frequency is too high, the low resolution mesh I use loses a lot of the high frequencies whenever they reach the slow material, because the mesh becomes larger than the frequencies. I bet this simulation would work better with a much finer mesh.

6. LargestImpedanceMismatch inputs even more energy overall, which after a point seems suspect, because we would expect that waves trapped in the faster material until the temporal switch would lose a significant amount of the energy. Perhaps there is a bug in the code. Moreover it is strange to see energy changig between a switch in material parameters in all of these simulations.
