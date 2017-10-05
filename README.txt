Version 2 History:

10-3-2017: Changes made
1. Removed Turbomachinery mode and Cp, ds/dx Cartesian plot.
2. Created new script to evaluate velocity values at hub, tip and mean. Run the
script in CFView and generate values manually.
3. Completed loss modelling correlations for Rotor 2.
4. Using the same correlations to compute losses across Rotor 1 and Stator 1.
5. Evaluated Diffusion factors for all stages.

10-4-2017: Changes made
1. Added Flow Coefficient of compressor and Vt to output into text file.
2. Moved velocity data to beginning.
3. Removed loop calculating for diffusion factor.
4. Incidence Loss is considered to be included in skin friction loss term.
5. Removed errors in calculating absolute flow angle at diffuser exit.
6. Added macro for single-stage baseline model.
