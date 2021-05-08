# lauesim
Simulate Laue backreflection diffraction patterns

Requires VESTA to calculate hkl planes. Download here: https://jp-minerals.org/vesta/en/download.html

To prepare .csv files containing hkl info that the simulation programs can read:
1. Open a .cif file in VESTA.
2. Go to: Utilities > Powder Diffraction Pattern.
3. In the tab 'Conditions', select only 1 wavelength.
4. Click 'Calculate'.
5. In the tab 'Reflections', check that all of the reflections you are interested in appear here. If not, you may need to decrease the value of the wavelength.
6. Once you are satisfied with the list of reflections, go to: File > Export Reflection Table, and save the .txt file.
7. Use the program 'makecsv.m' to convert this .txt file to a .csv file. At this stage, you can specify a minimum structure factor |F|.

*An example .txt and .csv file are included for reference

Use 'simulatepattern.m' to simulate a single pattern. Specify the parameters in the block labeled 'INPUT PARAMETERS', and run the script to perform the simulation.

You can also use the 'simulatepattern_pseudosymcompare.m' and 'simulatepattern_pseudosymcompare2.m' programs to compare the patterns for 2 or 3 pseudosymmetric variants, respectively, on top of each other.


Feel free to email me at epang@mit.edu if you have any questions/difficulties/suggestions.
