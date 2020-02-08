# IAST
Ideal Adsorbed Solution Theory

Python code for computing adsorption equilibrium of mixtures in porous materials. Currently, only the analytical integration of the Langmuir isotherm is used in calculating the spreading pressure.

Installation
============
```bash
git clone https://github.com/Heinen1/iast.git 
```

Usage
=====
```
python main.py 
```

In main.py two filename have to specified (file_name1 and file_name2) that refer to the single component adsorption isotherms that are calculated using the default RASPA perl scripts. In code, the first 22 lines of these first are being skipped, and column 1 and 8 correspond to the pressure [Pa] and absolute loading [mol/kg], respectively.


