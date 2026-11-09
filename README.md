# IAST

Ideal Adsorbed Solution Theory for binary adsorption equilibria using Langmuir isotherms.

## Quick start

```bash
python -m pip install -r requirements.txt
python main.py
```

This calculates mixture loadings for a binary system using the bundled example isotherm data in `src/Results.dat-ZIF-8-433K-methane` and `src/Results.dat-ZIF-8-433K-ethane`.

## CLI usage

```bash
python main.py --y1 0.5 --pressure-start 1e5 --pressure-end 2e6 --points 20
python main.py --file1 /path/to/adsorption-methane.dat --file2 /path/to/adsorption-ethane.dat --output results.csv
```

The project fits each single-component isotherm to a Langmuir expression and then solves the coupled IAST equations for the mixture loading.


