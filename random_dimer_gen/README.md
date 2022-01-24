# Random Dimer Structure Generation

An end-to-end example of generating a set of dimer structures from the example methyl viologen.

## Files

* ``random_dimer_gen.py``: the main file that generates the random dimer structures.
* ``monomer.xyz``: initial geometry passed for dimer genertaion.
* ``random_pairs``: a folder includes each randomly generated dimer structure geometry_{index}.xyz along with the rotataion matrix and translation vector that generate the dimer.

## Usage

```
python random_dimer_gen.py -n number_of_dimer_structure
```
## Testing environment

* Python==3.7.10
* Numpy==1.20.1
