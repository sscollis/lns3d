# Collis PhD Thesis

README:		September 5, 1995

Finalized:	March 10, 1997

Author:		S. Scott Collis

This directory contains my Ph.D. thesis in latex format.

# Note

I have made some very minor changes to the Table of Contents and 
List of Figures to make these compatible with the hyperref LaTeX
package.   Now all references are hyperlinks making the PDF much
easier to navigate.

The utility script `clean` automatically removes all
LaTeX intermediate files.

# Building

```bash
latex thesis && bibtex thesis && latex thesis && latex thesis
dvipdf thesis
```
or use the script `build`.

To clean derived files
```bash
./clean
```
----
S. Scott Collis
flow.physics.simulation@gmail.com
