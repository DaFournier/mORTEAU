# mORTEAU

David Fournier 2020

A tool to find 1:1 orthologues accross species. 
In the example, mouse orthologues in 13 other species
are tested. 

Orthologues are found in the following order:
- if proteome available, test protein sequence of mouse
against the proteome of the target species. 
- if no 1:1 orthologue found, a tBLASTn is performed 
against DNA genome.
- if no 1:1 orthologue found, a tBLASTn against
shotgun sequencing contigs is performed. 
- if no 1:1 orthologue found, absence of 1:1 orthologue. 

As a control, the decrease of orthology is much better
from mouse to lower forms if analysis goes further along
the 4 steps. 

