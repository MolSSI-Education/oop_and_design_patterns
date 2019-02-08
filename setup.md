---
title: Setup
---
We are utilizing two Molecular Dynamics packages to show examples of some of the design patterns.
[MDanalysis] and [MDTraj] are both available through conda-forge.

To start, open an Anaconda prompt.
Then enter `conda create --name design_patterns`.
Enter `Y` to continue.

We now need to activate our new environment.
~~~
	$ conda activate design_patterns
~~~
{: .bash}

Then we will install MDAnalysis and MDTraj.
~~~
	$ conda install -c conda-forge mdanalysis
~~~
{: .bash}
~~~
	$ conda install -c conda-forge mdtraj
~~~
{: .bash}
We have included a sample trajectory file to use with the example design patterns.
The [trajectory] was extracted from the [MDAnalysisTests]


{% include links.md %}
[trajectory]: ../data/protein.pdb
[MDAnalysisTests]: https://pypi.org/project/MDAnalysisTests/
[MDAnalysis]: https://www.mdanalysis.org/
[MDTraj]: http://mdtraj.org/1.9.0/
[NumPy]: http://www.numpy.org/