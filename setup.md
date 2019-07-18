---
title: Setup
---
## Installing Python through Anaconda
[Python](https://python.org/) is a popular language for scientific computing, and great for general-purpose programming as well. Installing all of its scientific packages individually can be a bit difficult, however, so we recommend the all-in-one installer Anaconda.

Regardless of how you choose to install it, *please make sure you install Python version 3.x (e.g., 3.4 is fine, 2.7 is not)*.  Also, please set up your python environment at least a day in advance of the workshop. If you encounter problems with the installation procedure, the instructors will be available 30 minutes before the workshop begins to help you.

## Windows - [Video Tutorial](https://www.youtube.com/watch?v=xxQ0mzZ8UvA)

1. Open the [Anaconda Windows download page](https://www.anaconda.com/download/#windows).
2. Download the installer.  **Be sure you get the Python 3 version.**
3. Double-click the installer icon and follow the setup instructions on screen.  You can use MOST of the default options.  The only exception is to check the **Make Anaconda the default Python** option.

## Mac OS X - [Video Tutorial](https://www.youtube.com/watch?v=TcSAln46u9U)

1. Open the [Anaconda MacOS download page](https://www.anaconda.com/download/#macos).
2. Download the installer. **Be sure you get the Python 3 version.**
3. Double click the installer icon and follow the setup instructions.  You can use all of the default options.

## Obtain lesson materials
We are utilizing two Molecular Dynamics packages to show examples of some of the design patterns.
[MDanalysis] and [MDTraj] are both available through conda-forge.

To start, open an Anaconda prompt.
Then enter `conda create --name design_patterns python=3.7`.
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
[trajectory]: ./data/protein.pdb
[MDAnalysisTests]: https://pypi.org/project/MDAnalysisTests/
[MDAnalysis]: https://www.mdanalysis.org/
[MDTraj]: http://mdtraj.org/1.9.0/
[NumPy]: http://www.numpy.org/