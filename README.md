# Internship Genova
## Max Royer

There is a fork of the CDT repo in the `code` folder.

## CDT - Constrained Delaunay Tetrahedrization made robust and practical
This code implements an algorithm to calculate a Constrained Delaunay Tetrahedrization (CDT) of an input PLC represented by on OFF file.
Steiner points are possibly added to make the input admit a CDT.
Details of the algorithm are described in "**Constrained Delaunay Tetrahedrization: A robust and practical approach**" by L. Diazzi, D. Panozzo, A. Vaxman and M. Attene (ACM Trans Graphics Vol 42, N. 6, Procs of SIGGRAPH Asia 2023). 
You may download a copy here: http://arxiv.org/abs/2309.09805

### Usage

Build the executable as follows:
```
cmake -B build -S .
```

This will produce an appropriate building configuration for your system.
On Windows MSVC, this will produce a cdt.sln file.
On Linux/MacOS, this will produce a Makefile. 
Use it as usual to compile cdt. Alternatively, you can use the command line:
```
cmake --build build --config Release
```

When compiled, the code generates an executable called ``cdt``.
Launch it with no command line parameters to have a list of supported options.

Example:

```
cdt input_file.off
```
creates a file called ``input_file.off.tet`` representing the constrained tetrahedrization.

