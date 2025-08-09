# Internship Genova
## Max Royer

This is the repository of my M1 internship at Genova CNR Imati in 2025 supervised by Marco Attene and Marco
Livesu. See the report in the `report` folder for further explanation of my internship.

## Usage

### Build

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

### Utilisation

When compiled, the code generates an executable called ``cdt``.
Launch it with no command line parameters to have a list of supported options.

Example:

```
cdt $input_file.off
```

creates a file called ``input_file.off.tet`` representing the constrained tetrahedrization. Some input files
are available in the Input_file folder.

**The option to switch the optimization on is `-o`**

### Tests

To perform tests you should run the `testing.py` file in the `tests` folder, unfortunately it was not
convenient to put the big dataset on the repository, so it will only run on the example models.

