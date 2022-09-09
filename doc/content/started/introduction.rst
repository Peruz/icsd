Introduction
------------

.. Should show a very short overview of what can be done with the product, using one or two extremely simplified use cases. This is the thirty-second pitch for your project.


Most of the codes available for the interpretation of geoelectrical survey focus on revovering the subsurface electrical resistivity. 
For some specific cases, when a direct excitation of the a conductive body (or mass) is applied, it is more relevant to map the electrical current density within the subsurface. 
The `pyGeoCSD` open-source generic algorithm provide a scheme for the inversion of subsurface current source density. 
It targets the geophysical community in primary and provide readers and parsers for two of the most widely used libraries for forward/inversion geoelectrical libraries i.e. pyGIMLI and Resipy. 
The `pyGeoCSD` package keep the dependencies to the bare minimum to encourage other libraries to depend on. 
The technique adopted is a logical extension of the work of [@binley_detecting_1997; @binley_detecting_1999]. 
The mathematical formulation is also borrowed from the neurosciences community who develop the same approach for the medical imaging.
