Introduction
------------

.. Should show a very short overview of what can be done with the product, using one or two extremely simplified use cases. This is the thirty-second pitch for your project.

Exemple on how to run ICSD code for a 3d prospection on a landfill.

The simpliest processing can be achieved with the python API::

    from icsd3d_class import iCSD3d_Class as i3d
    icsd3d_landfill=i3d(dirName=path2files)   
	icsd3d_landfill.type='3d' # Prospection dimension
	icsd3d_landfill.sim="SAno.txt" # Simulated green functions filename
	icsd3d_landfill.obs="OAno_synt.txt" # Observations file
	icsd3d_landfill.coord_file="VRTeCoord.txt" # Coordin
	icsd3d_landfill.regMesh='unstrc' # spatial regularisation VRTE mesh type

More examples are available in the Example gallery section.
