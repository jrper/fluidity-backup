name = vbox

input: clean
	gmsh -2 src/$(name).geo
	../../scripts/gmsh2triangle --2D src/$(name).msh

clean:
	rm -rf *.d.* *.stat *.vtu $(name).node $(name).ele $(name).face *checkpoint* *.convergence
