name = vbox

input: clean
	gmsh -2 src/$(name).geo -o src/$(name).msh
	cd src; ../../../scripts/gmsh2triangle --2D $(name).msh

clean:
	rm -rf *.d.* *.stat *.vtu $(name).node $(name).ele $(name).face *checkpoint* *.convergence
