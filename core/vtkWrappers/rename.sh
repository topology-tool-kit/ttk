#!/bin/bash

for f in `ls -d vtk*`; do
	ttk_dir=`echo $f | sed s/vtk/ttk/`
	sed -ri "s/vtk([A-Z])/ttk\1/g" $f/CMakeLists.txt
	mv $f/CMakeLists.txt $ttk_dir/CMakeLists.txt
done

