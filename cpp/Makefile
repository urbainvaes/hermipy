all:
	mkdir -p build
	cd build; cmake -D CMAKE_BUILD_TYPE=release ..; make; 
	cd build; ctest

debug:
	mkdir -p build
	cd build; cmake -D CMAKE_BUILD_TYPE=debug ..; make
	cd build; ctest

clean:
	rm -rf build
