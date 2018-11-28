all:
	mkdir -p build
	cd build; cmake -D CMAKE_BUILD_TYPE=release ..; make -j4; 
	cd build; ctest

debug:
	mkdir -p build
	cd build; cmake -D CMAKE_BUILD_TYPE=debug ..; make -j4;
	cd build; ctest

clean:
	rm -rf build
