CCXX = g++

CXXFLAGS = -std=c++11 

OPENGL_LFLAGS = -lglut -lGL -lGLEW -lGLU 

test: test.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

double_buffer: double_buffer.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

cube: cube.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

planet: planet.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

unproject: unproject.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

smooth_triangle: smooth_triangle.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

light_sphere: light_sphere.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

movelight: movelight.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

alpha: alpha.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

big_chunk: big_chunk.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@
