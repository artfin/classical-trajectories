CCXX = g++

CXXFLAGS = -std=c++11 

OPENGL_LFLAGS = -lglut -lGL -lGLEW -lGLU 

test: test.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@

double_buffer: double_buffer.cpp
	$(CCXX) $(CXXFLAGS) $(OPENGL_LFLAGS) $^ -o $@
