MPICCXX = mpic++

sender.o: sender.cpp
	$(MPICCXX) -c $< -o $@

receiver.o: receiver.cpp
	$(MPICCXX) -c $< -o $@

test: test.cpp sender.o receiver.o
	$(MPICCXX) $^ -o $@
